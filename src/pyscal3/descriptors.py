"""
Structural descriptors for ASE Atoms, backed by pyscal's C++ routines.

All functions take an ASE Atoms object (with neighbors already computed
via pyscal3.find_neighbors) and return computed values while also storing
results in atoms.arrays / atoms.info.

Example
-------
>>> from ase.build import bulk
>>> import pyscal3
>>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
>>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)
>>> q = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
>>> print(atoms.arrays["pyscal_q4"])
"""
import numpy as np
import itertools
from ase import Atoms

import pyscal3.csystem as pc
from pyscal3._bridge import (
    get_box_params, atoms_to_dict, dict_to_atoms,
    ensure_neighbors, create_attribute,
    pad_atoms_for_neighbor_finding,
)
from pyscal3.neighbors import find_neighbors


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _get_dict_with_neighbors(atoms: Atoms) -> dict:
    """Build the C++ dict, ensuring neighbors exist."""
    ensure_neighbors(atoms)
    return atoms_to_dict(atoms)


def _sync_back(d: dict, atoms: Atoms, keys: list):
    """Write specific keys from dict back into atoms."""
    n = len(atoms)
    for key in keys:
        if key not in d:
            continue
        store_key = "pyscal_" + key
        try:
            arr = np.asarray(d[key])
            if arr.ndim >= 1 and len(arr) == n and arr.dtype.kind != 'O':
                atoms.arrays[store_key] = arr
                continue
        except (ValueError, TypeError):
            pass
        atoms.info[store_key] = d[key]


# ---------------------------------------------------------------------------
# Steinhardt Parameters
# ---------------------------------------------------------------------------

def steinhardt_parameter(atoms: Atoms, l, averaged=False):
    """
    Calculate Steinhardt bond order parameters q_l.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    l : int or list of int
        Steinhardt parameter order(s), e.g. 6 or [4, 6].
    averaged : bool, optional
        If True, compute neighbor-averaged q values. Default False.

    Returns
    -------
    list of numpy arrays
        One array per requested l value, each of shape (natoms,).
    """
    if isinstance(l, int):
        ll = [l]
    else:
        ll = list(l)

    d = _get_dict_with_neighbors(atoms)

    if averaged:
        # Need base q values first
        for val in ll:
            pc.calculate_q_single(d, val)
        for val in ll:
            pc.calculate_aq_single(d, val)
        result_keys = ["avg_q%d" % v for v in ll]
    else:
        for val in ll:
            pc.calculate_q_single(d, val)
        result_keys = ["q%d" % v for v in ll]

    # Sync all q-related keys back
    sync_keys = []
    for val in ll:
        sync_keys.extend([
            "q%d" % val, "q%d_real" % val, "q%d_imag" % val,
            "avg_q%d" % val,
        ])
    _sync_back(d, atoms, sync_keys)

    return [np.array(d[k]) for k in result_keys]


# ---------------------------------------------------------------------------
# Moment Tensor Potentials (MTP) Descriptors
# ---------------------------------------------------------------------------

def _mtp_cutoff(r, r_cut):
    """Smooth polynomial cutoff for MTP: (1 - r/r_cut)^2."""
    mask = r < r_cut
    result = np.zeros_like(r)
    result[mask] = (1.0 - r[mask] / r_cut)**2
    return result


def _chebyshev_basis(r, n_basis, r_min, r_max):
    """
    Evaluate Chebyshev polynomials T_0 to T_{n_basis-1}.
    
    Maps r to [-1, 1] via: x = (2r - r_min - r_max) / (r_max - r_min)
    """
    # Clip r to [r_min, r_max] for stability
    r_clipped = np.clip(r, r_min, r_max)
    
    # Map to [-1, 1]
    x = (2*r_clipped - r_min - r_max) / (r_max - r_min + 1e-10)
    
    basis = np.zeros((len(r), n_basis))
    basis[:, 0] = 1.0
    if n_basis > 1:
        basis[:, 1] = x
    for n in range(2, n_basis):
        basis[:, n] = 2*x*basis[:, n-1] - basis[:, n-2]
    
    return basis


def moment_tensor_descriptors(atoms: Atoms, level=8, n_radial=4, r_cut=5.0,
                               normalized=True):
    """
    Calculate MTP (Moment Tensor Potentials) descriptors.
    
    MTP descriptors provide systematically improvable invariant descriptions
    of local atomic environments through moment tensor contractions.
    
    Parameters
    ----------
    atoms : Atoms
        Structure with neighbors already computed via find_neighbors.
    level : int, optional
        Maximum level for basis functions. Higher = more descriptors.
        Typical values: 8, 10, 12. Default 8.
    n_radial : int, optional
        Number of radial (Chebyshev) basis functions. Default 4.
    r_cut : float, optional
        Cutoff radius. Default 5.0.
    normalized : bool, optional
        If True, L2-normalize descriptors. Default True.
    
    Returns
    -------
    dict
        "descriptors": (N, D) array of MTP descriptors
        "level": level used
        "n_radial": number of radial functions
        "descriptor_size": D
    
    Notes
    -----
    The moment tensor of rank nu is:
    
        M_{mu,nu} = sum_j f_mu(r_ij) * r_ij ⊗ r_ij ⊗ ... ⊗ r_ij (nu times)
    
    Scalar invariants are obtained by contracting these tensors.
    
    Level controls which basis functions are included:
    - Level 2: M_{0,0} (coordination-like)
    - Level 6: M_{1,0}, tr(M_{0,2})
    - Level 8: M_{0,1} · M_{0,1} (angular correlations)
    - Higher: More body-order terms
    
    Results stored in atoms.arrays["pyscal_mtp"].
    
    References
    ----------
    Shapeev, Multiscale Model. Simul. 14, 1153-1173 (2016)
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
    >>> result = pyscal3.moment_tensor_descriptors(atoms, level=8)
    >>> print(result["descriptors"].shape)
    """
    d = _get_dict_with_neighbors(atoms)
    n_atoms = len(atoms)
    
    # Collect all distances to determine r_min
    all_dists = []
    for i in range(n_atoms):
        neighbors_i = d["neighbors"][i]
        if hasattr(neighbors_i, '__len__') and len(neighbors_i) > 0:
            dists_i = np.array(d["neighbordist"][i])
            mask = dists_i < r_cut
            if np.any(mask):
                all_dists.extend(dists_i[mask].tolist())
    
    r_min = min(all_dists) * 0.5 if all_dists else 0.5
    r_max = r_cut
    
    # Determine which basis functions to include based on level
    # Format: list of (function_name, mu_list, nu_list)
    # Level formula: lev = 2 + 4*mu + nu for single tensor
    basis_specs = []
    
    # Level 2: M_{0,0}
    if level >= 2:
        for mu in range(n_radial):
            if 2 + 4*mu <= level:
                basis_specs.append(("M", mu, 0))
    
    # Level 6+: tr(M_{mu,2}) 
    if level >= 6:
        for mu in range(n_radial):
            if 2 + 4*mu + 2 <= level:
                basis_specs.append(("trM2", mu, 2))
    
    # Level 8+: M_{mu1,1} · M_{mu2,1} (dot product of rank-1 tensors)
    if level >= 8:
        for mu1 in range(n_radial):
            for mu2 in range(mu1, n_radial):  # mu2 >= mu1 by symmetry
                if 2 + 4*mu1 + 1 + 4*mu2 + 1 <= level:
                    basis_specs.append(("M1M1", mu1, mu2))
    
    # Level 10+: tr(tr(M_{mu,4})) - trace twice on rank-4 tensor
    if level >= 10:
        for mu in range(n_radial):
            if 2 + 4*mu + 4 <= level:
                basis_specs.append(("trtrM4", mu, 4))
    
    desc_size = len(basis_specs)
    if desc_size == 0:
        desc_size = 1  # Ensure at least one descriptor
        basis_specs.append(("M", 0, 0))
    
    descriptors = np.zeros((n_atoms, desc_size))
    
    for i in range(n_atoms):
        neighbors_i = d["neighbors"][i]
        n_neigh = len(neighbors_i) if hasattr(neighbors_i, '__len__') else 0
        
        if n_neigh == 0:
            continue
        
        # Get neighbor displacements and distances
        diffs = np.array(d["diff"][i])  # (n_neigh, 3)
        dists = np.array(d["neighbordist"][i])  # (n_neigh,)
        
        # Filter to within r_cut
        mask = dists < r_cut
        if not np.any(mask):
            continue
        
        diffs = diffs[mask]
        dists = dists[mask]
        n_neigh = len(dists)
        
        # Cutoff weights
        f_cut = _mtp_cutoff(dists, r_cut)
        
        # Radial basis values: (n_neigh, n_radial)
        phi = _chebyshev_basis(dists, n_radial, r_min, r_max)
        
        # Weighted radial basis
        f_mu = phi * f_cut[:, np.newaxis]  # (n_neigh, n_radial)
        
        # Compute each basis function
        for idx, spec in enumerate(basis_specs):
            if spec[0] == "M":
                # M_{mu,0} = sum_j f_mu(r_ij) - just radial sum
                mu = spec[1]
                descriptors[i, idx] = np.sum(f_mu[:, mu])
            
            elif spec[0] == "trM2":
                # tr(M_{mu,2}) = sum_j f_mu * |r_ij|^2
                mu = spec[1]
                r2 = dists**2
                descriptors[i, idx] = np.sum(f_mu[:, mu] * r2)
            
            elif spec[0] == "M1M1":
                # M_{mu1,1} · M_{mu2,1} = sum_{j,k} f_mu1(j) f_mu2(k) (r_j · r_k)
                mu1, mu2 = spec[1], spec[2]
                # M_{mu,1} = sum_j f_mu(j) * r_j  -> (3,) vector
                M1 = np.sum(f_mu[:, mu1, np.newaxis] * diffs, axis=0)
                M2 = np.sum(f_mu[:, mu2, np.newaxis] * diffs, axis=0)
                descriptors[i, idx] = np.dot(M1, M2)
            
            elif spec[0] == "trtrM4":
                # Double trace of rank-4 tensor: sum_j f_mu |r_j|^4
                # This is trace(trace(M_{mu,4})) which simplifies to scalar
                mu = spec[1]
                r4 = dists**4
                descriptors[i, idx] = np.sum(f_mu[:, mu] * r4)
    
    if normalized:
        norms = np.linalg.norm(descriptors, axis=1, keepdims=True)
        norms[norms < 1e-15] = 1.0
        descriptors = descriptors / norms
    
    # Store results
    atoms.arrays["pyscal_mtp"] = descriptors
    
    return {
        "descriptors": descriptors,
        "level": level,
        "n_radial": n_radial,
        "r_cut": r_cut,
        "descriptor_size": desc_size,
        "basis_specs": basis_specs,
    }


# ---------------------------------------------------------------------------
# Disorder Parameter
# ---------------------------------------------------------------------------

def disorder(atoms: Atoms, q=6, averaged=False):
    """
    Calculate the disorder parameter.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    q : int, optional
        Steinhardt parameter order for disorder calc. Default 6.
    averaged : bool, optional
        If True, average disorder over neighbors. Default False.

    Returns
    -------
    numpy array
        Per-atom disorder values.
    """
    d = _get_dict_with_neighbors(atoms)

    # Ensure q values exist
    keys_needed = ["q%d_real" % q, "q%d_imag" % q]
    need_calc = False
    for k in keys_needed:
        if k not in d:
            need_calc = True
            break
        v = d[k]
        if not hasattr(v, '__len__') or len(v) == 0:
            need_calc = True
            break
        # Check if it looks like a 2-D container (list-of-lists or 2-D array)
        first = v[0]
        if not (hasattr(first, '__len__') and not isinstance(first, str)):
            need_calc = True
            break
    if need_calc:
        pc.calculate_q_single(d, q)

    pc.calculate_disorder(d, q)

    sync_keys = ["disorder", "q%d" % q, "q%d_real" % q, "q%d_imag" % q]

    if averaged:
        # Average disorder over neighbors in C++
        pc.calculate_average_disorder(d)
        sync_keys.append("avg_disorder")

    _sync_back(d, atoms, sync_keys)

    if averaged:
        return np.array(d["avg_disorder"])
    return np.array(d["disorder"])


# ---------------------------------------------------------------------------
# Common Neighbor Analysis
# ---------------------------------------------------------------------------

def common_neighbor_analysis(atoms: Atoms, lattice_constant=None):
    """
    Calculate Common Neighbor Analysis (CNA) or Adaptive CNA.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure.
    lattice_constant : float, optional
        If given, use conventional CNA. If None, use adaptive CNA.

    Returns
    -------
    dict
        Counts: {"fcc": n, "hcp": n, "bcc": n, "ico": n, "others": n}
    """
    d, (triclinic, rot, rotinv, boxdims), nreal = pad_atoms_for_neighbor_finding(atoms)
    n = len(d["positions"])

    # Create structure attribute
    d["structure"] = [0] * n

    # Find temp neighbors (by number, nmax=14)
    _reset_and_find_temp_neighbors(d, triclinic, rot, rotinv, boxdims, nmax=14)

    if lattice_constant is None:
        # Adaptive CNA
        pc.get_acna_neighbors_cn12(d, triclinic, rot, rotinv, boxdims)
        pc.identify_cn12(d, triclinic, rot, rotinv, boxdims)
        pc.get_acna_neighbors_cn14(d, triclinic, rot, rotinv, boxdims)
        pc.identify_cn14(d, triclinic, rot, rotinv, boxdims)
    else:
        pc.get_cna_neighbors(d, triclinic, rot, rotinv, boxdims, lattice_constant, 1)
        pc.identify_cn12(d, triclinic, rot, rotinv, boxdims)
        pc.get_cna_neighbors(d, triclinic, rot, rotinv, boxdims, lattice_constant, 2)
        pc.identify_cn14(d, triclinic, rot, rotinv, boxdims)

    structure = np.array(d["structure"][:nreal])
    atoms.arrays["pyscal_structure"] = structure

    return {
        "others": int(np.sum(structure == 0)),
        "fcc": int(np.sum(structure == 1)),
        "hcp": int(np.sum(structure == 2)),
        "bcc": int(np.sum(structure == 3)),
        "ico": int(np.sum(structure == 4)),
    }


def diamond_structure(atoms: Atoms):
    """
    Identify diamond structure using extended CNA.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure.

    Returns
    -------
    dict
        Counts per structure type.
    """
    d, (triclinic, rot, rotinv, boxdims), nreal = pad_atoms_for_neighbor_finding(atoms)
    n = len(d["positions"])

    d["structure"] = [0] * n
    _reset_and_find_temp_neighbors(d, triclinic, rot, rotinv, boxdims, nmax=4)

    pc.identify_diamond_cna(d, triclinic, rot, rotinv, boxdims)

    structure = np.array(d["structure"][:nreal])
    atoms.arrays["pyscal_structure"] = structure

    return {
        "others": int(np.sum(structure == 0)),
        "cubic diamond": int(np.sum(structure == 1)),
        "cubic diamond 1NN": int(np.sum(structure == 2)),
        "cubic diamond 2NN": int(np.sum(structure == 3)),
        "hex diamond": int(np.sum(structure == 4)),
        "hex diamond 1NN": int(np.sum(structure == 5)),
        "hex diamond 2NN": int(np.sum(structure == 6)),
    }


# ---------------------------------------------------------------------------
# Centrosymmetry Parameter
# ---------------------------------------------------------------------------

def centrosymmetry(atoms: Atoms, nmax=12):
    """
    Calculate the centrosymmetry parameter.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure.
    nmax : int, optional
        Number of neighbors (must be positive even integer). Default 12.

    Returns
    -------
    numpy array
        Per-atom centrosymmetry values.
    """
    if nmax <= 0:
        raise ValueError("nmax must be positive")
    if nmax % 2 != 0:
        raise ValueError("nmax must be even")

    # Find neighbors by number
    find_neighbors(atoms, method='number', nmax=nmax, assign_neighbor=True)
    d = _get_dict_with_neighbors(atoms)
    d["centrosymmetry"] = [0.0] * len(atoms)

    pc.calculate_centrosymmetry(d, nmax)

    cs = np.array(d["centrosymmetry"])
    atoms.arrays["pyscal_centrosymmetry"] = cs
    return cs


# ---------------------------------------------------------------------------
# Voronoi Vector
# ---------------------------------------------------------------------------

def voronoi_vector(atoms: Atoms, edge_cutoff=0.05, area_cutoff=0.01):
    """
    Calculate the Voronoi structure identification vector (n3, n4, n5, n6).

    Parameters
    ----------
    atoms : ase.Atoms
        Structure (must have Voronoi neighbors computed).
    edge_cutoff : float, optional
        Minimum edge length fraction. Default 0.05.
    area_cutoff : float, optional
        Minimum face area fraction. Default 0.01.

    Returns
    -------
    numpy array of shape (natoms, 4)
        Voronoi vectors [n3, n4, n5, n6] per atom.
    """
    d = _get_dict_with_neighbors(atoms)

    if "face_vertices" not in d:
        raise ValueError(
            "Voronoi analysis required. Call find_neighbors(atoms, method='voronoi') first."
        )

    pc.calculate_voronoi_vector(d, edge_cutoff, area_cutoff)

    vv = np.array(d["vorovector"])
    atoms.arrays["pyscal_vorovector"] = vv
    return vv


# ---------------------------------------------------------------------------
# Entropy Parameter
# ---------------------------------------------------------------------------

def entropy(atoms: Atoms, rm, sigma=0.2, rstart=0.001, h=0.001,
            local=False, average=False):
    """
    Calculate the entropy parameter for each atom.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    rm : float
        Cutoff distance for integration.
    sigma : float, optional
        Broadening parameter. Default 0.2.
    rstart : float, optional
        Integration start. Default 0.001.
    h : float, optional
        Integration step (trapezoidal). Default 0.001.
    local : bool, optional
        Use local density instead of global. Default False.
    average : bool, optional
        Compute neighbor-averaged entropy. Default False.

    Returns
    -------
    numpy array
        Per-atom entropy (or averaged entropy) values.
    """
    d = _get_dict_with_neighbors(atoms)

    n = len(atoms)
    volume = abs(np.linalg.det(atoms.cell))
    kb = 1

    if local:
        rho = 0
    else:
        rho = n / volume

    pc.calculate_entropy(d, sigma, rho, rstart, rm, h, kb)

    sync_keys = ["entropy"]

    if average:
        pc.calculate_average_entropy(d)
        sync_keys.append("average_entropy")

    _sync_back(d, atoms, sync_keys)

    if average:
        return np.array(d["average_entropy"])
    return np.array(d["entropy"])


# ---------------------------------------------------------------------------
# Short-Range Order
# ---------------------------------------------------------------------------

def short_range_order(atoms: Atoms, reference_type=1, compare_type=2, average=True):
    """
    Calculate Warren-Cowley short-range order parameter.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    reference_type : int, optional
        Atomic number of reference type. Default 1.
    compare_type : int, optional
        Atomic number to compare. Default 2.
    average : bool, optional
        If True, return system average. Default True.

    Returns
    -------
    numpy array or float
        Per-atom SRO values, or system average.
    """
    d = _get_dict_with_neighbors(atoms)

    pc.calculate_short_range_order(d, reference_type, compare_type)

    sro = np.array(d["sro"])
    atoms.arrays["pyscal_sro"] = sro

    if average:
        return float(np.mean(sro))
    return sro


# ---------------------------------------------------------------------------
# Radial Distribution Function
# ---------------------------------------------------------------------------

def radial_distribution_function(atoms: Atoms, rmin=0, rmax=5.0, bins=100):
    """
    Calculate radial distribution function g(r).

    Parameters
    ----------
    atoms : ase.Atoms
        Structure.
    rmin, rmax : float
        Distance range.
    bins : int
        Number of histogram bins.

    Returns
    -------
    (rdf, r) : tuple of numpy arrays
    """
    find_neighbors(atoms, method="cutoff", cutoff=rmax)
    d = atoms_to_dict(atoms)

    distances = np.concatenate(d["neighbordist"])
    hist, bin_edges = np.histogram(distances, bins=bins, range=(rmin, rmax), density=True)

    edgewidth = abs(bin_edges[1] - bin_edges[0])
    r = bin_edges[:-1]
    n = len(atoms)
    volume = abs(np.linalg.det(atoms.cell))
    rho = n / volume

    shell_vols = (4.0 / 3.0) * np.pi * ((r + edgewidth) ** 3 - r ** 3)
    rdf = (hist / shell_vols) / rho

    return rdf, r


# ---------------------------------------------------------------------------
# Angular Criteria
# ---------------------------------------------------------------------------

def angular_criteria(atoms: Atoms):
    """
    Calculate angular criteria for diamond structure identification.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.

    Returns
    -------
    numpy array
        Per-atom angular parameter A values.
    """
    d = _get_dict_with_neighbors(atoms)

    pc.calculate_angular_criteria(d)

    ang = np.array(d["angular"])
    atoms.arrays["pyscal_angular"] = ang
    return ang


# ---------------------------------------------------------------------------
# Chi Parameters
# ---------------------------------------------------------------------------

def chi_params(atoms: Atoms, angles=False):
    """
    Calculate chi-parameter vector for structure identification.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    angles : bool, optional
        If True, also return cosine angles. Default False.

    Returns
    -------
    numpy array of shape (natoms, 9)
        Chi parameter vectors.
    """
    d = _get_dict_with_neighbors(atoms)

    pc.calculate_chi_params(d)

    cp = np.array(d["chiparams"])
    atoms.arrays["pyscal_chiparams"] = cp

    if angles:
        cosines_list = d["cosines"]
        atoms.info["pyscal_cosines"] = cosines_list
        return cp, cosines_list
    return cp


# ---------------------------------------------------------------------------
# Solid/Liquid Identification
# ---------------------------------------------------------------------------

def find_solids(atoms: Atoms, bonds=0.5, threshold=0.5, avgthreshold=0.6,
                cluster=True, q=6, cutoff=0, right=True):
    """
    Distinguish solid and liquid atoms.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    bonds : int or float
        Min solid bonds (int) or min fraction of solid neighbors (float 0-1).
    threshold : float
        Bond correlation cutoff. Default 0.5.
    avgthreshold : float
        Average bond cutoff. Default 0.6.
    cluster : bool
        If True, cluster solid atoms and return largest cluster size.
    q : int
        Steinhardt parameter order. Default 6.
    cutoff : float
        Cluster cutoff (0 = use neighbor cutoff).
    right : bool
        If True, use > comparison. Default True.

    Returns
    -------
    int or None
        Largest cluster size if cluster=True.
    """
    d = _get_dict_with_neighbors(atoms)

    if isinstance(bonds, int):
        criteria = 0
    elif isinstance(bonds, float) and 0 <= bonds <= 1:
        criteria = 1
    else:
        raise TypeError("bonds must be int or float in [0,1]")

    compare_criteria = 0 if right else 1

    # Calculate Steinhardt parameters
    pc.calculate_q_single(d, q)

    # Calculate bonds/solid classification
    pc.calculate_bonds(d, q, threshold, avgthreshold, bonds, compare_criteria, criteria)

    _sync_back(d, atoms, ["solid", "bonds", "sij", "avg_sij",
                           "q%d" % q, "q%d_real" % q, "q%d_imag" % q])

    if cluster:
        return find_clusters(atoms, condition=np.array(d["solid"]) > 0,
                            cutoff=cutoff, d=d)
    return None


def find_clusters(atoms: Atoms, condition, largest=True, cutoff=0, d=None):
    """
    Cluster atoms based on a boolean condition.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    condition : array-like of bool
        Which atoms to include in clustering.
    largest : bool
        If True, return largest cluster size.
    cutoff : float
        Cluster cutoff (0 = use neighbor cutoff).

    Returns
    -------
    int
        Largest cluster size if largest=True.
    """
    if d is None:
        d = _get_dict_with_neighbors(atoms)

    condition = np.asarray(condition, dtype=bool)
    d["condition"] = condition.tolist()

    pc.find_clusters(d, cutoff)

    cluster_ids = np.array(d["cluster"])
    atoms.arrays["pyscal_cluster"] = cluster_ids

    if largest:
        valid = cluster_ids[cluster_ids >= 0]
        if len(valid) > 0:
            unique, counts = np.unique(valid, return_counts=True)
            largest_size = int(counts.max())
            largest_id = unique[counts.argmax()]
            atoms.arrays["pyscal_largest_cluster"] = (cluster_ids == largest_id)
            return largest_size
        return 0
    return None


# ---------------------------------------------------------------------------
# Average over neighbors (utility)
# ---------------------------------------------------------------------------

def average_over_neighbors(atoms: Atoms, key: str, include_self=True):
    """
    Average a per-atom property over each atom's neighbors.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    key : str
        Key in atoms.arrays (with or without 'pyscal_' prefix).
    include_self : bool
        Include the atom itself in the average. Default True.

    Returns
    -------
    numpy array
    """
    d = _get_dict_with_neighbors(atoms)

    # Find the data
    lookup = key
    if key not in d and "pyscal_" + key in atoms.arrays:
        lookup = key
        values = atoms.arrays["pyscal_" + key]
    elif key in d:
        values = np.array(d[key])
    else:
        raise KeyError(f"Property '{key}' not found")

    # 1-D values: use fast C++ averaging
    values = np.asarray(values)
    if values.ndim == 1:
        result = pc.calculate_average_over_neighbors(
            d, values.tolist(), include_self)
        return np.array(result)

    # Multi-dimensional: fall back to Python loop
    result = []
    for i in range(len(atoms)):
        vals = [values[i]] if include_self else []
        for j in d["neighbors"][i]:
            vals.append(values[j])
        result.append(np.mean(vals))

    return np.array(result)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _reset_and_find_temp_neighbors(d, triclinic, rot, rotinv, boxdims, nmax=14):
    """Reset neighbors and find by number (for CNA/diamond)."""
    n = len(d["positions"])
    d["neighbors"] = [[] for _ in range(n)]
    d["neighbordist"] = [[] for _ in range(n)]
    d["temp_neighbors"] = [[] for _ in range(n)]
    d["temp_neighbordist"] = [[] for _ in range(n)]
    d["neighborweight"] = [[] for _ in range(n)]
    d["diff"] = [[] for _ in range(n)]
    d["r"] = [[] for _ in range(n)]
    d["theta"] = [[] for _ in range(n)]
    d["phi"] = [[] for _ in range(n)]
    d["cutoff"] = [0.0] * n

    pc.get_all_neighbors_bynumber(
        d, 0.0, triclinic, rot, rotinv, boxdims,
        2, nmax, (n > 250), False)
