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
# SOAP (Smooth Overlap of Atomic Positions)
# ---------------------------------------------------------------------------

def _soap_cutoff(r, r_cut):
    """Smooth cutoff function for SOAP: 0.5*(cos(pi*r/r_cut) + 1) for r < r_cut."""
    mask = r < r_cut
    result = np.zeros_like(r)
    result[mask] = 0.5 * (np.cos(np.pi * r[mask] / r_cut) + 1.0)
    return result


def _soap_radial_basis(r, n_max, r_cut, sigma):
    """
    Compute orthonormalized radial basis functions evaluated at distances r.
    
    Uses polynomial basis phi_n(r) = r^n * exp(-alpha*r^2), orthonormalized
    via Cholesky decomposition.
    
    Parameters
    ----------
    r : array
        Distances to evaluate at
    n_max : int
        Number of radial basis functions
    r_cut : float
        Cutoff radius
    sigma : float
        Gaussian width
        
    Returns
    -------
    array of shape (len(r), n_max)
        Radial basis values
    """
    alpha = 1.0 / (2 * sigma**2)
    
    # Build primitive basis and overlap matrix on a fine grid
    n_grid = 500
    r_grid = np.linspace(0, r_cut, n_grid)
    dr = r_grid[1] - r_grid[0]
    
    # Primitive basis on grid: phi_n(r) = r^n * exp(-alpha*r^2)
    phi_grid = np.zeros((n_grid, n_max))
    for n in range(n_max):
        phi_grid[:, n] = (r_grid ** n) * np.exp(-alpha * r_grid**2)
    
    # Overlap matrix S_nn' = integral of phi_n * phi_n' * r^2 dr
    # Numerical integration
    weight = r_grid**2 * dr
    S = phi_grid.T @ (phi_grid * weight[:, np.newaxis])
    
    # Regularize if needed
    S += 1e-10 * np.eye(n_max)
    
    # Orthonormalization via Cholesky: S = L L^T, W = inv(L)
    try:
        L = np.linalg.cholesky(S)
        W = np.linalg.inv(L)
    except np.linalg.LinAlgError:
        # Fallback to SVD-based orthonormalization
        U, s, Vt = np.linalg.svd(S)
        W = U @ np.diag(1.0 / np.sqrt(s + 1e-10)) @ Vt
    
    # Evaluate primitive basis at actual r values
    phi_r = np.zeros((len(r), n_max))
    for n in range(n_max):
        phi_r[:, n] = (r ** n) * np.exp(-alpha * r**2)
    
    # Apply orthonormalization: g_n = sum_n' W_nn' * phi_n'
    g_r = phi_r @ W.T
    
    return g_r


def soap(atoms: Atoms, r_cut=5.0, n_max=8, l_max=6, sigma=0.5, normalized=True):
    """
    Calculate SOAP (Smooth Overlap of Atomic Positions) descriptors.
    
    SOAP provides a rotationally invariant representation of the local atomic
    environment, widely used in machine learning potentials.
    
    Parameters
    ----------
    atoms : Atoms
        Structure with neighbors already computed via find_neighbors.
        The cutoff used for find_neighbors should be >= r_cut.
    r_cut : float, optional
        Cutoff radius in Angstrom. Default 5.0.
    n_max : int, optional
        Number of radial basis functions. Default 8.
    l_max : int, optional
        Maximum angular momentum. Default 6.
    sigma : float, optional
        Gaussian width for smearing. Default 0.5.
    normalized : bool, optional
        If True, L2-normalize the power spectrum. Default True.
    
    Returns
    -------
    dict
        "power_spectrum": (N, D) array of SOAP power spectra
        "n_max": n_max used
        "l_max": l_max used
        "descriptor_size": D
        
    Notes
    -----
    The density of each environment is expanded as:
    
        rho(r) = sum_j f_cut(r_ij) * exp(-|r - r_ij|^2 / 2*sigma^2)
    
    Expansion coefficients c_nlm are computed, and the rotationally invariant
    power spectrum is:
    
        p_nn'l = sqrt(8/(2l+1)) * sum_m c*_nlm c_n'lm
    
    Results stored in atoms.arrays["pyscal_soap"].
    
    References
    ----------
    Bartók, Kondor, Csányi, Phys. Rev. B 87, 184115 (2013)
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
    >>> result = pyscal3.soap(atoms, r_cut=5.0, n_max=6, l_max=4)
    >>> print(result["power_spectrum"].shape)
    """
    from scipy.special import sph_harm_y
    
    d = _get_dict_with_neighbors(atoms)
    n_atoms = len(atoms)
    positions = np.array(d["positions"])
    
    # Descriptor dimension: n_max*(n_max+1)/2 * (l_max+1)
    # (using symmetry p_nn'l = p_n'nl)
    n_pairs = n_max * (n_max + 1) // 2
    desc_size = n_pairs * (l_max + 1)
    
    power_spectrum = np.zeros((n_atoms, desc_size))
    
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
        f_cut = _soap_cutoff(dists, r_cut)
        
        # Spherical coordinates
        # theta = arccos(z/r), phi = arctan2(y, x)
        theta = np.arccos(np.clip(diffs[:, 2] / dists, -1.0, 1.0))
        phi = np.arctan2(diffs[:, 1], diffs[:, 0])
        
        # Radial basis values: (n_neigh, n_max)
        g_nj = _soap_radial_basis(dists, n_max, r_cut, sigma)
        
        # Weight by cutoff
        g_nj = g_nj * f_cut[:, np.newaxis]
        
        # Compute c_nlm coefficients by summing over neighbors
        # c_nlm = sum_j g_n(r_j) * Y_l^m(theta_j, phi_j)*
        c_nlm = np.zeros((n_max, l_max + 1, 2 * l_max + 1), dtype=complex)
        
        for l in range(l_max + 1):
            for m in range(-l, l + 1):
                # scipy.special.sph_harm_y uses (n, m, theta, phi) where n=l
                ylm = sph_harm_y(l, m, theta, phi)  # (n_neigh,)
                for n in range(n_max):
                    c_nlm[n, l, m + l_max] = np.sum(g_nj[:, n] * np.conj(ylm))
        
        # Compute power spectrum p_nn'l
        idx = 0
        for l in range(l_max + 1):
            for n in range(n_max):
                for np_ in range(n, n_max):  # n' >= n (exploit symmetry)
                    # p_nn'l = sqrt(8/(2l+1)) * sum_m c*_nlm c_n'lm
                    val = 0.0
                    for m in range(-l, l + 1):
                        val += np.conj(c_nlm[n, l, m + l_max]) * c_nlm[np_, l, m + l_max]
                    power_spectrum[i, idx] = np.real(val) * np.sqrt(8.0 / (2 * l + 1))
                    idx += 1
    
    if normalized:
        # L2 normalize each atom's descriptor
        norms = np.linalg.norm(power_spectrum, axis=1, keepdims=True)
        norms[norms < 1e-15] = 1.0
        power_spectrum = power_spectrum / norms
    
    # Store results
    atoms.arrays["pyscal_soap"] = power_spectrum
    
    return {
        "power_spectrum": power_spectrum,
        "n_max": n_max,
        "l_max": l_max,
        "r_cut": r_cut,
        "sigma": sigma,
        "descriptor_size": desc_size,
    }


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
