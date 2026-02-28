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
# Behler-Parrinello Symmetry Functions (ACSFs)
# ---------------------------------------------------------------------------

# Default parameter sets for common systems
ACSF_PARAMS_SIMPLE = {
    "cutoff": 6.0,
    "G2": [(0.001, 0.0), (0.01, 0.0), (0.03, 0.0), (0.06, 0.0),
           (0.01, 1.5), (0.01, 2.5), (0.01, 3.5), (0.01, 4.5)],
    "G4": [(0.001, 1, 1), (0.001, 1, -1), (0.001, 2, 1), (0.001, 2, -1),
           (0.01, 1, 1), (0.01, 1, -1), (0.01, 4, 1), (0.01, 4, -1)],
}


def _cutoff_cosine(r: np.ndarray, rc: float) -> np.ndarray:
    """
    Cosine cutoff function.
    
    f_c(r) = 0.5 * (cos(π r / r_c) + 1) for r < r_c, else 0
    """
    result = np.zeros_like(r)
    mask = r < rc
    result[mask] = 0.5 * (np.cos(np.pi * r[mask] / rc) + 1.0)
    return result


def symmetry_functions(
    atoms: Atoms,
    cutoff: float = 6.0,
    g2_params: list = None,
    g4_params: list = None,
    g5_params: list = None,
    element_filter: list = None,
) -> dict:
    """
    Compute Behler-Parrinello atom-centered symmetry functions (ACSFs).
    
    These are smooth, many-body descriptors useful as input features for
    machine learning potentials. This is a pedagogical implementation;
    for production use, consider specialized packages like DScribe.
    
    Parameters
    ----------
    atoms : Atoms
        ASE Atoms object (neighbors must be found via find_neighbors first
        with cutoff >= the ACSF cutoff).
    cutoff : float, optional
        Cutoff radius in Angstroms (default: 6.0).
    g2_params : list of (eta, Rs) tuples, optional
        Parameters for G2 radial functions. Default uses 8 standard values.
        - eta: Gaussian width (Å⁻²)
        - Rs: Gaussian center (Å)
    g4_params : list of (eta, zeta, lambda) tuples, optional
        Parameters for G4 angular functions. Default uses 8 standard values.
        - eta: radial decay (Å⁻²)
        - zeta: angular resolution exponent
        - lambda: +1 or -1 (preferred angle direction)
    g5_params : list of (eta, zeta, lambda) tuples, optional
        Parameters for G5 angular functions (like G4 but without j-k distance).
        Default: None (not computed).
    element_filter : list of int, optional
        Atomic numbers to include. Default: all unique elements in atoms.
    
    Returns
    -------
    dict with keys:
        "G2": (N, n_g2 * n_elements) array of G2 values
        "G4": (N, n_g4 * n_element_pairs) array of G4 values (if g4_params)
        "G5": (N, n_g5 * n_element_pairs) array of G5 values (if g5_params)
        "features": concatenation of all features
        "element_types": list of atomic numbers used
        "feature_names": list of feature descriptions
    
    Notes
    -----
    Results also stored in:
    - atoms.arrays["pyscal_acsf"]: concatenated feature vector
    - atoms.info["pyscal_acsf_params"]: parameter metadata
    
    References
    ----------
    Behler & Parrinello, PRL 98, 146401 (2007)
    Behler, J. Chem. Phys. 134, 074106 (2011)
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
    >>> sf = pyscal3.symmetry_functions(atoms, cutoff=6.0)
    >>> print(sf["features"].shape)  # (108, n_features)
    """
    # Default parameters
    if g2_params is None:
        g2_params = ACSF_PARAMS_SIMPLE["G2"]
    if g4_params is None:
        g4_params = ACSF_PARAMS_SIMPLE["G4"]
    
    d = _get_dict_with_neighbors(atoms)
    n = len(atoms)
    positions = np.array(d["positions"])
    atomic_numbers = np.array(atoms.get_atomic_numbers())
    
    # Determine element types
    if element_filter is not None:
        element_types = sorted(element_filter)
    else:
        element_types = sorted(set(atomic_numbers))
    n_elements = len(element_types)
    
    # Prepare neighbor data
    neighbor_lists = d["neighbors"]
    neighbor_dists = d["neighbordist"]
    neighbor_diffs = d["diff"]
    
    results = {}
    feature_names = []
    
    # -------------------------------------------------------------------------
    # G2 Radial Symmetry Functions
    # -------------------------------------------------------------------------
    n_g2 = len(g2_params)
    g2_features = np.zeros((n, n_g2 * n_elements))
    
    for i in range(n):
        neigh_idx = neighbor_lists[i]
        if len(neigh_idx) == 0:
            continue
        
        neigh_dists = np.array(neighbor_dists[i])
        neigh_types = atomic_numbers[neigh_idx]
        
        # Apply cutoff
        mask = neigh_dists < cutoff
        neigh_dists = neigh_dists[mask]
        neigh_types = neigh_types[mask]
        
        if len(neigh_dists) == 0:
            continue
        
        fc = _cutoff_cosine(neigh_dists, cutoff)
        
        # For each element type
        for ie, elem_z in enumerate(element_types):
            elem_mask = (neigh_types == elem_z)
            if not np.any(elem_mask):
                continue
            
            r_elem = neigh_dists[elem_mask]
            fc_elem = fc[elem_mask]
            
            # For each G2 parameter set
            for ip, (eta, rs) in enumerate(g2_params):
                g2_val = np.sum(np.exp(-eta * (r_elem - rs)**2) * fc_elem)
                g2_features[i, ip * n_elements + ie] = g2_val
    
    results["G2"] = g2_features
    for ip, (eta, rs) in enumerate(g2_params):
        for elem_z in element_types:
            feature_names.append(f"G2(η={eta},Rs={rs},Z={elem_z})")
    
    # -------------------------------------------------------------------------
    # G4 Angular Symmetry Functions
    # -------------------------------------------------------------------------
    if g4_params:
        # Element pairs (with ordering for j <= k)
        element_pairs = []
        for ie1, z1 in enumerate(element_types):
            for ie2, z2 in enumerate(element_types):
                if ie2 >= ie1:
                    element_pairs.append((z1, z2))
        n_pairs = len(element_pairs)
        
        n_g4 = len(g4_params)
        g4_features = np.zeros((n, n_g4 * n_pairs))
        
        for i in range(n):
            neigh_idx = neighbor_lists[i]
            if len(neigh_idx) < 2:
                continue
            
            neigh_dists = np.array(neighbor_dists[i])
            neigh_diffs = np.array(neighbor_diffs[i])
            neigh_types = atomic_numbers[neigh_idx]
            
            # Apply cutoff
            mask = neigh_dists < cutoff
            neigh_idx_filt = np.array(neigh_idx)[mask]
            neigh_dists = neigh_dists[mask]
            neigh_diffs = neigh_diffs[mask]
            neigh_types = neigh_types[mask]
            
            n_neigh = len(neigh_dists)
            if n_neigh < 2:
                continue
            
            fc = _cutoff_cosine(neigh_dists, cutoff)
            
            # Loop over unique pairs of neighbors
            for j_idx in range(n_neigh):
                r_ij = neigh_dists[j_idx]
                vec_ij = neigh_diffs[j_idx]
                type_j = neigh_types[j_idx]
                fc_ij = fc[j_idx]
                
                for k_idx in range(j_idx + 1, n_neigh):
                    r_ik = neigh_dists[k_idx]
                    vec_ik = neigh_diffs[k_idx]
                    type_k = neigh_types[k_idx]
                    fc_ik = fc[k_idx]
                    
                    # Compute angle via dot product
                    cos_theta = np.dot(vec_ij, vec_ik) / (r_ij * r_ik + 1e-10)
                    cos_theta = np.clip(cos_theta, -1.0, 1.0)
                    
                    # Compute j-k distance
                    # vec_jk = pos_k - pos_j = (pos_k - pos_i) - (pos_j - pos_i) = vec_ik - vec_ij
                    vec_jk = vec_ik - vec_ij
                    r_jk = np.linalg.norm(vec_jk)
                    fc_jk = _cutoff_cosine(np.array([r_jk]), cutoff)[0]
                    
                    # Determine element pair index (sorted)
                    if type_j <= type_k:
                        pair = (type_j, type_k)
                    else:
                        pair = (type_k, type_j)
                    
                    try:
                        pair_idx = element_pairs.index(pair)
                    except ValueError:
                        continue
                    
                    # For each G4 parameter set
                    for ip, (eta, zeta, lam) in enumerate(g4_params):
                        angular = (1.0 + lam * cos_theta) ** zeta
                        radial = np.exp(-eta * (r_ij**2 + r_ik**2 + r_jk**2))
                        prefactor = 2.0 ** (1 - zeta)
                        g4_val = prefactor * angular * radial * fc_ij * fc_ik * fc_jk
                        g4_features[i, ip * n_pairs + pair_idx] += g4_val
        
        results["G4"] = g4_features
        for ip, (eta, zeta, lam) in enumerate(g4_params):
            for z1, z2 in element_pairs:
                feature_names.append(f"G4(η={eta},ζ={zeta},λ={lam},Z={z1}-{z2})")
    
    # -------------------------------------------------------------------------
    # G5 Angular Symmetry Functions (no j-k term)
    # -------------------------------------------------------------------------
    if g5_params:
        element_pairs = []
        for ie1, z1 in enumerate(element_types):
            for ie2, z2 in enumerate(element_types):
                if ie2 >= ie1:
                    element_pairs.append((z1, z2))
        n_pairs = len(element_pairs)
        
        n_g5 = len(g5_params)
        g5_features = np.zeros((n, n_g5 * n_pairs))
        
        for i in range(n):
            neigh_idx = neighbor_lists[i]
            if len(neigh_idx) < 2:
                continue
            
            neigh_dists = np.array(neighbor_dists[i])
            neigh_diffs = np.array(neighbor_diffs[i])
            neigh_types = atomic_numbers[neigh_idx]
            
            # Apply cutoff
            mask = neigh_dists < cutoff
            neigh_dists = neigh_dists[mask]
            neigh_diffs = neigh_diffs[mask]
            neigh_types = neigh_types[mask]
            
            n_neigh = len(neigh_dists)
            if n_neigh < 2:
                continue
            
            fc = _cutoff_cosine(neigh_dists, cutoff)
            
            for j_idx in range(n_neigh):
                r_ij = neigh_dists[j_idx]
                vec_ij = neigh_diffs[j_idx]
                type_j = neigh_types[j_idx]
                fc_ij = fc[j_idx]
                
                for k_idx in range(j_idx + 1, n_neigh):
                    r_ik = neigh_dists[k_idx]
                    vec_ik = neigh_diffs[k_idx]
                    type_k = neigh_types[k_idx]
                    fc_ik = fc[k_idx]
                    
                    cos_theta = np.dot(vec_ij, vec_ik) / (r_ij * r_ik + 1e-10)
                    cos_theta = np.clip(cos_theta, -1.0, 1.0)
                    
                    if type_j <= type_k:
                        pair = (type_j, type_k)
                    else:
                        pair = (type_k, type_j)
                    
                    try:
                        pair_idx = element_pairs.index(pair)
                    except ValueError:
                        continue
                    
                    for ip, (eta, zeta, lam) in enumerate(g5_params):
                        angular = (1.0 + lam * cos_theta) ** zeta
                        radial = np.exp(-eta * (r_ij**2 + r_ik**2))
                        prefactor = 2.0 ** (1 - zeta)
                        g5_val = prefactor * angular * radial * fc_ij * fc_ik
                        g5_features[i, ip * n_pairs + pair_idx] += g5_val
        
        results["G5"] = g5_features
        for ip, (eta, zeta, lam) in enumerate(g5_params):
            for z1, z2 in element_pairs:
                feature_names.append(f"G5(η={eta},ζ={zeta},λ={lam},Z={z1}-{z2})")
    
    # -------------------------------------------------------------------------
    # Combine all features
    # -------------------------------------------------------------------------
    feature_list = [results["G2"]]
    if "G4" in results:
        feature_list.append(results["G4"])
    if "G5" in results:
        feature_list.append(results["G5"])
    
    features = np.hstack(feature_list)
    results["features"] = features
    results["element_types"] = element_types
    results["feature_names"] = feature_names
    
    # Store on atoms
    atoms.arrays["pyscal_acsf"] = features
    atoms.info["pyscal_acsf_params"] = {
        "cutoff": cutoff,
        "g2_params": g2_params,
        "g4_params": g4_params,
        "g5_params": g5_params,
        "element_types": element_types,
        "n_features": features.shape[1],
    }
    
    return results


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
