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


# ---------------------------------------------------------------------------
# Topological / Graph-based Descriptors
# ---------------------------------------------------------------------------

def coordination_numbers(atoms: Atoms):
    """
    Compute coordination number for each atom.
    
    The coordination number is the count of neighbors within the cutoff
    used during neighbor finding.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
        
    Returns
    -------
    ndarray
        Coordination number for each atom.
        
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    >>> cn = pyscal3.coordination_numbers(atoms)
    >>> print(f"Mean coordination: {cn.mean():.1f}")
    """
    d = _get_dict_with_neighbors(atoms)
    
    cn = np.array([len(n) for n in d["neighbors"]])
    atoms.arrays["pyscal_coordination"] = cn
    
    return cn


def coordination_stats(atoms: Atoms):
    """
    Compute coordination number statistics.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'coordination': per-atom coordination numbers
        - 'mean': mean coordination
        - 'std': standard deviation
        - 'distribution': dict mapping coordination -> fraction
        
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    >>> stats = pyscal3.coordination_stats(atoms)
    >>> print(f"Mean: {stats['mean']:.2f}, Std: {stats['std']:.2f}")
    """
    cn = coordination_numbers(atoms)
    
    values, counts = np.unique(cn, return_counts=True)
    distribution = dict(zip(values.tolist(), (counts / len(cn)).tolist()))
    
    return {
        'coordination': cn,
        'mean': float(np.mean(cn)),
        'std': float(np.std(cn)),
        'distribution': distribution
    }


def bond_angle_distribution(atoms: Atoms, bins=90, range_deg=(0, 180)):
    """
    Compute the bond angle distribution.
    
    For each atom with at least 2 neighbors, computes all bond angles
    formed by pairs of neighbors. This is useful for characterizing
    local structure in amorphous materials.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    bins : int, default 90
        Number of histogram bins.
    range_deg : tuple, default (0, 180)
        Range of angles in degrees.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'angles': all bond angles (degrees)
        - 'bin_centers': histogram bin centers
        - 'histogram': normalized histogram values
        - 'mean': mean angle
        - 'std': standard deviation
        
    Notes
    -----
    Common peaks in crystalline structures:
    - FCC: 60°, 90°, 120°, 180° (for 1st shell neighbors)
    - BCC: ~70.5°, ~109.5° (for 1st shell)
    - Diamond: 109.5° (tetrahedral)
    
    References
    ----------
    .. [1] Stillinger, F.H. & Weber, T.A. (1985). "Computer simulation of 
           local order in condensed phases of silicon." Phys. Rev. B 31, 5262.
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Si", "diamond", cubic=True).repeat(2)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    >>> bad = pyscal3.bond_angle_distribution(atoms)
    >>> peak_angle = bad['bin_centers'][np.argmax(bad['histogram'])]
    >>> print(f"Peak angle: {peak_angle:.1f} degrees")
    """
    d = _get_dict_with_neighbors(atoms)
    positions = np.array(d['positions'])
    neighbors = d['neighbors']
    diffs = d['diff']
    
    angles = []
    
    for i in range(len(atoms)):
        neighs = neighbors[i]
        if len(neighs) < 2:
            continue
        
        diff_i = diffs[i]
        
        # All pairs of neighbors
        for j_idx in range(len(neighs)):
            for k_idx in range(j_idx + 1, len(neighs)):
                rij = np.array(diff_i[j_idx])
                rik = np.array(diff_i[k_idx])
                
                rij_norm = np.linalg.norm(rij)
                rik_norm = np.linalg.norm(rik)
                
                if rij_norm < 1e-10 or rik_norm < 1e-10:
                    continue
                
                cos_theta = np.dot(rij, rik) / (rij_norm * rik_norm)
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta_deg = np.degrees(np.arccos(cos_theta))
                
                angles.append(theta_deg)
    
    angles = np.array(angles)
    
    if len(angles) > 0:
        hist, bin_edges = np.histogram(angles, bins=bins, range=range_deg, density=True)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        mean_angle = float(np.mean(angles))
        std_angle = float(np.std(angles))
    else:
        bin_centers = np.linspace(range_deg[0], range_deg[1], bins)
        hist = np.zeros(bins)
        mean_angle = 0.0
        std_angle = 0.0
    
    result = {
        'angles': angles,
        'bin_centers': bin_centers,
        'histogram': hist,
        'mean': mean_angle,
        'std': std_angle
    }
    
    atoms.info["pyscal_bond_angle_distribution"] = {
        'mean': mean_angle, 'std': std_angle
    }
    
    return result


def clustering_coefficient(atoms: Atoms):
    """
    Compute local and global clustering coefficients.
    
    The clustering coefficient measures how connected the neighbors
    of each atom are to each other. High values indicate dense local
    networks, low values indicate sparse or chain-like connections.
    
    For atom i with k neighbors, the local clustering coefficient is:
    C_i = 2 * e_i / (k * (k-1))
    
    where e_i is the number of bonds between neighbors of i.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'local': per-atom clustering coefficients
        - 'global': mean clustering coefficient
        - 'transitivity': ratio of triangles to connected triples
        
    Notes
    -----
    Typical values:
    - FCC perfect crystal: 0.5-0.6 (depends on cutoff)  
    - Random network: ~0.1-0.2
    - Linear chain: 0.0
    
    References
    ----------
    .. [1] Watts, D.J. & Strogatz, S.H. (1998). "Collective dynamics of 
           'small-world' networks." Nature 393, 440.
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    >>> cc = pyscal3.clustering_coefficient(atoms)
    >>> print(f"Global clustering: {cc['global']:.3f}")
    """
    d = _get_dict_with_neighbors(atoms)
    neighbors = d['neighbors']
    n_atoms = len(atoms)
    
    # Convert to sets for O(1) lookup
    neighbor_sets = [set(n) for n in neighbors]
    
    local_cc = np.zeros(n_atoms)
    triangles = 0
    connected_triples = 0
    
    for i in range(n_atoms):
        neighs = list(neighbor_sets[i])
        k = len(neighs)
        
        if k < 2:
            local_cc[i] = 0.0
            continue
        
        # Count edges between neighbors of i
        edges_between = 0
        for j_idx in range(len(neighs)):
            for k_idx in range(j_idx + 1, len(neighs)):
                if neighs[k_idx] in neighbor_sets[neighs[j_idx]]:
                    edges_between += 1
        
        # Local clustering coefficient
        possible_edges = k * (k - 1) // 2
        local_cc[i] = edges_between / possible_edges if possible_edges > 0 else 0.0
        
        triangles += edges_between
        connected_triples += possible_edges
    
    global_cc = float(np.mean(local_cc))
    transitivity = triangles / connected_triples if connected_triples > 0 else 0.0
    
    atoms.arrays["pyscal_clustering_coefficient"] = local_cc
    
    return {
        'local': local_cc,
        'global': global_cc,
        'transitivity': float(transitivity)
    }


def ring_statistics(atoms: Atoms, max_ring_size=8):
    """
    Compute ring statistics using shortest-path method.
    
    Counts primitive rings of different sizes in the bond network.
    Important for characterizing network topology in glasses, zeolites,
    and other network-forming materials.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    max_ring_size : int, default 8
        Maximum ring size to search for.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'counts': dict mapping ring size to total count
        - 'per_atom': dict mapping ring size to count per atom
        - 'mean_size': mean ring size
        - 'rings': list of ring atom indices (if not too many)
        
    Notes
    -----
    The algorithm uses Franzblau's shortest-path method:
    For each bond (i,j), find the shortest path from i to j 
    that doesn't use the direct bond, forming a ring.
    
    Common ring sizes:
    - Silica glass: predominantly 5, 6, 7-membered rings
    - Crystalline quartz: 6-membered rings only
    - Diamond: 6-membered (chair) rings
    
    References
    ----------
    .. [1] Franzblau, D.S. (1991). "Computation of ring statistics for 
           network models of solids." Phys. Rev. B 44, 4925.
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Si", "diamond", cubic=True).repeat(2)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.8)
    >>> rings = pyscal3.ring_statistics(atoms, max_ring_size=8)
    >>> print(f"Ring counts: {rings['counts']}")
    """
    from collections import deque
    
    d = _get_dict_with_neighbors(atoms)
    neighbors = d['neighbors']
    n_atoms = len(atoms)
    
    # Convert to sets for O(1) lookup
    neighbor_sets = [set(n) for n in neighbors]
    
    def bfs_shortest_path(start, target, exclude_direct=True):
        """Find shortest path from start to target."""
        if start == target:
            return None
            
        visited = {start: None}
        queue = deque([start])
        
        while queue:
            current = queue.popleft()
            
            for neighbor in neighbor_sets[current]:
                # Skip direct bond if requested
                if exclude_direct and current == start and neighbor == target:
                    continue
                    
                if neighbor == target:
                    # Found target - reconstruct path
                    path = [neighbor]
                    node = current
                    while node is not None:
                        path.append(node)
                        node = visited[node]
                    return path[::-1]
                
                if neighbor not in visited:
                    visited[neighbor] = current
                    queue.append(neighbor)
        
        return None
    
    ring_counts = {n: 0 for n in range(3, max_ring_size + 1)}
    found_rings = set()
    rings_list = []
    
    # For each bond, find shortest path not using that bond
    for i in range(n_atoms):
        for j in neighbor_sets[i]:
            if j <= i:  # Process each bond once
                continue
            
            path = bfs_shortest_path(i, j, exclude_direct=True)
            
            if path and len(path) >= 3 and len(path) <= max_ring_size:
                # Canonical form for deduplication
                ring_canonical = tuple(sorted(path))
                
                if ring_canonical not in found_rings:
                    found_rings.add(ring_canonical)
                    ring_size = len(path)
                    ring_counts[ring_size] += 1
                    
                    if len(rings_list) < 10000:  # Limit stored rings
                        rings_list.append(list(path))
    
    # Per-atom ring counts
    per_atom = {n: c / n_atoms for n, c in ring_counts.items()}
    
    # Mean ring size
    total_rings = sum(ring_counts.values())
    if total_rings > 0:
        mean_size = sum(n * c for n, c in ring_counts.items()) / total_rings
    else:
        mean_size = 0.0
    
    atoms.info["pyscal_ring_statistics"] = {
        'counts': ring_counts,
        'mean_size': mean_size
    }
    
    return {
        'counts': ring_counts,
        'per_atom': per_atom,
        'mean_size': float(mean_size),
        'rings': rings_list
    }


def topological_descriptors(atoms: Atoms, compute_rings=False, max_ring_size=6):
    """
    Compute a comprehensive set of topological descriptors.
    
    This is a convenience function that computes multiple graph-based
    descriptors in one call.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    compute_rings : bool, default False
        Whether to compute ring statistics (can be slow for large systems).
    max_ring_size : int, default 6
        Maximum ring size if computing rings.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'coordination': coordination statistics
        - 'clustering': clustering coefficients
        - 'bond_angles': bond angle distribution statistics
        - 'rings': ring statistics (if compute_rings=True)
        
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    >>> topo = pyscal3.topological_descriptors(atoms)
    >>> print(f"Mean coordination: {topo['coordination']['mean']:.1f}")
    >>> print(f"Clustering: {topo['clustering']['global']:.3f}")
    """
    result = {
        'coordination': coordination_stats(atoms),
        'clustering': clustering_coefficient(atoms),
        'bond_angles': bond_angle_distribution(atoms)
    }
    
    if compute_rings:
        result['rings'] = ring_statistics(atoms, max_ring_size)
    
    return result
