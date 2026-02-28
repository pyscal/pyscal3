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
# Coulomb Matrix
# ---------------------------------------------------------------------------

def coulomb_matrix(atoms: Atoms, representation='eigenvalues', max_atoms=None, 
                   use_atomic_numbers=True):
    """
    Compute the Coulomb matrix descriptor.
    
    The Coulomb matrix encodes electrostatic interactions between atoms,
    originally introduced by Rupp et al. (2012) for predicting molecular
    atomization energies.
    
    Matrix elements:
    - Diagonal: M_ii = 0.5 * Z_i^2.4 (atomic self-energy proxy)
    - Off-diagonal: M_ij = Z_i * Z_j / |R_i - R_j| (Coulomb repulsion)
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure to compute descriptor for.
    representation : str, default 'eigenvalues'
        How to convert matrix to fixed-size descriptor:
        - 'eigenvalues': Sorted eigenvalues (permutation invariant)
        - 'sorted': Sorted by row norm (upper triangle)
        - 'matrix': Raw matrix (only for same-size systems)
    max_atoms : int, optional
        Maximum number of atoms for padding. Required for 'eigenvalues'
        and 'sorted' representations. If None, uses len(atoms).
    use_atomic_numbers : bool, default True
        If True, use atomic numbers. If False, use 1 for all atoms.
        
    Returns
    -------
    ndarray
        Coulomb matrix descriptor:
        - 'eigenvalues': shape (max_atoms,)
        - 'sorted': shape (max_atoms * (max_atoms + 1) // 2,)
        - 'matrix': shape (N, N)
        
    Notes
    -----
    The Coulomb matrix is:
    - Translationally invariant
    - Rotationally invariant (eigenvalues)
    - NOT permutationally invariant unless sorted
    
    References
    ----------
    .. [1] Rupp, M. et al. (2012). "Fast and Accurate Modeling of Molecular 
           Atomization Energies with Machine Learning." Phys. Rev. Lett. 108, 058301.
    
    Examples
    --------
    >>> from ase.build import molecule
    >>> import pyscal3
    >>> mol = molecule('H2O')
    >>> desc = pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=10)
    >>> print(desc.shape)
    (10,)
    """
    positions = atoms.get_positions()
    N = len(atoms)
    
    if use_atomic_numbers:
        Z = atoms.get_atomic_numbers().astype(float)
    else:
        Z = np.ones(N)
    
    # Compute distance matrix
    diff = positions[:, None, :] - positions[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    
    # Build Coulomb matrix
    # Off-diagonal: Z_i * Z_j / r_ij
    with np.errstate(divide='ignore', invalid='ignore'):
        M = np.outer(Z, Z) / dist
    
    # Diagonal: 0.5 * Z^2.4
    np.fill_diagonal(M, 0.5 * Z ** 2.4)
    
    # Handle inf on diagonal (already set above)
    M = np.nan_to_num(M, nan=0.0, posinf=0.0, neginf=0.0)
    
    if representation == 'matrix':
        atoms.arrays["pyscal_coulomb_matrix"] = M
        return M
    
    if max_atoms is None:
        max_atoms = N
    
    if representation == 'eigenvalues':
        eigenvalues = np.linalg.eigvalsh(M)
        eigenvalues = np.sort(eigenvalues)[::-1]  # Descending order
        
        # Pad with zeros
        result = np.zeros(max_atoms)
        n_eig = min(len(eigenvalues), max_atoms)
        result[:n_eig] = eigenvalues[:n_eig]
        
        atoms.arrays["pyscal_coulomb_eigenvalues"] = result
        return result
    
    elif representation == 'sorted':
        # Sort rows/cols by row norm
        row_norms = np.linalg.norm(M, axis=1)
        order = np.argsort(row_norms)[::-1]
        M_sorted = M[order][:, order]
        
        # Pad to max_atoms
        M_padded = np.zeros((max_atoms, max_atoms))
        n = min(N, max_atoms)
        M_padded[:n, :n] = M_sorted[:n, :n]
        
        # Extract upper triangle
        i_upper, j_upper = np.triu_indices(max_atoms)
        result = M_padded[i_upper, j_upper]
        
        atoms.info["pyscal_coulomb_sorted"] = result
        return result
    
    else:
        raise ValueError(f"Unknown representation: {representation}")


# ---------------------------------------------------------------------------
# MBTR (Many-Body Tensor Representation)
# ---------------------------------------------------------------------------

def _mbtr_k1(atomic_numbers, grid, sigma):
    """
    Compute k=1 MBTR term: atomic species distribution.
    
    Parameters
    ----------
    atomic_numbers : ndarray
        Atomic numbers of all atoms.
    grid : ndarray
        Grid points for histogram.
    sigma : float
        Gaussian broadening width.
        
    Returns
    -------
    ndarray
        k=1 MBTR values on grid.
    """
    Z = atomic_numbers.astype(float)
    
    # Gaussian broadening: sum over atoms
    diff = grid[:, None] - Z[None, :]
    gaussians = np.exp(-0.5 * (diff / sigma) ** 2)
    
    return gaussians.sum(axis=1)


def _mbtr_k2(positions, atomic_numbers, grid, sigma, weighting='exponential', 
             alpha=0.5, species_filter=None):
    """
    Compute k=2 MBTR term: inverse distance distribution.
    
    Parameters
    ----------
    positions : ndarray
        Atomic positions.
    atomic_numbers : ndarray
        Atomic numbers.
    grid : ndarray
        Grid points (1/r values).
    sigma : float
        Gaussian broadening width.
    weighting : str
        'unity' or 'exponential'.
    alpha : float
        Decay constant for exponential weighting.
    species_filter : tuple, optional
        (Z1, Z2) to filter specific pairs.
        
    Returns
    -------
    ndarray
        k=2 MBTR values on grid.
    """
    N = len(positions)
    if N < 2:
        return np.zeros(len(grid))
    
    # All pairwise distances
    diff = positions[:, None, :] - positions[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    
    # Upper triangle indices (i < j)
    i_idx, j_idx = np.triu_indices(N, k=1)
    r_ij = dist[i_idx, j_idx]
    
    # Species filter
    if species_filter is not None:
        Z1, Z2 = species_filter
        mask = (((atomic_numbers[i_idx] == Z1) & (atomic_numbers[j_idx] == Z2)) |
                ((atomic_numbers[i_idx] == Z2) & (atomic_numbers[j_idx] == Z1)))
        r_ij = r_ij[mask]
        i_idx = i_idx[mask]
        j_idx = j_idx[mask]
    
    if len(r_ij) == 0:
        return np.zeros(len(grid))
    
    # Filter out very small distances
    valid = r_ij > 1e-10
    r_ij = r_ij[valid]
    
    if len(r_ij) == 0:
        return np.zeros(len(grid))
    
    # Inverse distances
    inv_r = 1.0 / r_ij
    
    # Weighting
    if weighting == 'exponential':
        weights = np.exp(-alpha * r_ij)
    else:
        weights = np.ones_like(r_ij)
    
    # Gaussian broadening
    diff_grid = grid[:, None] - inv_r[None, :]
    gaussians = np.exp(-0.5 * (diff_grid / sigma) ** 2)
    
    return (gaussians * weights[None, :]).sum(axis=1)


def _mbtr_k3(positions, atomic_numbers, grid, sigma, weighting='exponential',
             alpha=0.5, cutoff=5.0, species_filter=None):
    """
    Compute k=3 MBTR term: angular distribution.
    
    Parameters
    ----------
    positions : ndarray
        Atomic positions.
    atomic_numbers : ndarray
        Atomic numbers.
    grid : ndarray
        Grid points (angles in radians from 0 to pi).
    sigma : float
        Gaussian broadening width (radians).
    weighting : str
        'unity' or 'exponential'.
    alpha : float
        Decay constant for exponential weighting.
    cutoff : float
        Distance cutoff for considering triplets.
    species_filter : tuple, optional
        (Z_center, Z1, Z2) to filter specific triplets.
        
    Returns
    -------
    ndarray
        k=3 MBTR values on grid.
    """
    N = len(positions)
    if N < 3:
        return np.zeros(len(grid))
    
    # Distance matrix
    diff = positions[:, None, :] - positions[None, :, :]
    dist = np.linalg.norm(diff, axis=-1)
    
    angles = []
    weights = []
    
    # Loop over central atom i
    for i in range(N):
        # Check species filter for central atom
        if species_filter is not None:
            Zc, Z1, Z2 = species_filter
            if atomic_numbers[i] != Zc:
                continue
        
        # Find neighbors within cutoff
        neighbors = np.where((dist[i] < cutoff) & (dist[i] > 1e-10))[0]
        
        # Loop over pairs (j, k) with j < k
        for idx_j, j in enumerate(neighbors):
            for k in neighbors[idx_j + 1:]:
                # Species filter for j and k
                if species_filter is not None:
                    Z1, Z2 = species_filter[1], species_filter[2]
                    pair_Z = {atomic_numbers[j], atomic_numbers[k]}
                    filter_Z = {Z1, Z2}
                    if pair_Z != filter_Z:
                        continue
                
                # Vectors from i to j and i to k
                v_ij = positions[j] - positions[i]
                v_ik = positions[k] - positions[i]
                r_ij = dist[i, j]
                r_ik = dist[i, k]
                
                # Angle
                cos_theta = np.dot(v_ij, v_ik) / (r_ij * r_ik)
                cos_theta = np.clip(cos_theta, -1.0, 1.0)
                theta = np.arccos(cos_theta)
                
                angles.append(theta)
                
                # Weighting
                if weighting == 'exponential':
                    r_jk = dist[j, k]
                    w = np.exp(-alpha * (r_ij + r_ik + r_jk))
                else:
                    w = 1.0
                weights.append(w)
    
    if len(angles) == 0:
        return np.zeros(len(grid))
    
    angles = np.array(angles)
    weights = np.array(weights)
    
    # Gaussian broadening
    diff_grid = grid[:, None] - angles[None, :]
    gaussians = np.exp(-0.5 * (diff_grid / sigma) ** 2)
    
    return (gaussians * weights[None, :]).sum(axis=1)


def mbtr(atoms: Atoms, k=(1, 2), n_grid=100, sigma_k1=0.5, sigma_k2=0.05,
         sigma_k3=0.1, alpha=0.5, cutoff_k3=5.0, weighting='exponential',
         normalize=True, species_resolved=False):
    """
    Compute Many-Body Tensor Representation (MBTR) descriptors.
    
    MBTR provides a smooth, continuous representation of atomic environments
    using Gaussian-broadened distributions for different body orders.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure to compute descriptors for.
    k : tuple of int, default (1, 2)
        Body orders to include:
        - k=1: Atomic species distribution
        - k=2: Inverse distance distribution (RDF-like)
        - k=3: Angular distribution
    n_grid : int, default 100
        Number of grid points for each term.
    sigma_k1 : float, default 0.5
        Gaussian width for k=1 (atomic number units).
    sigma_k2 : float, default 0.05
        Gaussian width for k=2 (1/Angstrom units).
    sigma_k3 : float, default 0.1
        Gaussian width for k=3 (radians).
    alpha : float, default 0.5
        Decay constant for exponential weighting.
    cutoff_k3 : float, default 5.0
        Distance cutoff for k=3 triplets (Angstrom).
    weighting : str, default 'exponential'
        Weighting scheme: 'unity' or 'exponential'.
    normalize : bool, default True
        If True, normalize descriptor to unit norm.
    species_resolved : bool, default False
        If True, compute separate channels per species combination.
        If False, aggregate all species together.
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'k1': ndarray for k=1 term (if included)
        - 'k2': ndarray for k=2 term (if included)
        - 'k3': ndarray for k=3 term (if included)
        - 'full': concatenated descriptor
        - 'grids': dict of grid arrays
        
    Notes
    -----
    MBTR is:
    - Translationally invariant
    - Rotationally invariant
    - Permutationally invariant
    - Fixed-size regardless of system size
    
    References
    ----------
    .. [1] Huo, H. & Rupp, M. (2017). "Unified Representation of Molecules and 
           Crystals for Machine Learning." arXiv:1704.06439.
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
    >>> result = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
    >>> print(result['full'].shape)
    """
    positions = atoms.get_positions()
    atomic_numbers = atoms.get_atomic_numbers()
    species = np.unique(atomic_numbers)
    
    result = {'grids': {}}
    all_descriptors = []
    
    # Grid definitions
    # k=1: atomic number range
    grid_k1 = np.linspace(0, max(100, max(species) + 10), n_grid)
    # k=2: inverse distance range (0.01 to 2.0 1/Angstrom typical)
    grid_k2 = np.linspace(0.01, 2.0, n_grid)
    # k=3: angle range (0 to pi)
    grid_k3 = np.linspace(0.01, np.pi - 0.01, n_grid)
    
    if 1 in k:
        result['grids']['k1'] = grid_k1
        if species_resolved:
            k1_list = []
            for Z in species:
                mask = atomic_numbers == Z
                k1_term = _mbtr_k1(atomic_numbers[mask], grid_k1, sigma_k1)
                k1_list.append(k1_term)
            k1 = np.concatenate(k1_list)
        else:
            k1 = _mbtr_k1(atomic_numbers, grid_k1, sigma_k1)
        result['k1'] = k1
        all_descriptors.append(k1)
    
    if 2 in k:
        result['grids']['k2'] = grid_k2
        if species_resolved:
            k2_list = []
            for i, Z1 in enumerate(species):
                for Z2 in species[i:]:
                    k2_term = _mbtr_k2(positions, atomic_numbers, grid_k2,
                                       sigma_k2, weighting, alpha,
                                       species_filter=(Z1, Z2))
                    k2_list.append(k2_term)
            k2 = np.concatenate(k2_list)
        else:
            k2 = _mbtr_k2(positions, atomic_numbers, grid_k2, sigma_k2,
                          weighting, alpha)
        result['k2'] = k2
        all_descriptors.append(k2)
    
    if 3 in k:
        result['grids']['k3'] = grid_k3
        if species_resolved:
            k3_list = []
            for Zc in species:
                for i, Z1 in enumerate(species):
                    for Z2 in species[i:]:
                        k3_term = _mbtr_k3(positions, atomic_numbers, grid_k3,
                                           sigma_k3, weighting, alpha, cutoff_k3,
                                           species_filter=(Zc, Z1, Z2))
                        k3_list.append(k3_term)
            k3 = np.concatenate(k3_list)
        else:
            k3 = _mbtr_k3(positions, atomic_numbers, grid_k3, sigma_k3,
                          weighting, alpha, cutoff_k3)
        result['k3'] = k3
        all_descriptors.append(k3)
    
    # Concatenate all
    full = np.concatenate(all_descriptors)
    
    if normalize:
        norm = np.linalg.norm(full)
        if norm > 1e-10:
            full = full / norm
            for key in ['k1', 'k2', 'k3']:
                if key in result:
                    result[key] = result[key] / norm
    
    result['full'] = full
    
    # Store in atoms
    atoms.info["pyscal_mbtr"] = full
    
    return result