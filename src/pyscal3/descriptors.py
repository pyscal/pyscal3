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
from scipy.spatial import cKDTree
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
# Wigner-Seitz defect analysis
# ---------------------------------------------------------------------------

def wigner_seitz_analysis(
    atoms: Atoms,
    reference: Atoms,
    affine_mapping: str = "none",
    per_type_occupancies: bool = False,
) -> dict:
    """
    Wigner-Seitz cell analysis for vacancy/interstitial detection.
    
    This assigns each atom in `atoms` to the nearest lattice site in `reference`
    and counts site occupancies. Sites with occupancy 0 are vacancies; sites
    with occupancy > 1 contain interstitials.
    
    Parameters
    ----------
    atoms : Atoms
        The configuration to analyze (containing defects).
    reference : Atoms
        The perfect/reference configuration defining lattice sites.
    affine_mapping : str, optional
        How to handle cell distortion:
        - "none" (default): Use positions as-is
        - "to_reference": Rescale `atoms` positions to match reference cell
    per_type_occupancies : bool, optional
        If True, track occupancy by atom type (for antisite detection).
    
    Returns
    -------
    dict
        Keys:
        - "vacancy_count": Number of sites with occupancy = 0
        - "interstitial_count": Total excess atoms (sum of max(occ-1, 0))
        - "occupancy": (N_ref,) array of site occupancies
        - "site_index": (N_atoms,) array mapping each atom to its assigned site
        - "vacancy_indices": indices of reference sites with vacancies
        - "interstitial_sites": indices of sites with excess atoms
        
        If per_type_occupancies=True, also includes:
        - "occupancy_by_type": dict mapping type → (N_ref,) occupancy array
    
    Notes
    -----
    Results are also stored on `atoms`:
    - atoms.arrays["pyscal_ws_site_index"]: site assignment
    - atoms.arrays["pyscal_ws_occupancy"]: occupancy of assigned site
    - atoms.info["pyscal_ws_vacancy_count"]
    - atoms.info["pyscal_ws_interstitial_count"]
    
    The algorithm uses nearest-neighbor search via cKDTree. For periodic
    systems, reference sites near cell boundaries are replicated to handle
    atoms that may have wrapped to different periodic images.
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> # Perfect FCC reference
    >>> ref = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> # Create vacancy by deleting atom
    >>> defected = ref.copy()
    >>> del defected[0]
    >>> result = pyscal3.wigner_seitz_analysis(defected, ref)
    >>> print(f"Vacancies: {result['vacancy_count']}")  # 1
    >>> print(f"Interstitials: {result['interstitial_count']}")  # 0
    """
    ref_pos = reference.get_positions()
    disp_pos = atoms.get_positions().copy()
    n_ref = len(reference)
    n_atoms = len(atoms)
    
    # Apply affine mapping if requested
    if affine_mapping == "to_reference":
        # Transform displaced positions to fractional coords using displaced cell,
        # then to Cartesian using reference cell
        disp_cell = atoms.get_cell()
        ref_cell = reference.get_cell()
        if disp_cell.any() and ref_cell.any():
            # Convert to fractional, then to reference Cartesian
            disp_frac = np.linalg.solve(disp_cell.T, disp_pos.T).T
            disp_pos = disp_frac @ ref_cell
    elif affine_mapping != "none":
        raise ValueError(f"affine_mapping must be 'none' or 'to_reference', got '{affine_mapping}'")
    
    # Handle periodicity by replicating reference sites near boundaries
    # We'll use the reference cell for PBC handling
    ref_cell = reference.get_cell()
    pbc = reference.get_pbc()
    
    if any(pbc) and ref_cell.any():
        # Build expanded reference with periodic images
        expanded_ref, expanded_indices = _expand_for_pbc(ref_pos, ref_cell, pbc)
    else:
        expanded_ref = ref_pos
        expanded_indices = np.arange(n_ref)
    
    # Build KD-tree on expanded reference
    tree = cKDTree(expanded_ref)
    
    # Find nearest reference site for each atom
    distances, nearest_expanded = tree.query(disp_pos, k=1)
    
    # Map back to original reference indices
    site_index = expanded_indices[nearest_expanded]
    
    # Count occupancies
    occupancy = np.bincount(site_index, minlength=n_ref)
    
    # Compute summary statistics
    vacancy_count = int(np.sum(occupancy == 0))
    interstitial_count = int(np.sum(np.maximum(occupancy - 1, 0)))
    vacancy_indices = np.where(occupancy == 0)[0]
    interstitial_sites = np.where(occupancy > 1)[0]
    
    result = {
        "occupancy": occupancy,
        "site_index": site_index,
        "vacancy_count": vacancy_count,
        "interstitial_count": interstitial_count,
        "vacancy_indices": vacancy_indices,
        "interstitial_sites": interstitial_sites,
    }
    
    # Per-type occupancies for antisite detection
    if per_type_occupancies:
        atom_types = atoms.get_chemical_symbols()
        unique_types = sorted(set(atom_types))
        occupancy_by_type = {}
        for t in unique_types:
            mask = np.array([s == t for s in atom_types])
            type_sites = site_index[mask]
            occupancy_by_type[t] = np.bincount(type_sites, minlength=n_ref)
        result["occupancy_by_type"] = occupancy_by_type
    
    # Store results on atoms
    atoms.arrays["pyscal_ws_site_index"] = site_index
    atoms.arrays["pyscal_ws_occupancy"] = occupancy[site_index]  # Per-atom view
    atoms.info["pyscal_ws_vacancy_count"] = vacancy_count
    atoms.info["pyscal_ws_interstitial_count"] = interstitial_count
    
    return result


def _expand_for_pbc(positions, cell, pbc, skin=3.0):
    """
    Expand reference positions with periodic images near boundaries.
    
    Parameters
    ----------
    positions : ndarray (N, 3)
        Original positions
    cell : ndarray (3, 3)
        Cell vectors (rows)
    pbc : array-like of bool
        Periodic boundary conditions
    skin : float
        Distance from boundary to include images (Angstroms)
    
    Returns
    -------
    expanded_positions : ndarray (M, 3)
        Positions including relevant periodic images
    original_indices : ndarray (M,)
        Index into original positions for each expanded position
    """
    cell = np.asarray(cell)
    pbc = np.asarray(pbc)
    n = len(positions)
    
    # Convert to fractional coordinates
    try:
        cell_inv = np.linalg.inv(cell)
    except np.linalg.LinAlgError:
        # Degenerate cell, return as-is
        return positions, np.arange(n)
    
    frac_pos = positions @ cell_inv
    
    # Determine which shifts to apply
    shifts = []
    for ix in ([-1, 0, 1] if pbc[0] else [0]):
        for iy in ([-1, 0, 1] if pbc[1] else [0]):
            for iz in ([-1, 0, 1] if pbc[2] else [0]):
                shifts.append([ix, iy, iz])
    shifts = np.array(shifts)
    
    # Just apply all 27 (or fewer) shifts and let KDTree handle it
    # This is simpler and the overhead is small for typical systems
    all_positions = []
    all_indices = []
    
    for shift in shifts:
        shifted_frac = frac_pos + shift
        shifted_cart = shifted_frac @ cell
        all_positions.append(shifted_cart)
        all_indices.append(np.arange(n))
    
    expanded_positions = np.vstack(all_positions)
    original_indices = np.concatenate(all_indices)
    
    return expanded_positions, original_indices


def identify_defect_atoms(
    atoms: Atoms,
    reference: Atoms,
    affine_mapping: str = "none",
) -> dict:
    """
    Identify which atoms are at vacancies, interstitials, or antisites.
    
    This is a convenience wrapper around wigner_seitz_analysis that
    returns masks for different defect types.
    
    Parameters
    ----------
    atoms : Atoms
        The configuration to analyze.
    reference : Atoms 
        The reference configuration.
    affine_mapping : str, optional
        Affine mapping mode (see wigner_seitz_analysis).
    
    Returns
    -------
    dict
        - "perfect_mask": bool array, atoms at singly-occupied sites
        - "interstitial_mask": bool array, atoms at multiply-occupied sites
        - "vacancy_positions": (N_vac, 3) positions of empty reference sites
        - "defect_summary": string description of defects found
    """
    result = wigner_seitz_analysis(atoms, reference, affine_mapping, 
                                    per_type_occupancies=True)
    
    occupancy = result["occupancy"]
    site_index = result["site_index"]
    
    # Atoms at singly-occupied sites are "perfect"
    perfect_mask = (occupancy[site_index] == 1)
    
    # Atoms at multiply-occupied sites are interstitials (or share with one)
    interstitial_mask = (occupancy[site_index] > 1)
    
    # Vacancy positions
    ref_pos = reference.get_positions()
    vacancy_positions = ref_pos[result["vacancy_indices"]]
    
    # Summary
    summary_parts = []
    if result["vacancy_count"] > 0:
        summary_parts.append(f"{result['vacancy_count']} vacancies")
    if result["interstitial_count"] > 0:
        summary_parts.append(f"{result['interstitial_count']} interstitials")
    
    # Check for antisites in multi-component systems
    if "occupancy_by_type" in result and len(result["occupancy_by_type"]) > 1:
        ref_types = reference.get_chemical_symbols()
        antisite_count = 0
        for site_idx in range(len(reference)):
            if occupancy[site_idx] == 1:
                # Single atom at this site - check type match
                site_type = ref_types[site_idx]
                atom_at_site = np.where(site_index == site_idx)[0]
                if len(atom_at_site) == 1:
                    atom_type = atoms.get_chemical_symbols()[atom_at_site[0]]
                    if atom_type != site_type:
                        antisite_count += 1
        if antisite_count > 0:
            summary_parts.append(f"{antisite_count} antisites")
    
    defect_summary = ", ".join(summary_parts) if summary_parts else "No defects"
    
    return {
        "perfect_mask": perfect_mask,
        "interstitial_mask": interstitial_mask,
        "vacancy_positions": vacancy_positions,
        "defect_summary": defect_summary,
        **result,  # Include all WS results
    }

