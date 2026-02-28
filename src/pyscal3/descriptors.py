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
# Translational Order Parameters
# ---------------------------------------------------------------------------

def _get_reciprocal_vectors(cell):
    """
    Compute reciprocal lattice vectors from real-space cell.
    
    G_i = 2π (a_j × a_k) / (a_i · (a_j × a_k))
    """
    cell = np.asarray(cell)
    a1, a2, a3 = cell[0], cell[1], cell[2]
    volume = np.dot(a1, np.cross(a2, a3))
    
    if abs(volume) < 1e-10:
        raise ValueError("Cell has zero or near-zero volume")
    
    G1 = 2 * np.pi * np.cross(a2, a3) / volume
    G2 = 2 * np.pi * np.cross(a3, a1) / volume
    G3 = 2 * np.pi * np.cross(a1, a2) / volume
    
    return np.array([G1, G2, G3])


def _estimate_nn_distance(positions, n_samples=100):
    """Estimate nearest-neighbor distance from positions."""
    from scipy.spatial import cKDTree
    
    if len(positions) < 2:
        return 1.0
    
    tree = cKDTree(positions)
    # Sample atoms and find nearest neighbor
    n_samples = min(n_samples, len(positions))
    indices = np.random.choice(len(positions), n_samples, replace=False)
    
    nn_dists = []
    for i in indices:
        dist, _ = tree.query(positions[i], k=2)  # k=2: self and nearest
        nn_dists.append(dist[1])
    
    return np.mean(nn_dists)


def translational_order(atoms: Atoms, method='displacement', G_vectors=None,
                       reference=None, sigma=None):
    """
    Compute translational order parameter τ.
    
    The translational order parameter measures how well atoms maintain
    positions relative to an ideal crystalline lattice.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure to analyze.
    method : str, default 'displacement'
        Calculation method:
        - 'displacement': Uses displacement from reference positions (default)
        - 'fourier': Uses reciprocal lattice vectors (advanced)
    G_vectors : array-like, optional
        Custom reciprocal lattice vectors for Fourier method.
        Should be fundamental reciprocal lattice vectors of the 
        crystal (not supercell). Required for 'fourier' method.
    reference : ase.Atoms, optional
        Reference structure for displacement method. Must have same
        number of atoms.
    sigma : float, optional
        Gaussian width for displacement method. If None, uses 0.1
        times the estimated nearest-neighbor distance.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'tau': per-atom translational order (N,)
        - 'tau_global': global order parameter
        - 'mean': mean of per-atom τ
        - 'std': standard deviation
        - 'method': calculation method used
        
    Notes
    -----
    **Displacement method** (recommended): τ_i = exp(-u_i²/2σ²) where 
    u_i is the displacement from reference position. Values range from 0 to 1.
    A perfect match gives τ = 1, large displacements give τ → 0.
    
    **Fourier method**: τ = (1/N) |Σ_j exp(iG·r_j)| computed using 
    reciprocal lattice vectors G. This is structure-dependent and 
    requires appropriate G vectors for the crystal type.
    
    Typical values:
    - Perfect crystal: τ ≈ 1
    - Crystal near melting: τ ≈ 0.3-0.7
    - Liquid: τ ≈ 0
    
    References
    ----------
    .. [1] Hansen, J.P. & McDonald, I.R. "Theory of Simple Liquids"
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    >>> reference = atoms.copy()
    >>> # Add thermal displacements
    >>> atoms.positions += 0.1 * np.random.randn(*atoms.positions.shape)
    >>> result = pyscal3.translational_order(atoms, reference=reference)
    >>> print(f"Global τ: {result['tau_global']:.3f}")
    """
    positions = atoms.get_positions()
    cell = np.array(atoms.get_cell())
    N = len(atoms)
    
    if method == 'fourier':
        if G_vectors is None:
            raise ValueError(
                "G_vectors required for Fourier method. "
                "Provide reciprocal lattice vectors of the crystal "
                "(not supercell). Use method='displacement' for a "
                "simpler approach."
            )
        G_vectors = np.asarray(G_vectors)
        if G_vectors.ndim == 1:
            G_vectors = G_vectors.reshape(1, -1)
        
        # Compute global τ for each G vector (|⟨exp(iG·r)⟩|)
        tau_global_components = []
        tau_per_atom_components = []
        
        for G in G_vectors:
            phases = np.dot(positions, G)
            
            # Global: |sum of complex exponentials| / N
            complex_sum = np.sum(np.exp(1j * phases))
            tau_G_global = np.abs(complex_sum) / N
            tau_global_components.append(tau_G_global)
            
            # Per-atom: cos²(phase/2) = (1 + cos(phase))/2
            # This gives 1 when atom is on a lattice plane (phase = 0 mod 2π)
            tau_G_local = (1 + np.cos(phases)) / 2
            tau_per_atom_components.append(tau_G_local)
        
        # Average over reciprocal vectors
        tau_global = np.mean(tau_global_components)
        tau = np.mean(tau_per_atom_components, axis=0)
        
    elif method == 'displacement':
        if reference is None:
            raise ValueError(
                "reference structure required for displacement method."
            )
        
        ref_positions = reference.get_positions()
        
        if len(ref_positions) != N:
            raise ValueError(
                f"Reference has {len(ref_positions)} atoms but structure "
                f"has {N}. They must match."
            )
        
        # Compute displacements with minimum image convention
        u = positions - ref_positions
        
        # Apply minimum image convention for PBC
        if np.any(atoms.pbc):
            cell_inv = np.linalg.inv(cell)
            fractional_u = u @ cell_inv
            fractional_u = fractional_u - np.round(fractional_u)
            u = fractional_u @ cell
        
        u_sq = np.sum(u**2, axis=1)
        
        if sigma is None:
            nn_dist = _estimate_nn_distance(ref_positions)
            sigma = 0.1 * nn_dist
        
        tau = np.exp(-u_sq / (2 * sigma**2))
        tau_global = float(np.mean(tau))
        
    else:
        raise ValueError(f"Unknown method: {method}. Use 'fourier' or 'displacement'.")
    
    atoms.arrays['pyscal_translational_order'] = tau
    
    return {
        'tau': tau,
        'tau_global': float(tau_global),
        'mean': float(np.mean(tau)),
        'std': float(np.std(tau)),
        'method': method
    }


def lindemann_parameter(atoms: Atoms, reference: Atoms = None,
                       nn_distance: float = None):
    """
    Compute Lindemann parameter δ.
    
    The Lindemann parameter measures RMS displacement relative to the
    nearest-neighbor distance. The Lindemann criterion states that
    melting occurs when δ ≈ 0.1-0.15.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure to analyze.
    reference : ase.Atoms
        Reference structure (ideal lattice positions).
    nn_distance : float, optional
        Nearest-neighbor distance. If None, estimated from reference.
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'delta': per-atom Lindemann parameter (N,)
        - 'global': global Lindemann parameter √⟨u²⟩/a
        - 'msd': mean squared displacement ⟨u²⟩
        - 'nn_distance': nearest-neighbor distance used
        
    Notes
    -----
    The Lindemann parameter is defined as:
    
    δ_i = |u_i| / a
    
    where u_i is the displacement from the ideal position and a is
    the nearest-neighbor distance.
    
    The global parameter is:
    δ = √⟨u²⟩ / a
    
    Lindemann criterion for melting: δ ≈ 0.1-0.15
    
    References
    ----------
    .. [1] Lindemann, F.A. (1910). "Über die Berechnung molekularer 
           Eigenfrequenzen." Physikalische Zeitschrift 11: 609.
    
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> import numpy as np
    >>> 
    >>> # Create structure with thermal displacements
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    >>> reference = atoms.copy()
    >>> atoms.positions += 0.1 * np.random.randn(*atoms.positions.shape)
    >>> 
    >>> result = pyscal3.lindemann_parameter(atoms, reference)
    >>> print(f"Global δ: {result['global']:.3f}")
    """
    if reference is None:
        raise ValueError("reference structure required for Lindemann parameter")
    
    positions = atoms.get_positions()
    ref_positions = reference.get_positions()
    cell = np.array(atoms.get_cell())
    N = len(atoms)
    
    if len(ref_positions) != N:
        raise ValueError(
            f"Reference has {len(ref_positions)} atoms but structure "
            f"has {N}. They must match."
        )
    
    # Compute displacements
    u = positions - ref_positions
    
    # Apply minimum image convention for PBC
    if np.any(atoms.pbc):
        cell_inv = np.linalg.inv(cell)
        fractional_u = u @ cell_inv
        fractional_u = fractional_u - np.round(fractional_u)
        u = fractional_u @ cell
    
    u_mag = np.linalg.norm(u, axis=1)
    u_sq = u_mag ** 2
    msd = np.mean(u_sq)
    
    if nn_distance is None:
        nn_distance = _estimate_nn_distance(ref_positions)
    
    delta = u_mag / nn_distance
    delta_global = np.sqrt(msd) / nn_distance
    
    atoms.arrays['pyscal_lindemann'] = delta
    
    return {
        'delta': delta,
        'global': float(delta_global),
        'msd': float(msd),
        'nn_distance': float(nn_distance)
    }


def mean_squared_displacement(atoms: Atoms, reference: Atoms):
    """
    Compute mean squared displacement from reference structure.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Current structure.
    reference : ase.Atoms
        Reference structure (e.g., initial positions).
        
    Returns
    -------
    dict
        Dictionary containing:
        - 'displacement': per-atom displacement magnitude (N,)
        - 'msd': global mean squared displacement
        - 'rmsd': root mean squared displacement
        
    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> import numpy as np
    >>> 
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    >>> reference = atoms.copy()
    >>> atoms.positions += 0.5 * np.random.randn(*atoms.positions.shape)
    >>> 
    >>> result = pyscal3.mean_squared_displacement(atoms, reference)
    >>> print(f"MSD: {result['msd']:.3f} Å²")
    """
    positions = atoms.get_positions()
    ref_positions = reference.get_positions()
    cell = np.array(atoms.get_cell())
    N = len(atoms)
    
    if len(ref_positions) != N:
        raise ValueError(
            f"Reference has {len(ref_positions)} atoms but structure "
            f"has {N}. They must match."
        )
    
    # Compute displacements
    u = positions - ref_positions
    
    # Apply minimum image convention for PBC
    if np.any(atoms.pbc):
        cell_inv = np.linalg.inv(cell)
        fractional_u = u @ cell_inv
        fractional_u = fractional_u - np.round(fractional_u)
        u = fractional_u @ cell
    
    u_mag = np.linalg.norm(u, axis=1)
    msd = np.mean(u_mag ** 2)
    
    atoms.arrays['pyscal_displacement'] = u_mag
    
    return {
        'displacement': u_mag,
        'msd': float(msd),
        'rmsd': float(np.sqrt(msd))
    }


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
