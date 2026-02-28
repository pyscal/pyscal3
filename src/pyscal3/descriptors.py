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
    get_box_params,
    atoms_to_dict,
    dict_to_atoms,
    ensure_neighbors,
    create_attribute,
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
            if arr.ndim >= 1 and len(arr) == n and arr.dtype.kind != "O":
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
        sync_keys.extend(
            [
                "q%d" % val,
                "q%d_real" % val,
                "q%d_imag" % val,
                "avg_q%d" % val,
            ]
        )
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
        if not hasattr(v, "__len__") or len(v) == 0:
            need_calc = True
            break
        # Check if it looks like a 2-D container (list-of-lists or 2-D array)
        first = v[0]
        if not (hasattr(first, "__len__") and not isinstance(first, str)):
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
    find_neighbors(atoms, method="number", nmax=nmax, assign_neighbor=True)
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


def entropy(
    atoms: Atoms, rm, sigma=0.2, rstart=0.001, h=0.001, local=False, average=False
):
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
    hist, bin_edges = np.histogram(
        distances, bins=bins, range=(rmin, rmax), density=True
    )

    edgewidth = abs(bin_edges[1] - bin_edges[0])
    r = bin_edges[:-1]
    n = len(atoms)
    volume = abs(np.linalg.det(atoms.cell))
    rho = n / volume

    shell_vols = (4.0 / 3.0) * np.pi * ((r + edgewidth) ** 3 - r**3)
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


def find_solids(
    atoms: Atoms,
    bonds=0.5,
    threshold=0.5,
    avgthreshold=0.6,
    cluster=True,
    q=6,
    cutoff=0,
    right=True,
):
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

    _sync_back(
        d,
        atoms,
        ["solid", "bonds", "sij", "avg_sij", "q%d" % q, "q%d_real" % q, "q%d_imag" % q],
    )

    if cluster:
        return find_clusters(
            atoms, condition=np.array(d["solid"]) > 0, cutoff=cutoff, d=d
        )
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
            atoms.arrays["pyscal_largest_cluster"] = cluster_ids == largest_id
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
        result = pc.calculate_average_over_neighbors(d, values.tolist(), include_self)
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
# Density Correlation Functions
# ---------------------------------------------------------------------------


def structure_factor(
    atoms: Atoms,
    k_max: float = 10.0,
    n_k: int = 100,
    n_samples: int = 50,
    method: str = "direct",
):
    """
    Compute static structure factor S(k).

    The structure factor is the Fourier transform of density correlations,
    directly comparable to X-ray/neutron scattering experiments.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure to analyze.
    k_max : float, default 10.0
        Maximum wavevector magnitude (Å⁻¹).
    n_k : int, default 100
        Number of k points.
    n_samples : int, default 50
        Number of random directions per k magnitude (for 'direct' method).
    method : str, default 'direct'
        Calculation method:
        - 'direct': Direct sum over k-vectors (more accurate for small systems)

    Returns
    -------
    dict
        Dictionary containing:
        - 'k': wavevector magnitudes (Å⁻¹)
        - 'S': structure factor values
        - 'S_0': extrapolated S(k→0), related to compressibility

    Notes
    -----
    The structure factor is defined as:

    S(k) = (1/N) |Σⱼ exp(ik·rⱼ)|²

    For isotropic systems, averaged over all directions with |k| = k.

    Physical interpretation:
    - S(k→0) ∝ isothermal compressibility
    - First peak position k₁ ≈ 2π/d where d is interatomic spacing
    - Peak height indicates degree of short-range order

    References
    ----------
    .. [1] Hansen, J.P. & McDonald, I.R. "Theory of Simple Liquids"

    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    >>> result = pyscal3.structure_factor(atoms, k_max=12)
    >>> print(f"First peak at k ≈ {result['k'][np.argmax(result['S'][5:])+5]:.2f} Å⁻¹")
    """
    positions = atoms.get_positions()
    N = len(atoms)

    k_values = np.linspace(0.1, k_max, n_k)
    S_k = np.zeros(n_k)

    for i, k_mag in enumerate(k_values):
        S_acc = 0.0

        for _ in range(n_samples):
            # Random direction on unit sphere
            theta = np.arccos(2 * np.random.random() - 1)
            phi = 2 * np.pi * np.random.random()
            k_vec = k_mag * np.array(
                [
                    np.sin(theta) * np.cos(phi),
                    np.sin(theta) * np.sin(phi),
                    np.cos(theta),
                ]
            )

            # Compute |ρ(k)|² / N
            phase = np.dot(positions, k_vec)
            rho_k = np.sum(np.exp(1j * phase))
            S_acc += np.abs(rho_k) ** 2 / N

        S_k[i] = S_acc / n_samples

    # Extrapolate S(k→0) using linear fit of small-k region
    k_fit_mask = k_values < 1.5
    if np.sum(k_fit_mask) >= 3:
        coeffs = np.polyfit(k_values[k_fit_mask], S_k[k_fit_mask], 1)
        S_0 = coeffs[1]  # y-intercept
    else:
        S_0 = S_k[0]

    atoms.info["pyscal_structure_factor"] = {"k_max": k_max, "n_k": n_k}

    return {
        "k": k_values,
        "S": S_k,
        "S_0": float(max(0, S_0)),  # S(0) should be non-negative
    }


def local_density(
    atoms: Atoms, method: str = "voronoi", cutoff: float = None, sigma: float = None
):
    """
    Compute local atomic density.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure to analyze. For 'voronoi' method, must have
        Voronoi tessellation computed.
    method : str, default 'voronoi'
        Calculation method:
        - 'voronoi': ρ_i = 1/V_i (inverse Voronoi volume)
        - 'neighbor': ρ_i = n_i / (4πr³/3) (neighbors in sphere)
        - 'gaussian': Gaussian-smoothed density
    cutoff : float, optional
        Cutoff radius for 'neighbor' method. Required if method='neighbor'.
    sigma : float, optional
        Smoothing width for 'gaussian' method. Required if method='gaussian'.

    Returns
    -------
    dict
        Dictionary containing:
        - 'density': per-atom density (N,)
        - 'mean': mean density
        - 'std': standard deviation
        - 'method': calculation method used

    Notes
    -----
    **Voronoi method**: Uses the reciprocal of Voronoi cell volume.
    Higher values indicate more densely packed regions.

    **Neighbor method**: Counts neighbors in a sphere of given radius.
    Simple but cutoff-dependent.

    **Gaussian method**: Smoothed density using Gaussian kernel.
    Good for continuous density fields.

    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    >>> pyscal3.find_neighbors(atoms, method="voronoi")
    >>> result = pyscal3.local_density(atoms, method='voronoi')
    >>> print(f"Mean density: {result['mean']:.4f} atoms/Å³")
    """
    N = len(atoms)

    if method == "voronoi":
        if "pyscal_voronoi_volume" not in atoms.arrays:
            raise ValueError(
                "Voronoi tessellation required. Run "
                "pyscal3.find_neighbors(atoms, method='voronoi') first."
            )
        volumes = atoms.arrays["pyscal_voronoi_volume"]
        density = 1.0 / volumes

    elif method == "neighbor":
        if cutoff is None:
            raise ValueError("cutoff required for 'neighbor' method")

        d = _get_dict_with_neighbors(atoms)
        neighbors = d["neighbors"]
        neighbordist = d["neighbordist"]

        # Count neighbors within cutoff
        counts = np.zeros(N)
        for i in range(N):
            count = 0
            for dist in neighbordist[i]:
                if dist <= cutoff:
                    count += 1
            counts[i] = count

        sphere_vol = (4 / 3) * np.pi * cutoff**3
        density = counts / sphere_vol

    elif method == "gaussian":
        if sigma is None:
            raise ValueError("sigma required for 'gaussian' method")

        positions = atoms.get_positions()
        cell = np.array(atoms.get_cell())

        cutoff_gauss = 4 * sigma  # Truncate at 4σ
        norm = 1.0 / (2 * np.pi * sigma**2) ** 1.5

        density = np.zeros(N)

        for i in range(N):
            for j in range(N):
                if i == j:
                    continue

                r_ij = positions[j] - positions[i]

                # Minimum image convention
                if np.any(atoms.pbc):
                    cell_inv = np.linalg.inv(cell)
                    frac = r_ij @ cell_inv
                    frac = frac - np.round(frac)
                    r_ij = frac @ cell

                dist = np.linalg.norm(r_ij)

                if dist < cutoff_gauss:
                    density[i] += norm * np.exp(-(dist**2) / (2 * sigma**2))
    else:
        raise ValueError(
            f"Unknown method: {method}. Use 'voronoi', 'neighbor', or 'gaussian'."
        )

    atoms.arrays["pyscal_local_density"] = density

    return {
        "density": density,
        "mean": float(np.mean(density)),
        "std": float(np.std(density)),
        "method": method,
    }


def density_fluctuations(atoms: Atoms, n_blocks: int = 5):
    """
    Compute density fluctuations by block analysis.

    Divides the simulation cell into n_blocks³ subcells and
    measures the variance in atom count.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure to analyze. Must have periodic boundary conditions.
    n_blocks : int, default 5
        Number of blocks per dimension.

    Returns
    -------
    dict
        Dictionary containing:
        - 'mean_N': average atoms per block
        - 'var_N': variance in atom count
        - 'normalized_variance': ⟨(ΔN)²⟩/⟨N⟩ (related to compressibility)
        - 'block_counts': atom counts per block

    Notes
    -----
    The normalized variance is related to the static structure factor:

    ⟨(ΔN)²⟩/⟨N⟩ = S(k→0) = ρ k_B T κ_T

    where κ_T is the isothermal compressibility.

    For ideal gas: normalized_variance = 1
    For liquids: typically 0.02-0.05
    For crystals: → 0 (hyperuniform)

    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(5)
    >>> result = pyscal3.density_fluctuations(atoms, n_blocks=4)
    >>> print(f"Normalized variance: {result['normalized_variance']:.4f}")
    """
    positions = atoms.get_positions()
    cell = np.array(atoms.get_cell())

    # Convert to fractional coordinates
    cell_inv = np.linalg.inv(cell)
    frac_pos = positions @ cell_inv

    # Wrap to [0, 1)
    frac_pos = frac_pos % 1.0

    # Scale to block indices
    block_indices = (frac_pos * n_blocks).astype(int)
    block_indices = np.clip(block_indices, 0, n_blocks - 1)

    # Count atoms per block
    block_counts = np.zeros((n_blocks, n_blocks, n_blocks))
    for ix, iy, iz in block_indices:
        block_counts[ix, iy, iz] += 1

    counts = block_counts.flatten()
    mean_N = np.mean(counts)
    var_N = np.var(counts)

    normalized_var = var_N / mean_N if mean_N > 0 else 0.0

    atoms.info["pyscal_density_fluctuations"] = {
        "n_blocks": n_blocks,
        "normalized_variance": normalized_var,
    }

    return {
        "mean_N": float(mean_N),
        "var_N": float(var_N),
        "normalized_variance": float(normalized_var),
        "block_counts": counts,
    }


def hyperuniformity(
    atoms: Atoms, k_max: float = 5.0, n_k: int = 50, k_fit_max: float = 1.0
):
    """
    Analyze hyperuniformity from structure factor.

    A system is hyperuniform if S(k) → 0 as k → 0, indicating
    suppressed large-scale density fluctuations.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure to analyze.
    k_max : float, default 5.0
        Maximum k for structure factor calculation.
    n_k : int, default 50
        Number of k points.
    k_fit_max : float, default 1.0
        Maximum k for power-law fit.

    Returns
    -------
    dict
        Dictionary containing:
        - 'S_k': structure factor result
        - 'alpha': power-law exponent (S ~ k^α)
        - 'A': power-law prefactor
        - 'is_hyperuniform': True if system appears hyperuniform
        - 'hyperuniform_class': 'I' (α>1), 'II' (α=1), 'III' (0<α<1), or None

    Notes
    -----
    Hyperuniform materials have unusual properties:
    - Density fluctuations scale as surface area, not volume
    - Include crystals, certain disordered systems
    - Important for photonic materials, optimal packing

    Hyperuniformity classes:
    - Class I: S(k) ~ k^α with α > 1
    - Class II: S(k) ~ k (linear)
    - Class III: S(k) ~ k^α with 0 < α < 1

    References
    ----------
    .. [1] Torquato, S. & Stillinger, F.H. (2003). "Local density
           fluctuations, hyperuniformity, and order metrics."
           Phys. Rev. E 68, 041113.

    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(5)
    >>> result = pyscal3.hyperuniformity(atoms)
    >>> print(f"Hyperuniform: {result['is_hyperuniform']}")
    >>> print(f"Class: {result['hyperuniform_class']}")
    """
    # Compute structure factor
    S_result = structure_factor(atoms, k_max=k_max, n_k=n_k)
    k_values = S_result["k"]
    S_k = S_result["S"]

    # Fit S(k) ~ A * k^α for small k
    mask = k_values <= k_fit_max
    k_fit = k_values[mask]
    S_fit = S_k[mask]

    # Log-log fit (avoid log(0))
    valid = S_fit > 1e-10
    if np.sum(valid) >= 3:
        log_k = np.log(k_fit[valid])
        log_S = np.log(S_fit[valid])
        coeffs = np.polyfit(log_k, log_S, 1)
        alpha = coeffs[0]
        A = np.exp(coeffs[1])
    else:
        alpha = 0.0
        A = S_fit[0] if len(S_fit) > 0 else 1.0

    # Determine hyperuniformity
    is_hyperuniform = alpha > 0.3 and S_result["S_0"] < 0.1

    if is_hyperuniform:
        if alpha > 1.5:
            hu_class = "I"
        elif 0.7 < alpha <= 1.5:
            hu_class = "II"
        else:
            hu_class = "III"
    else:
        hu_class = None

    atoms.info["pyscal_hyperuniformity"] = {
        "alpha": alpha,
        "is_hyperuniform": is_hyperuniform,
    }

    return {
        "S_k": S_result,
        "alpha": float(alpha),
        "A": float(A),
        "is_hyperuniform": is_hyperuniform,
        "hyperuniform_class": hu_class,
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
        d, 0.0, triclinic, rot, rotinv, boxdims, 2, nmax, (n > 250), False
    )
