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

import math
import functools
import numbers
import warnings
import numpy as np
import itertools
from scipy.spatial import cKDTree
from ase import Atoms
from scipy.special import sph_harm_y

import pyscal3.csystem as pc
from pyscal3._bridge import (
    atoms_to_dict,
    dict_to_atoms,
    ensure_neighbors,
    create_attribute,
    neighbor_arrays,
    rows_from_flat,
    stored_per_atom,
    gc_paused,
)
from pyscal3.neighbors import find_neighbors, _geometry


def _padded(atoms, key, fill, dtype):
    """Per-atom rows of a per-bond quantity as a 2-D (N, max_nn) array."""
    offsets, nb = neighbor_arrays(atoms, key)
    counts = np.diff(offsets)
    n = len(atoms)
    max_nn = int(counts.max()) if n else 0
    out = np.full((n, max_nn), fill, dtype=dtype)
    rows = np.repeat(np.arange(n), counts)
    cols = np.arange(len(nb[key])) - np.repeat(offsets[:-1], counts)
    out[rows, cols] = nb[key]
    return out


def _get_neighbor_dists_padded(atoms):
    """Return per-atom neighbor distances as a 2-D (N, max_nn) array.

    Uses ``atoms.arrays["pyscal_neighbordist"]`` when the rows are stored
    there (all atoms have the same number of neighbors), otherwise the flat
    neighbor data. Padded entries are zero.
    """
    if "pyscal_neighbordist" in atoms.arrays:
        return atoms.arrays["pyscal_neighbordist"]
    return _padded(atoms, "neighbordist", 0.0, float)


def _get_neighbor_indices_padded(atoms):
    """Return per-atom neighbor indices as a 2-D (N, max_nn) array.

    Uses ``atoms.arrays["pyscal_neighbors"]`` when the rows are stored there,
    otherwise the flat neighbor data. Padded entries are -1.
    """
    if "pyscal_neighbors" in atoms.arrays:
        return atoms.arrays["pyscal_neighbors"]
    return _padded(atoms, "neighbors", -1, int)


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
            if 1 <= arr.ndim <= 2 and len(arr) == n and arr.dtype.kind != "O":
                atoms.arrays[store_key] = arr
                continue
        except (ValueError, TypeError):
            pass
        atoms.info[store_key] = d[key]


def _periodic_volume(atoms, name):
    """Cell volume for density-based descriptors; requires full periodicity."""
    volume = abs(np.linalg.det(np.asarray(atoms.cell)))
    if not np.all(atoms.pbc) or volume <= 0:
        raise ValueError(
            f"{name} needs the global density N/V and therefore a fully "
            "periodic cell (all pbc True, non-zero volume)."
        )
    return volume


def _as_int_list(l):
    """Normalise an int / numpy integer / iterable of them to a list of ints."""
    if isinstance(l, numbers.Integral):
        return [int(l)]
    return [int(v) for v in l]


def _compute_qlm(atoms: Atoms, d: dict, l: int, bonds=None):
    """Compute q_l and the q_lm parts of every atom into ``d``.

    ``bonds`` is the result of neighbor_arrays with theta, phi and
    neighborweight; pass it when computing several l to flatten only once.
    """
    if bonds is None:
        bonds = neighbor_arrays(atoms, "theta", "phi", "neighborweight")
    offsets, nb = bonds
    q, real, imag = pc.calculate_q_single(
        offsets, nb["theta"], nb["phi"], nb["neighborweight"], l
    )
    d["q%d" % l] = q
    d["q%d_real" % l] = real
    d["q%d_imag" % l] = imag


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
    ll = _as_int_list(l)

    d = {}
    bonds = neighbor_arrays(atoms, "theta", "phi", "neighborweight")
    for val in ll:
        _compute_qlm(atoms, d, val, bonds)
    if averaged:
        offsets, nb = neighbor_arrays(atoms, "neighbors")
        for val in ll:
            d["avg_q%d" % val] = pc.calculate_aq_single(
                offsets, nb["neighbors"], d["q%d_real" % val], d["q%d_imag" % val], val
            )
        result_keys = ["avg_q%d" % v for v in ll]
    else:
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
# Wigner W_l Parameter (Third-order rotational invariant)
# ---------------------------------------------------------------------------


def wigner_w_parameter(atoms: Atoms, l, averaged=False, normalized=True):
    """
    Calculate the third-order Steinhardt invariant W_l.

    W_l is the third-order rotational invariant of the bond-orientational
    order parameters, constructed by contracting q_lm with Wigner 3j
    symbols. It distinguishes crystal structures that have similar q_l
    values (e.g., FCC and HCP, which differ in the sign of W_4).

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed via find_neighbors.
    l : int or list of int
        Order(s) for the W parameter. Only even values give nonzero
        results (W_l = 0 for odd l when j1=j2=j3=l).
    averaged : bool, optional
        If True, compute neighbor-averaged W_l (Lechner-Dellago).
        Default False.
    normalized : bool, optional
        If True (default), return the normalized
        ``hat{W}_l = W_l / (sum_m |q_lm|^2)^(3/2)``. If False, return raw W_l.

    Returns
    -------
    list of numpy arrays
        One array per requested l value, each of shape (natoms,).

    References
    ----------
    .. [1] Steinhardt, Nelson & Ronchetti, Phys. Rev. B 28, 784 (1983).
    .. [2] Lechner & Dellago, J. Chem. Phys. 129, 114707 (2008).

    Notes
    -----
    Known values for hat{W}_6:
      - FCC: −0.01316
      - HCP: −0.01244
      - BCC: +0.01316
      - ICO (Mackay): −0.16975
      - Liquid: ≈ 0
    """
    ll = _as_int_list(l)

    d = {}

    # Ensure q_lm are computed first (W_l requires them)
    bonds = neighbor_arrays(atoms, "theta", "phi", "neighborweight")
    for val in ll:
        _compute_qlm(atoms, d, val, bonds)

    if averaged:
        offsets, nb = neighbor_arrays(atoms, "neighbors")
        for val in ll:
            d["avg_w%d" % val], d["avg_what%d" % val] = pc.calculate_aw_single(
                offsets, nb["neighbors"], d["q%d_real" % val], d["q%d_imag" % val], val
            )
        if normalized:
            result_keys = ["avg_what%d" % v for v in ll]
        else:
            result_keys = ["avg_w%d" % v for v in ll]
    else:
        for val in ll:
            d["w%d" % val], d["what%d" % val] = pc.calculate_w_single(
                d["q%d_real" % val], d["q%d_imag" % val], val
            )
        if normalized:
            result_keys = ["what%d" % v for v in ll]
        else:
            result_keys = ["w%d" % v for v in ll]

    # Sync all W-related keys back
    sync_keys = []
    for val in ll:
        sync_keys.extend(
            [
                "q%d" % val,
                "q%d_real" % val,
                "q%d_imag" % val,
                "w%d" % val,
                "what%d" % val,
                "avg_w%d" % val,
                "avg_what%d" % val,
            ]
        )
    _sync_back(d, atoms, sync_keys)

    return [np.array(d[k]) for k in result_keys]


# ---------------------------------------------------------------------------
# Minkowski Structure Metrics (Voronoi-weighted Steinhardt)
# ---------------------------------------------------------------------------

def minkowski_parameter(atoms: Atoms, l, voroexp=1, averaged=False):
    """
    Calculate Minkowski structure metrics :math:`q_l^{\\mathrm{Mink}}`.

    These are Voronoi face-area weighted Steinhardt parameters
    (Mickel *et al.*, J. Chem. Phys. **138**, 044501, 2013).  For each atom
    *i* and angular-momentum order *l* the metric is

    .. math::

        q_l^{\\mathrm{Mink}}(i) =
        \\sqrt{\\frac{4\\pi}{2l+1}
              \\sum_{m=-l}^{l}
              \\left|
                \\sum_j \\frac{A_j^\\alpha}{\\sum_k A_k^\\alpha}
                Y_{lm}(\\hat{\\mathbf{r}}_{ij})
              \\right|^2}

    where *A_j* is the Voronoi face area between atoms *i* and *j*, and
    *α* is the exponent ``voroexp``.

    This function performs Voronoi neighbor finding internally, so there is
    no need to call :func:`find_neighbors` beforehand (any existing neighbor
    data is overwritten).

    Parameters
    ----------
    atoms : ase.Atoms
        The atomic structure (periodic boundary conditions recommended).
    l : int or list of int
        Steinhardt parameter order(s), e.g. 6 or [4, 6].
    voroexp : int or float, optional
        Face-area weight exponent *α*.  Default 1.
    averaged : bool, optional
        If True, compute neighbor-averaged values. Default False.

    Returns
    -------
    list of numpy arrays
        One array per requested *l*, each of shape ``(natoms,)``.
        Values are also stored in ``atoms.arrays["pyscal_q{l}"]``
        (or ``"pyscal_avg_q{l}"`` when ``averaged=True``).

    Notes
    -----
    Unlike a plain ``steinhardt_parameter`` call after Voronoi neighbor
    finding, this function guarantees that Voronoi neighbors are
    (re)computed with the requested ``voroexp`` so results are
    reproducible regardless of prior state.

    References
    ----------
    W. Mickel, S. C. Kapfer, G. E. Schröder-Turk and K. Mecke,
    "Shortcomings of the bond orientational order parameters for the
    analysis of disordered particulate matter",
    *J. Chem. Phys.* **138**, 044501 (2013).
    `doi:10.1063/1.4774084 <https://doi.org/10.1063/1.4774084>`__
    """
    # Voronoi neighbor finding (overwrites any previous neighbors)
    find_neighbors(atoms, method='voronoi', voroexp=voroexp)

    # Delegate to regular Steinhardt — weights are already in neighborweight
    return steinhardt_parameter(atoms, l, averaged=averaged)


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

    # q_lm from the current neighbors: stored values may belong to an
    # earlier neighbor list
    _compute_qlm(atoms, d, q)

    offsets, nb = neighbor_arrays(atoms, "neighbors")
    real = np.ascontiguousarray(d["q%d_real" % q], dtype=float)
    imag = np.ascontiguousarray(d["q%d_imag" % q], dtype=float)
    d["disorder"] = pc.calculate_disorder(offsets, nb["neighbors"], real, imag, q)

    sync_keys = ["disorder", "q%d" % q, "q%d_real" % q, "q%d_imag" % q]

    if averaged:
        # Average disorder over neighbors in C++
        d["avg_disorder"] = pc.calculate_average_disorder(
            offsets, nb["neighbors"], d["disorder"]
        )
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
    # candidates: all atoms within 2 * (V / N)^(1/3); the 12 or 14 nearest are used
    structure, finished = pc.cna_structure(
        *_geometry(atoms, 2), 2.0,
        0.0 if lattice_constant is None else float(lattice_constant), 14,
    )
    if not finished:
        _warn_few_candidates(14)
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
    structure, finished = pc.diamond_structure_cna(*_geometry(atoms, 2), 2.0)
    if not finished:
        _warn_few_candidates(4)
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
    offsets, nb = neighbor_arrays(atoms, "diff")
    cs = pc.calculate_centrosymmetry(offsets, nb["diff"], nmax)
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

    if "neighborweight" not in d:
        # rows not stored (store_rows=False): rebuild them for the C++ routine
        offsets, nb = neighbor_arrays(atoms, "neighborweight")
        d["neighborweight"] = rows_from_flat(offsets, nb["neighborweight"])
    pc.calculate_voronoi_vector(d, edge_cutoff, area_cutoff)

    vv = np.array(d["vorovector"])
    atoms.arrays["pyscal_vorovector"] = vv
    return vv


# ---------------------------------------------------------------------------
# Entropy Parameter
# ---------------------------------------------------------------------------


def entropy(
    atoms: Atoms, rm, sigma=0.2, rstart=0.001, h=0.001, local=False, average=False,
    averaged=None,
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
        If True, use the local density of each atom,
        ``n_i / (4/3 pi r_c,i^3)`` with ``r_c,i`` its neighbor cutoff,
        instead of the global density N/V. Default False.
    average : bool, optional
        Compute neighbor-averaged entropy. Default False.
    averaged : bool, optional
        Alias of ``average`` for consistency with the other descriptors.

    Returns
    -------
    numpy array
        Per-atom entropy (or averaged entropy) values, also stored as
        ``atoms.arrays["pyscal_entropy"]`` / ``"pyscal_average_entropy"``.

    Notes
    -----
    Only atoms in the neighbor list contribute to :math:`g_m^i(r)`, so the
    neighbor cutoff should be at least ``rm``; a warning is issued
    otherwise.
    """
    if averaged is not None:
        average = averaged

    offsets, nb = neighbor_arrays(atoms, "neighbors", "neighbordist")
    d = {}

    n = len(atoms)
    volume = _periodic_volume(atoms, "entropy")
    kb = 1

    cutoffs = np.ascontiguousarray(stored_per_atom(atoms, "cutoff"), dtype=float)
    if cutoffs.size > 0 and np.max(cutoffs) > 0 and rm > np.max(cutoffs) * (1 + 1e-9):
        warnings.warn(
            "entropy: rm=%.3f is larger than the neighbor cutoff (%.3f). "
            "g(r) is zero beyond the cutoff, so the result depends on it; "
            "compute neighbors with a cutoff >= rm." % (rm, np.max(cutoffs)),
            UserWarning,
            stacklevel=2,
        )

    if local:
        rho = 0
    else:
        rho = n / volume

    d["entropy"] = pc.calculate_entropy(
        offsets, nb["neighbordist"], cutoffs, sigma, rho, rstart, rm, h, kb
    )

    sync_keys = ["entropy"]

    if average:
        d["average_entropy"] = pc.calculate_average_entropy(
            offsets, nb["neighbors"], d["entropy"]
        )
        sync_keys.append("average_entropy")

    _sync_back(d, atoms, sync_keys)

    if average:
        return np.array(d["average_entropy"])
    return np.array(d["entropy"])


# ---------------------------------------------------------------------------
# Short-Range Order
# ---------------------------------------------------------------------------


def _resolve_atomic_number(value):
    """Accept an atomic number or a chemical symbol."""
    if isinstance(value, str):
        from ase.data import atomic_numbers

        try:
            return int(atomic_numbers[value])
        except KeyError:
            raise ValueError(f"Unknown chemical symbol '{value}'") from None
    return int(value)


def short_range_order(atoms: Atoms, reference_type=None, compare_type=None, average=True):
    """
    Calculate the Warren-Cowley short-range order parameter.

    For atoms *i* of the reference type A,

    .. math::

        \\alpha_{AB}(i) = 1 - \\frac{p_{AB}(i)}{c_B}

    where :math:`p_{AB}(i)` is the fraction of neighbors of *i* that are of
    the compare type B and :math:`c_B` is the global concentration of B.
    :math:`\\alpha < 0` indicates chemical ordering (unlike neighbors
    preferred), :math:`\\alpha > 0` clustering, and :math:`\\alpha = 0`
    random mixing. Choosing B = A gives the like-pair parameter.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors computed.
    reference_type : int or str, optional
        Atomic number or chemical symbol of the reference species A.
        Default: the most abundant species.
    compare_type : int or str, optional
        Atomic number or chemical symbol of the species B counted among the
        neighbors. Default: the most abundant species other than A.
    average : bool, optional
        If True, return the mean over all atoms of type A. Default True.

    Returns
    -------
    float or numpy array
        System average, or per-atom values (NaN for atoms that are not of
        the reference type). Per-atom values are stored in
        ``atoms.arrays["pyscal_sro"]``.
    """
    offsets, nb = neighbor_arrays(atoms, "neighbors")

    numbers = atoms.get_atomic_numbers()
    unique, counts = np.unique(numbers, return_counts=True)
    by_abundance = [int(z) for z in unique[np.argsort(-counts, kind="stable")]]

    if reference_type is None:
        reference_type = by_abundance[0]
    ref = _resolve_atomic_number(reference_type)
    if compare_type is None:
        others = [z for z in by_abundance if z != ref]
        if not others:
            raise ValueError(
                "short_range_order needs at least two species; pass "
                "reference_type and compare_type explicitly for a single species."
            )
        compare_type = others[0]
    cmp = _resolve_atomic_number(compare_type)

    if ref not in unique:
        raise ValueError(f"No atoms of reference type {ref} in the structure")
    if cmp not in unique:
        raise ValueError(f"No atoms of compare type {cmp} in the structure")

    sro = pc.calculate_short_range_order(
        offsets, nb["neighbors"], np.ascontiguousarray(numbers, dtype=np.int64), ref, cmp
    )
    atoms.arrays["pyscal_sro"] = sro

    if average:
        return float(np.nanmean(sro))
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
        ``rdf`` is g(r) normalised such that it tends to 1 for an ideal gas
        and ``rho * int g(r) 4 pi r^2 dr`` is the number of neighbors;
        ``r`` holds the left edges of the bins.

    Notes
    -----
    This function recomputes the neighbor list with a fixed cutoff of
    ``rmax`` and overwrites any existing neighbor data on ``atoms``.
    """
    find_neighbors(atoms, method="cutoff", cutoff=rmax)

    distances = neighbor_arrays(atoms, "neighbordist")[1]["neighbordist"]
    counts, bin_edges = np.histogram(distances, bins=bins, range=(rmin, rmax))

    edgewidth = abs(bin_edges[1] - bin_edges[0])
    r = bin_edges[:-1]
    n = len(atoms)
    volume = _periodic_volume(atoms, "radial_distribution_function")
    rho = n / volume

    # g(r) = <number of pairs in shell> / (N * rho * V_shell)
    shell_vols = (4.0 / 3.0) * np.pi * ((r + edgewidth) ** 3 - r**3)
    rdf = counts / (n * rho * shell_vols)

    return rdf, r


# ---------------------------------------------------------------------------
# Angular Distribution Function
# ---------------------------------------------------------------------------

def angular_distribution_function(atoms: Atoms, bins=180):
    """
    Calculate the angular distribution function (ADF).

    The ADF is the histogram of all bond angles :math:`\\theta_{jik}` formed
    by pairs of neighbors (j, k) around each atom i.  It characterises the
    local angular environment independently of a specific order parameter.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    bins : int
        Number of histogram bins in the angle range [0, 180] degrees.
        Default 180 (1-degree resolution).

    Returns
    -------
    (adf, angles) : tuple of numpy.ndarray
        ``adf`` is the normalised probability density of bond angles, and
        ``angles`` is the array of bin left-edges in degrees.

    Notes
    -----
    Angles are computed from the cosine of the angle between displacement
    vectors, using the same C++ infrastructure as :func:`chi_params`.
    Results are also stored in ``atoms.info["pyscal_adf"]`` and
    ``atoms.info["pyscal_adf_angles"]``.
    """
    # Use chi_params(angles=True) to get all pairwise cosines
    _, cosines_list = chi_params(atoms, angles=True)

    # Flatten all cosines and convert to degrees
    all_cosines = np.concatenate([np.array(c) for c in cosines_list])
    # Clamp to [-1, 1] to avoid NaN from acos
    all_cosines = np.clip(all_cosines, -1.0, 1.0)
    all_angles = np.degrees(np.arccos(all_cosines))

    hist, bin_edges = np.histogram(all_angles, bins=bins, range=(0, 180),
                                   density=True)
    angles = bin_edges[:-1]

    atoms.info["pyscal_adf"] = hist
    atoms.info["pyscal_adf_angles"] = angles
    return hist, angles


def bond_length_distribution(atoms: Atoms, bins=100, rmin=None, rmax=None):
    """
    Calculate the bond-length distribution function (BLDF).

    Unlike the full radial distribution function, the BLDF histograms only
    the bonds defined by the current neighbor list (no shell-volume
    normalisation).

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    bins : int
        Number of histogram bins.  Default 100.
    rmin, rmax : float or None
        Distance range.  If None, determined from the data.

    Returns
    -------
    (bldf, r) : tuple of numpy.ndarray
        ``bldf`` is the normalised probability density, and ``r`` is the
        array of bin left-edges.  Results are also stored in
        ``atoms.info["pyscal_bldf"]`` and ``atoms.info["pyscal_bldf_r"]``.
    """
    ensure_neighbors(atoms)
    dists = _get_neighbor_dists_padded(atoms)
    mask = dists > 0
    all_dists = dists[mask].ravel()

    if rmin is None:
        rmin = all_dists.min() * 0.9
    if rmax is None:
        rmax = all_dists.max() * 1.1

    hist, bin_edges = np.histogram(all_dists, bins=bins, range=(rmin, rmax),
                                   density=True)
    r = bin_edges[:-1]

    atoms.info["pyscal_bldf"] = hist
    atoms.info["pyscal_bldf_r"] = r
    return hist, r


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
    offsets, nb = neighbor_arrays(atoms, "neighbordist", "diff")
    ang = pc.calculate_angular_criteria(offsets, nb["neighbordist"], nb["diff"])
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
    offsets, nb = neighbor_arrays(atoms, "diff")
    cp, cosines, cos_offsets = pc.calculate_chi_params(offsets, nb["diff"])
    atoms.arrays["pyscal_chiparams"] = cp

    if angles:
        with gc_paused():
            values = cosines.tolist()
            cosines_list = [
                values[a:b] for a, b in zip(cos_offsets[:-1].tolist(), cos_offsets[1:].tolist())
            ]
        atoms.info["pyscal_cosines"] = cosines_list
        return cp, cosines_list
    return cp


# ---------------------------------------------------------------------------
# Ackland-Jones Structure Classification
# ---------------------------------------------------------------------------

# Labels used by identify_ackland_jones
ACKLAND_OTHER = 0
ACKLAND_FCC = 1
ACKLAND_HCP = 2
ACKLAND_BCC = 3
ACKLAND_ICO = 4

_ACKLAND_NAMES = {0: "other", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico"}


def identify_ackland_jones(atoms: Atoms):
    """
    Classify atomic environments with the Ackland–Jones method.

    For each atom, :math:`r_0^2` is the mean squared distance of its six
    nearest atoms. The bond angles between the :math:`N_0` atoms with
    :math:`d^2 < 1.45\\,r_0^2` are counted in eight ranges of
    :math:`\\cos\\theta` (:math:`\\chi_0, \\dots, \\chi_7`), and
    :math:`N_1` is the number of atoms with :math:`d^2 < 1.55\\,r_0^2`.
    The structure is assigned from the deviations of the counts from those
    of perfect bcc, close packed, fcc and hcp environments, as in
    Ackland & Jones (2006) and the original implementation of LAMMPS
    ``compute ackland/atom``. Atoms that match no structure are labelled
    0 (other).

    The function finds its own neighbors. A neighbor list stored on
    ``atoms`` is not used and not changed.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure.

    Returns
    -------
    labels : numpy.ndarray of int, shape (natoms,)
        Per-atom structure label, with the codes of
        :func:`common_neighbor_analysis`:

        =====  ==========
        Value  Structure
        =====  ==========
        0      other / unknown
        1      FCC
        2      HCP
        3      BCC
        4      ICO (icosahedral)
        =====  ==========

    names : list of str
        Name for each atom: ``"fcc"``, ``"hcp"``, ``"bcc"``, ``"ico"`` or
        ``"other"``.

    Notes
    -----
    The labels are stored in ``atoms.arrays["pyscal_ackland_label"]`` and
    ``atoms.arrays["pyscal_structure"]`` (the key that
    :func:`common_neighbor_analysis` also uses), and the angle counts in
    ``atoms.arrays["pyscal_ackland_chi"]``, shape (natoms, 8).

    Bins of :math:`\\cos\\theta`: [-1, -0.945), [-0.945, -0.915),
    [-0.915, -0.755), [-0.755, -0.195), [-0.195, 0.195), [0.195, 0.245),
    [0.245, 0.795), [0.795, 1].

    References
    ----------
    G. J. Ackland and A. P. Jones, "Applications of local crystal
    structure measures in experiment and simulation",
    *Phys. Rev. B* **73**, 054104 (2006).
    `doi:10.1103/PhysRevB.73.054104
    <https://doi.org/10.1103/PhysRevB.73.054104>`__
    """
    labels, chi = pc.ackland_jones_structure(*_geometry(atoms, 2), 2.0)
    labels = np.asarray(labels, dtype=int)
    names = [_ACKLAND_NAMES[l] for l in labels]
    atoms.arrays["pyscal_ackland_label"] = labels
    atoms.arrays["pyscal_structure"] = labels.copy()
    atoms.arrays["pyscal_ackland_chi"] = chi
    return labels, names


# ---------------------------------------------------------------------------
# Deformation Descriptors (require reference configuration)
# ---------------------------------------------------------------------------

def _local_deformation(atoms_cur, atoms_ref):
    """Strain tensor, D^2_min and slip vector of every atom.

    Neighbors that appear in the neighbor lists of an atom in both
    configurations are paired (the first entry of each neighbor index in
    either list); the affine deformation gradient is fitted to the paired
    neighbor vectors.
    """
    ensure_neighbors(atoms_cur)
    ensure_neighbors(atoms_ref)
    off_cur, cur = neighbor_arrays(atoms_cur, "neighbors", "diff")
    off_ref, ref = neighbor_arrays(atoms_ref, "neighbors", "diff")
    if len(off_cur) != len(off_ref):
        raise ValueError("atoms and reference must have the same number of atoms")
    return pc.calculate_local_deformation(
        off_cur, cur["neighbors"], cur["diff"], off_ref, ref["neighbors"], ref["diff"]
    )


def atomic_strain(atoms: Atoms, reference: Atoms):
    """
    Calculate the atomic (Green-Lagrange) strain tensor for each atom.

    The local deformation gradient F is computed by minimizing the
    squared difference between reference and deformed neighbor vectors.
    The Lagrangian strain is then E = (F^T F - I) / 2.

    Parameters
    ----------
    atoms : ase.Atoms
        Deformed configuration with neighbors computed.
    reference : ase.Atoms
        Reference (undeformed) configuration with the same neighbor list.

    Returns
    -------
    numpy.ndarray of shape (natoms, 3, 3)
        Green-Lagrange strain tensor per atom.
        Also stored as ``atoms.arrays["pyscal_strain"]``.

    Notes
    -----
    Both configurations must have neighbors computed with identical
    neighbor lists (same method and cutoff).
    Reference: Falk & Langer, PRE 57 (1998) 7192 (D^2_min);
    Shimizu, Ogata, Li, Mat. Trans. 48 (2007) 2923 (atomic strain).
    """
    strain, _, _ = _local_deformation(atoms, reference)
    atoms.arrays["pyscal_strain"] = strain
    return strain


def von_mises_strain(atoms: Atoms, reference: Atoms):
    """
    Compute the von Mises shear strain invariant from the atomic strain.

    Following Shimizu, Ogata and Li (2007),

    .. math::

        \\eta^{\\text{Mises}} = \\sqrt{
            \\frac{1}{6} \\left[
                (E_{xx}-E_{yy})^2 + (E_{yy}-E_{zz})^2 + (E_{zz}-E_{xx})^2
            \\right]
            + E_{xy}^2 + E_{yz}^2 + E_{xz}^2
        }

    which equals :math:`\\sqrt{E'_{ij} E'_{ij} / 2}` with :math:`E'` the
    deviatoric part of the strain, so it does not depend on the
    orientation of the axes.

    Parameters
    ----------
    atoms : ase.Atoms
        Deformed configuration.
    reference : ase.Atoms
        Reference configuration.

    Returns
    -------
    numpy.ndarray of shape (natoms,)
        Scalar von Mises strain per atom.
        Also stored as ``atoms.arrays["pyscal_von_mises"]``.
    """
    E = atomic_strain(atoms, reference)
    exx, eyy, ezz = E[:, 0, 0], E[:, 1, 1], E[:, 2, 2]
    exy, eyz, exz = E[:, 0, 1], E[:, 1, 2], E[:, 0, 2]
    vm = np.sqrt(
        ((exx - eyy)**2 + (eyy - ezz)**2 + (ezz - exx)**2) / 6.0
        + exy**2 + eyz**2 + exz**2
    )
    vm[np.isnan(E).any(axis=(1, 2))] = np.nan
    atoms.arrays["pyscal_von_mises"] = vm
    return vm


def d2min(atoms: Atoms, reference: Atoms):
    """
    Compute the D^2_min non-affine displacement (Falk & Langer 1998).

    D^2_min is the residual mean-squared displacement after subtracting
    the best-fit affine deformation.

    Parameters
    ----------
    atoms : ase.Atoms
        Deformed configuration with neighbors computed.
    reference : ase.Atoms
        Reference configuration with neighbors computed.

    Returns
    -------
    numpy.ndarray of shape (natoms,)
        D^2_min per atom (Angstrom^2).
        Also stored as ``atoms.arrays["pyscal_d2min"]``.
    """
    _, d2, _ = _local_deformation(atoms, reference)
    atoms.arrays["pyscal_d2min"] = d2
    return d2


def slip_vector(atoms: Atoms, reference: Atoms):
    """
    Compute the slip vector for each atom.

    The slip vector is the mean change of the neighbor vectors between the
    reference and the current configuration,

        s_i = (1/N_i) sum_j (r_ij - R_ij),

    taken over all neighbors j that appear in both neighbor lists, without
    fitting an affine transformation. Unlike the original definition of
    Zimmerman et al., no threshold is applied to select "slipped"
    neighbors, and the sign convention is r - R. For a homogeneous affine
    deformation the contributions of symmetric neighbor pairs cancel and
    the slip vector vanishes; it is non-zero where neighbors have moved
    relative to each other (dislocation cores, stacking faults).

    Parameters
    ----------
    atoms : ase.Atoms
        Deformed configuration.
    reference : ase.Atoms
        Reference configuration.

    Returns
    -------
    numpy.ndarray of shape (natoms, 3)
        Slip vector per atom.
        Also stored as ``atoms.arrays["pyscal_slip_vector"]``.

    Notes
    -----
    Ref: Zimmerman, Kelchner, Klein, Hamilton, Foiles, PRL 87 (2001) 165507.
    """
    _, _, slip = _local_deformation(atoms, reference)
    atoms.arrays["pyscal_slip_vector"] = slip
    return slip


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
    d = {}

    if isinstance(bonds, int):
        criteria = 0
    elif isinstance(bonds, float) and 0 <= bonds <= 1:
        criteria = 1
    else:
        raise TypeError("bonds must be int or float in [0,1]")

    compare_criteria = 0 if right else 1

    # Calculate Steinhardt parameters
    _compute_qlm(atoms, d, q)
    offsets, nb = neighbor_arrays(atoms, "neighbors")

    # Calculate bonds/solid classification
    d["bonds"], sij, d["avg_sij"], d["solid"] = pc.calculate_bonds(
        offsets, nb["neighbors"], d["q%d_real" % q], d["q%d_imag" % q], q,
        threshold, avgthreshold, bonds, compare_criteria, criteria,
    )
    d["sij"] = rows_from_flat(offsets, sij)

    _sync_back(
        d,
        atoms,
        ["solid", "bonds", "sij", "avg_sij", "q%d" % q, "q%d_real" % q, "q%d_imag" % q],
    )

    if cluster:
        return find_clusters(atoms, condition=np.array(d["solid"]) > 0, cutoff=cutoff)
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
    d : dict, optional
        Ignored. Kept so that existing calls keep working.

    Returns
    -------
    int or None
        Largest cluster size if largest=True, else None. Cluster ids are
        stored in ``atoms.arrays["pyscal_cluster"]`` (-1 for atoms that do
        not satisfy the condition) and, if largest=True, a boolean mask of
        the largest cluster in ``atoms.arrays["pyscal_largest_cluster"]``.
    """
    offsets, nb = neighbor_arrays(atoms, "neighbors", "neighbordist")
    condition = np.ascontiguousarray(condition, dtype=bool)
    cluster_ids = pc.find_clusters(
        offsets, nb["neighbors"], nb["neighbordist"],
        np.ascontiguousarray(stored_per_atom(atoms, "cutoff"), dtype=float),
        condition, cutoff,
    )
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
        Name of a per-atom property: a pyscal result with or without the
        ``pyscal_`` prefix (``"q6"`` and ``"pyscal_q6"`` are equivalent) or
        any other key of ``atoms.arrays``.
    include_self : bool
        Include the atom itself in the average. Default True.

    Returns
    -------
    numpy array
    """
    d = _get_dict_with_neighbors(atoms)

    # Find the data: pyscal keys are stored in d without the prefix
    plain = key[len("pyscal_"):] if key.startswith("pyscal_") else key
    if plain in d:
        values = np.array(d[plain])
    elif key in atoms.arrays:
        values = atoms.arrays[key]
    else:
        raise KeyError(f"Property '{key}' not found")

    # 1-D values: use fast C++ averaging
    values = np.asarray(values)
    if values.ndim == 1:
        offsets, nb = neighbor_arrays(atoms, "neighbors")
        return pc.calculate_average_over_neighbors(
            offsets, nb["neighbors"], np.ascontiguousarray(values, dtype=float), include_self
        )

    # Multi-dimensional: fall back to Python loop
    offsets, nb = neighbor_arrays(atoms, "neighbors")
    result = []
    for i in range(len(atoms)):
        vals = [values[i]] if include_self else []
        for j in nb["neighbors"][offsets[i]:offsets[i + 1]]:
            vals.append(values[j])
        result.append(np.mean(vals))

    return np.array(result)


# ---------------------------------------------------------------------------
# Coordination Variants
# ---------------------------------------------------------------------------

def effective_coordination_number(atoms: Atoms):
    """
    Calculate the effective coordination number (ECoN).

    ECoN weights each neighbor by a continuous function of its distance
    relative to a weighted mean distance (Hoppe 1979):

    .. math::

        \\mathrm{ECoN}_i = \\sum_j \\exp\\!\\left[
            1 - \\left(\\frac{d_{ij}}{d_{\\mathrm{av},i}}\\right)^6
        \\right],
        \\qquad
        d_{\\mathrm{av},i} = \\frac{\\sum_j d_{ij}\\, w_{ij}}{\\sum_j w_{ij}},
        \\quad
        w_{ij} = \\exp\\!\\left[
            1 - \\left(\\frac{d_{ij}}{d_{\\mathrm{av},i}}\\right)^6
        \\right].

    :math:`d_{\\mathrm{av},i}` is found by iteration, starting from the
    shortest neighbor distance, until it changes by less than one part in
    :math:`10^{12}`. The sums run over the stored neighbors of atom i.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.

    Returns
    -------
    numpy.ndarray, shape (natoms,)
        Per-atom ECoN values (0 for atoms without neighbors).  Also stored
        in ``atoms.arrays["pyscal_econ"]``.

    References
    ----------
    R. Hoppe, "Effective coordination numbers (ECoN) and mean fictive
    ionic radii (MEFIR)", *Z. Kristallogr.* **150**, 23 (1979).
    """
    ensure_neighbors(atoms)
    dists = _get_neighbor_dists_padded(atoms)   # (N, max_nn), 0 where padded
    mask = dists > 0
    has_neighbors = mask.any(axis=1)

    def weights(dav):
        return np.where(mask, np.exp(1.0 - (dists / dav[:, None]) ** 6), 0.0)

    # start from the shortest distance and iterate the weighted mean
    dav = np.where(has_neighbors, np.where(mask, dists, np.inf).min(axis=1), 1.0)
    for _ in range(200):
        w = weights(dav)
        wsum = w.sum(axis=1)
        new = np.where(wsum > 0, (w * dists).sum(axis=1) / np.where(wsum > 0, wsum, 1.0), dav)
        converged = np.all(np.abs(new - dav) <= 1e-12 * dav)
        dav = new
        if converged:
            break
    else:
        warnings.warn("effective_coordination_number: the mean distance did not "
                      "converge in 200 iterations", stacklevel=2)

    econ = np.where(has_neighbors, weights(dav).sum(axis=1), 0.0)
    atoms.arrays["pyscal_econ"] = econ
    return econ


def coordination_number(atoms: Atoms):
    """
    Return the simple coordination number (integer neighbor count).

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.

    Returns
    -------
    numpy.ndarray of int, shape (natoms,)
        Per-atom coordination number.  Also stored in
        ``atoms.arrays["pyscal_cn"]``.
    """
    ensure_neighbors(atoms)
    dists = _get_neighbor_dists_padded(atoms)
    cn = (dists > 0).sum(axis=1)
    atoms.arrays["pyscal_cn"] = cn
    return cn


def generalized_coordination_number(atoms: Atoms, cn_max=None):
    """
    Calculate the generalized coordination number (GCN).

    GCN accounts for the coordination of each neighbor, penalizing
    under-coordinated surface atoms:

    .. math::

        \\mathrm{GCN}_i = \\frac{1}{N_{\\max}}
        \\sum_{j \\in \\mathrm{neigh}(i)} \\mathrm{CN}_j

    where :math:`N_{\\max}` is the maximum (bulk) coordination number.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    cn_max : int or None
        Bulk coordination number.  If None, uses the maximum CN observed
        in the system.

    Returns
    -------
    numpy.ndarray, shape (natoms,)
        Per-atom GCN values.  Also stored in
        ``atoms.arrays["pyscal_gcn"]``.

    References
    ----------
    F. Calle-Vallejo *et al.*, "Finding optimal surface sites on
    heterogeneous catalysts by counting nearest neighbors",
    *Science* **350**, 185 (2015).
    `doi:10.1126/science.aab3501
    <https://doi.org/10.1126/science.aab3501>`__
    """
    cn = coordination_number(atoms)
    if cn_max is None:
        cn_max = cn.max()
    if cn_max == 0:
        cn_max = 1  # avoid division by zero

    neighbors = _get_neighbor_indices_padded(atoms)
    n = len(atoms)
    gcn = np.zeros(n, dtype=float)
    for i in range(n):
        nn_sum = 0.0
        count = 0
        for j in neighbors[i]:
            if j >= 0 and j < n:
                nn_sum += cn[j]
                count += 1
        if count > 0:
            gcn[i] = nn_sum / cn_max

    atoms.arrays["pyscal_gcn"] = gcn
    return gcn


def local_density(atoms: Atoms):
    """
    Estimate the local atomic number density.

    For each atom, the local density is estimated from the mean neighbor
    distance:

    .. math::

        \\rho_i = \\frac{N_i}{\\frac{4}{3}\\pi \\bar{d}_i^3}

    where :math:`N_i` is the coordination number and :math:`\\bar{d}_i`
    is the mean neighbor distance.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.

    Returns
    -------
    numpy.ndarray, shape (natoms,)
        Per-atom local density (atoms per unit volume).  Also stored in
        ``atoms.arrays["pyscal_local_density"]``.
    """
    ensure_neighbors(atoms)
    dists = _get_neighbor_dists_padded(atoms)
    mask = dists > 0

    cn = mask.sum(axis=1).astype(float)
    mean_d = np.where(cn > 0,
                      np.where(mask, dists, 0).sum(axis=1) / np.maximum(cn, 1),
                      1.0)  # avoid div by zero

    density = cn / (4.0 / 3.0 * np.pi * mean_d**3)
    atoms.arrays["pyscal_local_density"] = density
    return density


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _warn_few_candidates(nmax):
    """Warn that some atoms had fewer than ``nmax`` neighbor candidates.

    Atoms with fewer than ``nmax`` candidates cannot be classified and are
    labelled "others". That is the right answer for a surface or a small
    cluster, so this warns instead of failing the whole analysis.
    """
    warnings.warn(
        "Could not find %d neighbor candidates for every atom; those atoms "
        "are reported as 'others'. The structure may be a small cluster or "
        "very sparse. If it is meant to be periodic, check that atoms.pbc "
        "is set and the cell is correct." % nmax,
        RuntimeWarning,
        stacklevel=3,
    )


# ---------------------------------------------------------------------------
# ACE (Atomic Cluster Expansion) Descriptors
# ---------------------------------------------------------------------------

def _ace_cutoff(r, cutoff):
    """Smooth cosine cutoff function for ACE basis.
    
    f_cut(r) = 0.5 * (cos(pi * r / r_cut) + 1) for r < r_cut, else 0
    
    Parameters
    ----------
    r : float or array
        Distance(s).
    cutoff : float
        Cutoff radius.
        
    Returns
    -------
    float or array
        Cutoff function value(s).
    """
    r = np.asarray(r)
    result = np.where(r < cutoff, 0.5 * (np.cos(np.pi * r / cutoff) + 1), 0.0)
    return result


def _ace_radial_basis(n, r, cutoff, rmin=0.5):
    """Chebyshev polynomial radial basis for ACE.
    
    R_n(r) = T_n(x) * f_cut(r)
    where x = 2*(r - rmin)/(cutoff - rmin) - 1 maps r to [-1, 1]
    
    Parameters
    ----------
    n : int
        Basis function index (0, 1, 2, ...).
    r : float or array
        Distance(s).
    cutoff : float
        Cutoff radius.
    rmin : float
        Inner cutoff (default 0.5 Angstrom).
        
    Returns
    -------
    float or array
        Radial basis function value(s).
    """
    r = np.asarray(r)
    # Map r to [-1, 1]
    x = 2 * (r - rmin) / (cutoff - rmin) - 1
    x = np.clip(x, -1, 1)
    # Chebyshev polynomial T_n(x) = cos(n * arccos(x))
    Tn = np.cos(n * np.arccos(x))
    return Tn * _ace_cutoff(r, cutoff)


def _ace_a_functions(d, nmax, lmax, cutoff):
    """Compute A-basis (single-particle) coefficients for ACE.
    
    A_{i,nlm} = sum_{j in neighbors(i)} R_n(r_ij) * Y_l^m(r_ij_hat)
    
    These are the fundamental building blocks from which higher-order
    correlations are constructed.
    
    Parameters
    ----------
    d : dict
        Neighbor data dictionary with 'diff', 'neighbordist'.
    nmax : int
        Number of radial basis functions.
    lmax : int
        Maximum angular momentum quantum number.
    cutoff : float
        Cutoff radius.
        
    Returns
    -------
    A : ndarray, shape (natoms, nmax, lmax+1, 2*lmax+1), complex
        A-basis coefficients. A[i, n, l, m+lmax] gives A_{i,nlm}.
    """
    natoms = len(d['positions'])
    # Complex array to hold A coefficients
    # Index mapping: m ranges from -l to +l, stored at index m + lmax
    A = np.zeros((natoms, nmax, lmax + 1, 2 * lmax + 1), dtype=np.complex128)

    # all bonds, flat and in the order of the neighbor lists
    atom, dists, diffs = _flat_bonds(d, natoms)
    keep = ~((dists < 1e-10) | (dists >= cutoff))
    atom, rij, vec = atom[keep], dists[keep], diffs[keep]
    if len(rij) == 0:
        return A

    # Spherical coordinates
    # theta = polar angle from z axis
    # phi = azimuthal angle in xy plane
    theta = np.arccos(np.clip(vec[:, 2] / rij, -1, 1))
    phi = np.arctan2(vec[:, 1], vec[:, 0])

    radial = [_ace_radial_basis(n, rij, cutoff) for n in range(nmax)]
    for l in range(lmax + 1):
        for m in range(-l, l + 1):
            # scipy sph_harm_y(l, m, theta, phi) uses physics convention
            Y_lm = sph_harm_y(l, m, theta, phi)
            for n in range(nmax):
                # per-atom sums in bond order, as a loop over the bonds would do
                term = radial[n] * Y_lm
                A[:, n, l, m + lmax] = (
                    np.bincount(atom, weights=term.real, minlength=natoms)
                    + 1j * np.bincount(atom, weights=term.imag, minlength=natoms)
                )

    return A


def _flat_bonds(d, natoms):
    """Atom index, distance and vector of every bond in the atom dict ``d``."""
    if "bond_offsets" in d:
        counts = np.diff(np.asarray(d["bond_offsets"]))
        return (np.repeat(np.arange(natoms), counts),
                np.asarray(d["bond_distance"], dtype=float),
                np.asarray(d["bond_vector"], dtype=float).reshape(-1, 3))
    counts = np.array([len(r) for r in d["neighbordist"]], dtype=np.int64) \
        if not isinstance(d["neighbordist"], np.ndarray) else \
        np.full(natoms, d["neighbordist"].shape[1] if d["neighbordist"].ndim == 2 else 0)
    atom = np.repeat(np.arange(natoms), counts)
    if isinstance(d["neighbordist"], np.ndarray):
        dists = np.asarray(d["neighbordist"], dtype=float).reshape(-1)
    else:
        dists = np.fromiter(itertools.chain.from_iterable(d["neighbordist"]), dtype=float,
                            count=int(counts.sum()))
    if isinstance(d["diff"], np.ndarray):
        diffs = np.asarray(d["diff"], dtype=float).reshape(-1, 3)
    else:
        diffs = np.fromiter(
            itertools.chain.from_iterable(itertools.chain.from_iterable(d["diff"])),
            dtype=float, count=3 * int(counts.sum()),
        ).reshape(-1, 3)
    return atom, dists, diffs


def _ace_b_basis_nu1(A, lmax):
    """Compute nu=1 B-basis (isotropic density).
    
    B^{(1)}_{i,n} = A_{i,n,0,0} (l=0, m=0 component only)
    
    This captures the radial neighbor density.
    
    Parameters
    ----------
    A : ndarray
        A-basis coefficients from _ace_a_functions.
    lmax : int
        Maximum angular momentum (needed for indexing).
        
    Returns
    -------
    B1 : ndarray, shape (natoms, nmax)
        Nu=1 B-basis descriptors.
    """
    # l=0, m=0 is stored at A[i, n, 0, lmax]
    return np.real(A[:, :, 0, lmax])


def _ace_b_basis_nu2(A, nmax, lmax):
    """Compute nu=2 B-basis (power spectrum / SOAP-like).
    
    B^{(2)}_{i,n1,n2,l} = sum_{m=-l}^{l} A*_{i,n1,l,m} * A_{i,n2,l,m}
    
    This is equivalent to the SOAP power spectrum and captures
    2-body angular correlations.
    
    Parameters
    ----------
    A : ndarray
        A-basis coefficients.
    nmax : int
        Number of radial basis functions.
    lmax : int
        Maximum angular momentum.
        
    Returns
    -------
    B2 : ndarray, shape (natoms, n_descriptors)
        Nu=2 B-basis descriptors (flattened).
    """
    natoms = A.shape[0]
    descriptors = []
    
    # Use symmetry: only n1 <= n2
    for n1 in range(nmax):
        for n2 in range(n1, nmax):
            for l in range(lmax + 1):
                # Sum over m: sum_m A*_{n1,l,m} * A_{n2,l,m}
                B_desc = np.zeros(natoms)
                for m in range(-l, l + 1):
                    # Re(conj(a) * b) in real arithmetic: numpy's complex
                    # product can round differently from run to run
                    a, b = A[:, n1, l, m + lmax], A[:, n2, l, m + lmax]
                    B_desc += a.real * b.real + a.imag * b.imag
                descriptors.append(B_desc)
    
    return np.column_stack(descriptors) if descriptors else np.zeros((natoms, 0))


@functools.lru_cache(maxsize=None)
def _wigner_3j(j1, j2, j3, m1, m2, m3):
    """Wigner 3j symbol (j1 j2 j3; m1 m2 m3) for integer arguments (Racah formula)."""
    if m1 + m2 + m3 != 0:
        return 0.0
    if j3 < abs(j1 - j2) or j3 > j1 + j2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m3) > j3:
        return 0.0
    f = math.factorial
    delta = f(j1 + j2 - j3) * f(j1 - j2 + j3) * f(-j1 + j2 + j3) / f(j1 + j2 + j3 + 1)
    pref = math.sqrt(
        delta * f(j1 + m1) * f(j1 - m1) * f(j2 + m2) * f(j2 - m2) * f(j3 + m3) * f(j3 - m3)
    )
    tmin = max(0, j2 - j3 - m1, j1 - j3 + m2)
    tmax = min(j1 + j2 - j3, j1 - m1, j2 + m2)
    total = 0.0
    for t in range(tmin, tmax + 1):
        total += (-1) ** t / (
            f(t) * f(j1 + j2 - j3 - t) * f(j1 - m1 - t) * f(j2 + m2 - t)
            * f(j3 - j2 + m1 + t) * f(j3 - j1 - m2 + t)
        )
    return (-1) ** (j1 - j2 - m3) * pref * total


def _ace_b_basis_nu3(A, nmax, lmax):
    """Compute nu=3 B-basis (bispectrum-like triplet correlations).
    
    B^{(3)}_{n1 n2 n3 l1 l2 l3} = sum_{m1+m2+m3=0}
        (l1 l2 l3; m1 m2 m3) * A_{n1,l1,m1} * A_{n2,l2,m2} * A_{n3,l3,m3}
    
    where (l1 l2 l3; m1 m2 m3) is the Wigner 3j symbol, which couples the
    three A-functions to total angular momentum L=0 and thereby makes the
    descriptor rotationally invariant.
    
    This captures 3-body angular correlations.
    
    Parameters
    ----------
    A : ndarray
        A-basis coefficients.
    nmax : int
        Number of radial basis functions.
    lmax : int
        Maximum angular momentum.
        
    Returns
    -------
    B3 : ndarray, shape (natoms, n_descriptors)
        Nu=3 B-basis descriptors.
    """
    natoms = A.shape[0]
    descriptors = []
    
    # Limit combinations to keep computation tractable
    # Use n1 <= n2 <= n3 for symmetry
    for n1 in range(min(nmax, 3)):  # Limit radial indices
        for n2 in range(n1, min(nmax, 3)):
            for n3 in range(n2, min(nmax, 3)):
                for l1 in range(min(lmax + 1, 3)):  # Limit angular momentum
                    for l2 in range(min(lmax + 1, 3)):
                        # Triangle rule: |l1-l2| <= l3 <= l1+l2
                        l3_min = abs(l1 - l2)
                        l3_max = min(l1 + l2, lmax, 2)  # Also limit l3
                        for l3 in range(l3_min, l3_max + 1):
                            # Parity rule: l1 + l2 + l3 must be even
                            if (l1 + l2 + l3) % 2 != 0:
                                continue
                            
                            B_desc = np.zeros(natoms)
                            for m1 in range(-l1, l1 + 1):
                                for m2 in range(-l2, l2 + 1):
                                    m3 = -(m1 + m2)  # Enforce m1+m2+m3=0
                                    if abs(m3) > l3:
                                        continue
                                    w3j = _wigner_3j(l1, l2, l3, m1, m2, m3)
                                    if w3j == 0.0:
                                        continue
                                    
                                    # 3j-coupled product of three A-functions,
                                    # real part in real arithmetic (see nu=2)
                                    a = A[:, n1, l1, m1 + lmax]
                                    b = A[:, n2, l2, m2 + lmax]
                                    c = A[:, n3, l3, m3 + lmax]
                                    ab_re = a.real * b.real - a.imag * b.imag
                                    ab_im = a.real * b.imag + a.imag * b.real
                                    B_desc += w3j * (ab_re * c.real - ab_im * c.imag)
                            
                            # Always append so the descriptor count is a
                            # deterministic function of (nmax, lmax) — needed
                            # to compare descriptors across different structures.
                            descriptors.append(B_desc)
    
    return np.column_stack(descriptors) if descriptors else np.zeros((natoms, 0))


def ace(atoms: Atoms, nmax=4, lmax=4, nu_max=2, cutoff=None, normalize=True):
    """
    Compute Atomic Cluster Expansion (ACE) descriptors.
    
    ACE provides a systematic and complete expansion of atomic environments,
    with SOAP (nu=2) and bispectrum (nu=3) as special cases. The descriptors
    are rotationally, translationally, and permutationally invariant.
    
    The implementation follows Drautz (2019) and computes B-basis descriptors
    by coupling A-functions (single-particle basis) to form rotationally
    invariant combinations at each correlation order.
    
    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.
    nmax : int, default 4
        Number of radial basis functions. Higher values capture finer
        radial resolution but increase computation.
    lmax : int, default 4
        Maximum angular momentum quantum number. Higher values capture
        more angular detail. Typically 3-6 for ML potentials.
    nu_max : int, default 2
        Maximum correlation order:
        - nu=1: Radial density (neighbor count per shell)
        - nu=2: Pair correlations (SOAP power spectrum)
        - nu=3: Triplet correlations (bispectrum)
        Higher orders rapidly increase descriptor count.
    cutoff : float, optional
        Radial cutoff of the basis. If None, uses the cutoff from
        find_neighbors. Neighbors must have been computed with
        :func:`find_neighbors` beforehand; only neighbors inside the
        neighbor list contribute.
    normalize : bool, default True
        If True, divide each atom's descriptor vector by its L2 norm
        (per-atom normalisation across features).
        
    Returns
    -------
    dict
        Dictionary with keys:
        - 'nu1': ndarray (natoms, nmax) - radial density descriptors
        - 'nu2': ndarray (natoms, n2) - power spectrum (if nu_max >= 2)
        - 'nu3': ndarray (natoms, n3) - bispectrum-like (if nu_max >= 3)
        - 'full': ndarray (natoms, n_total) - concatenated descriptors
        
    Notes
    -----
    The descriptor count scales as:
    - nu=1: O(nmax)
    - nu=2: O(nmax^2 * lmax)
    - nu=3: O(nmax^3 * lmax^3) but limited for tractability
    
    References
    ----------
    .. [1] Drautz, R. (2019). "Atomic cluster expansion for accurate and 
           transferable interatomic potentials." Phys. Rev. B 99, 014104.
    .. [2] Dusson et al. (2022). "Atomic cluster expansion: Completeness,
           efficiency and stability." J. Comput. Phys.

    Examples
    --------
    >>> from ase.build import bulk
    >>> import pyscal3
    >>> atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    >>> pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
    >>> desc = pyscal3.ace(atoms, nmax=4, lmax=3, nu_max=2)
    >>> print(desc['full'].shape)
    >>> print("nu=2 descriptors:", desc['nu2'].shape[1])
    """
    d = _get_dict_with_neighbors(atoms)
    natoms = len(atoms)
    
    # Determine cutoff
    if cutoff is None:
        cutoffs = d.get("cutoff", [])
        cutoffs_arr = np.asarray(cutoffs)
        if cutoffs_arr.size > 0 and np.max(cutoffs_arr) > 0:
            cutoff = float(np.max(cutoffs_arr))
        else:
            cutoff = 5.0  # Default fallback
    
    # Compute A-functions (single-particle basis)
    A = _ace_a_functions(d, nmax, lmax, cutoff)
    
    result = {}
    all_descriptors = []
    
    # Nu=1: Radial density
    B1 = _ace_b_basis_nu1(A, lmax)
    result['nu1'] = B1
    all_descriptors.append(B1)
    
    # Nu=2: Power spectrum (SOAP-like)
    if nu_max >= 2:
        B2 = _ace_b_basis_nu2(A, nmax, lmax)
        result['nu2'] = B2
        all_descriptors.append(B2)
    
    # Nu=3: Triplet correlations (bispectrum-like)
    if nu_max >= 3:
        B3 = _ace_b_basis_nu3(A, nmax, lmax)
        result['nu3'] = B3
        all_descriptors.append(B3)
    
    # Concatenate all descriptors
    full = np.hstack(all_descriptors)
    
    # Normalize if requested
    if normalize:
        norms = np.linalg.norm(full, axis=1, keepdims=True)
        norms = np.where(norms > 1e-10, norms, 1.0)
        full = full / norms
        result['nu1'] = result['nu1'] / norms
        if 'nu2' in result:
            result['nu2'] = result['nu2'] / norms
        if 'nu3' in result:
            result['nu3'] = result['nu3'] / norms
    
    result['full'] = full
    
    # Store in atoms
    atoms.arrays["pyscal_ace"] = full
    atoms.info["pyscal_ace_params"] = {
        'nmax': nmax, 'lmax': lmax, 'nu_max': nu_max, 'cutoff': cutoff
    }
    
    return result

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
    pbc = np.asarray(reference.get_pbc(), dtype=bool)
    
    if any(pbc) and ref_cell.any():
        # Wrap the displaced atoms into the reference cell along the periodic
        # directions so that atoms that drifted several cells away (unwrapped
        # trajectories) are still matched to the right site.
        frac = np.linalg.solve(np.asarray(ref_cell).T, disp_pos.T).T
        frac[:, pbc] -= np.floor(frac[:, pbc])
        disp_pos = frac @ np.asarray(ref_cell)
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
        ref_types = np.asarray(reference.get_chemical_symbols())
        atom_types = np.asarray(atoms.get_chemical_symbols())
        # singly occupied sites whose atom has a different species than the site
        single = occupancy[site_index] == 1
        antisite_count = int(np.sum(atom_types[single] != ref_types[site_index[single]]))
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

