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
from scipy.special import sph_harm_y

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
# Wigner W_l Parameter (Third-order rotational invariant)
# ---------------------------------------------------------------------------


def wigner_w_parameter(atoms: Atoms, l, averaged=False, normalized=True):
    """
    Calculate the third-order Steinhardt invariant W_l.

    W_l is the third-order rotational invariant of the bond-orientational
    order parameters, constructed by contracting q_lm with Wigner 3j
    symbols. It distinguishes crystal structures that have similar q_l
    values (e.g., FCC vs HCP via the sign of W_6).

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
        If True (default), return the normalized hat{W}_l =
        W_l / (sum |q_lm|^2)^(3/2). If False, return raw W_l.

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
    if isinstance(l, int):
        ll = [l]
    else:
        ll = list(l)

    d = _get_dict_with_neighbors(atoms)

    # Ensure q_lm are computed first (W_l requires them)
    for val in ll:
        pc.calculate_q_single(d, val)

    if averaged:
        for val in ll:
            pc.calculate_aw_single(d, val)
        if normalized:
            result_keys = ["avg_what%d" % v for v in ll]
        else:
            result_keys = ["avg_w%d" % v for v in ll]
    else:
        for val in ll:
            pc.calculate_w_single(d, val)
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
    dists = atoms.arrays["pyscal_neighbordist"]
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

    Uses the chi-parameter histogram (angular signature) to assign each
    atom a structure type via a decision tree.  The chi histogram bins
    cosines of all pairwise neighbor angles into 9 ranges; the counts in
    specific bins discriminate FCC, HCP, BCC, and icosahedral (ICO)
    environments.

    The algorithm follows Ackland & Jones, *Phys. Rev. B* **73**, 054104
    (2006), adapted for the 9-bin chi scheme used by pyscal3.
    Parameters
    ----------
    atoms : ase.Atoms        Structure with neighbors already computed (via
        :func:`pyscal3.find_neighbors`).

    Returns
    -------
    labels : numpy.ndarray of int, shape (natoms,)
        Per-atom structure label:

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
        Human-readable name for each atom (``"fcc"``, ``"hcp"``, etc.).
        Also stored in ``atoms.arrays["pyscal_structure"]``.

    Notes
    -----
    Chi-parameter bin edges (9 bins, cosine of the angle):

    =====  =========================
    Bin    Cosine range
    =====  =========================
    χ₀     [−1.000, −0.945)  ~180°
    χ₁     [−0.945, −0.915)  ~160°
    χ₂     [−0.915, −0.755)  ~139°
    χ₃     [−0.755, −0.705)  ~135°
    χ₄     [−0.705, −0.195)
    χ₅     [−0.195, +0.195)  ~90°
    χ₆     [+0.195, +0.245)
    χ₇     [+0.245, +0.795)  ~60°
    χ₈     [+0.795, +1.000]  ~0°
    =====  =========================

    Decision tree (ideal counts for perfect structures):

    * **BCC** (14 neighbours): χ₀ = 7, χ₅ = 12, χ₇ = 36
    * **FCC** (12 neighbours): χ₀ = 6, χ₅ = 12, χ₇ = 24, χ₁₊₂₊₃ = 0
    * **HCP** (12 neighbours): χ₀ = 3, χ₂ = 6, χ₇ = 24
    * **ICO** (12 neighbours): χ₀ = 6, χ₅ = 0, χ₇ = 30

    References
    ----------
    G. J. Ackland and A. P. Jones, "Applications of local crystal
    structure measures in experiment and simulation",
    *Phys. Rev. B* **73**, 054104 (2006).
    `doi:10.1103/PhysRevB.73.054104
    <https://doi.org/10.1103/PhysRevB.73.054104>`__
    """
    chi = chi_params(atoms)  # (N, 9) int array

    n = len(atoms)
    labels = np.zeros(n, dtype=int)

    # --- BCC: 7 or more antiparallel (180°) pairs ---
    bcc_mask = chi[:, 0] >= 7
    labels[bcc_mask] = ACKLAND_BCC

    # --- FCC or ICO: 6 antiparallel pairs, no intermediate angles ---
    remaining = labels == 0
    fcc_or_ico = remaining & (chi[:, 0] >= 5) & (
        chi[:, 1] + chi[:, 2] + chi[:, 3] == 0
    )
    # ICO has no ~90° pairs (χ₅ = 0); FCC has χ₅ > 0
    labels[fcc_or_ico & (chi[:, 5] == 0)] = ACKLAND_ICO
    labels[fcc_or_ico & (chi[:, 5] > 0)] = ACKLAND_FCC

    # --- HCP: has angles in the ~139° range (χ₂ > 0) ---
    remaining = labels == 0
    hcp_mask = remaining & (chi[:, 2] > 0)
    labels[hcp_mask] = ACKLAND_HCP

    # Store results
    names = [_ACKLAND_NAMES[l] for l in labels]
    atoms.arrays["pyscal_ackland_label"] = labels
    atoms.arrays["pyscal_structure"] = np.array(names)

    return labels, names

# Deformation Descriptors (require reference configuration)
# ---------------------------------------------------------------------------

def _match_neighbor_indices(atoms_cur, atoms_ref):
    """
    Match atoms across deformed/reference configs and return neighbor mapping.

    Returns dict mapping atom index → list of (neighbor_index_cur, neighbor_index_ref)
    for atoms that appear in both neighbor lists.
    """
    # Get current neighbor data
    neighbors_cur = atoms_cur.arrays["pyscal_neighbors"]
    neighbors_ref = atoms_ref.arrays["pyscal_neighbors"]
    diff_cur = atoms_cur.arrays["pyscal_diff"]
    diff_ref = atoms_ref.arrays["pyscal_diff"]

    n = len(atoms_cur)
    mapping = {}

    for i in range(n):
        nbrs_cur = neighbors_cur[i]
        nbrs_ref = neighbors_ref[i]
        valid_cur = nbrs_cur[nbrs_cur >= 0]
        valid_ref = nbrs_ref[nbrs_ref >= 0]
        # Find common neighbors
        common = set(valid_cur) & set(valid_ref)
        pairs = []
        for j in common:
            j_cur_idx = np.where(nbrs_cur == j)[0][0]
            j_ref_idx = np.where(nbrs_ref == j)[0][0]
            pairs.append((j_cur_idx, j_ref_idx))
        mapping[i] = pairs

    return mapping, diff_cur, diff_ref


def atomic_strain(atoms: Atoms, reference: Atoms):
    """
    Calculate the atomic (Green-Lagrange) strain tensor for each atom.

    The local deformation gradient F is computed by minimizing the
    squared difference between reference and deformed neighbor vectors.
    The Lagrangian strain is then E = (F^T F - I) / 2.
    Parameters
    ----------
    atoms : ase.Atoms        Deformed configuration with neighbors computed.
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
    ensure_neighbors(atoms)
    ensure_neighbors(reference)

    mapping, diff_cur, diff_ref = _match_neighbor_indices(atoms, reference)
    n = len(atoms)
    strain = np.zeros((n, 3, 3))

    for i in range(n):
        pairs = mapping[i]
        if len(pairs) < 3:
            strain[i] = np.nan
            continue

        # Build matrices: rows = neighbor displacement vectors
        X = np.array([diff_ref[i, j_ref] for (_, j_ref) in pairs])  # reference
        Y = np.array([diff_cur[i, j_cur] for (j_cur, _) in pairs])  # current

        # Deformation gradient F via least squares: Y = X @ F^T
        # => F^T = (X^T X)^-1 X^T Y
        XtX = X.T @ X
        try:
            XtX_inv = np.linalg.inv(XtX)
        except np.linalg.LinAlgError:
            strain[i] = np.nan
            continue
        F = (XtX_inv @ X.T @ Y).T

        # Green-Lagrange strain E = (F^T F - I) / 2
        C = F.T @ F
        E = 0.5 * (C - np.eye(3))
        strain[i] = E

    atoms.arrays["pyscal_strain"] = strain
    return strain


def von_mises_strain(atoms: Atoms, reference: Atoms):
    """
    Compute the von Mises shear strain invariant from the atomic strain.

    .. math::

        \\eta^{\\text{Mises}} = \\sqrt{
            \\frac{1}{2} \\left[
                (E_{xx}-E_{yy})^2 + (E_{yy}-E_{zz})^2 + (E_{zz}-E_{xx})^2
            \\right]
            + E_{xy}^2 + E_{yz}^2 + E_{xz}^2
        }

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
    n = len(atoms)
    vm = np.zeros(n)

    for i in range(n):
        if np.any(np.isnan(E[i])):
            vm[i] = np.nan
            continue
        exx, eyy, ezz = E[i, 0, 0], E[i, 1, 1], E[i, 2, 2]
        exy, eyz, exz = E[i, 0, 1], E[i, 1, 2], E[i, 0, 2]
        vm[i] = np.sqrt(
            0.5 * ((exx - eyy)**2 + (eyy - ezz)**2 + (ezz - exx)**2)
            + exy**2 + eyz**2 + exz**2
        )

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
    ensure_neighbors(atoms)
    ensure_neighbors(reference)

    mapping, diff_cur, diff_ref = _match_neighbor_indices(atoms, reference)
    n = len(atoms)
    d2 = np.zeros(n)

    for i in range(n):
        pairs = mapping[i]
        if len(pairs) < 3:
            d2[i] = np.nan
            continue

        X = np.array([diff_ref[i, j_ref] for (_, j_ref) in pairs])
        Y = np.array([diff_cur[i, j_cur] for (j_cur, _) in pairs])

        XtX = X.T @ X
        try:
            XtX_inv = np.linalg.inv(XtX)
        except np.linalg.LinAlgError:
            d2[i] = np.nan
            continue
        F = (XtX_inv @ X.T @ Y).T

        # Predicted positions from affine: Y_pred = X @ F^T
        Y_pred = X @ F.T
        residuals = Y - Y_pred
        d2[i] = np.mean(np.sum(residuals**2, axis=1))

    atoms.arrays["pyscal_d2min"] = d2
    return d2


def slip_vector(atoms: Atoms, reference: Atoms):
    """
    Compute the slip vector for each atom.

    The slip vector s_i = (1/N_s) sum_j (dX_ij - dx_ij) is the average
    difference between reference and current displacement vectors,
    without fitting an affine transformation.

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
    ensure_neighbors(atoms)
    ensure_neighbors(reference)

    mapping, diff_cur, diff_ref = _match_neighbor_indices(atoms, reference)
    n = len(atoms)
    slip = np.zeros((n, 3))

    for i in range(n):
        pairs = mapping[i]
        if len(pairs) == 0:
            slip[i] = np.nan
            continue

        deltas = []
        for j_cur, j_ref in pairs:
            d_cur = diff_cur[i, j_cur]
            d_ref = diff_ref[i, j_ref]
            deltas.append(d_cur - d_ref)

        slip[i] = np.mean(deltas, axis=0)

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
# Coordination Variants
# ---------------------------------------------------------------------------

def effective_coordination_number(atoms: Atoms):
    """
    Calculate the effective coordination number (ECoN).

    ECoN weights each neighbor by a continuous function of distance
    relative to the average neighbor distance, converging iteratively:

    .. math::

        \\mathrm{ECoN}_i = \\sum_j \\exp\\!\\left[
            1 - \\left(\\frac{d_{ij}}{d_{\\mathrm{av},i}}\\right)^6
        \\right]

    where :math:`d_{\\mathrm{av},i}` is the weighted average distance
    that is itself computed from the weights (self-consistent iteration).

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with neighbors already computed.

    Returns
    -------
    numpy.ndarray, shape (natoms,)
        Per-atom ECoN values.  Also stored in
        ``atoms.arrays["pyscal_econ"]``.

    References
    ----------
    R. Hoppe, "Effective coordination numbers (ECoN) and mean fictive
    ionic radii (MEFIR)", *Z. Kristallogr.* **150**, 23 (1979).
    """
    ensure_neighbors(atoms)
    dists = atoms.arrays["pyscal_neighbordist"]   # (N, max_nn)
    mask = dists > 0  # valid neighbors

    n = len(atoms)
    econ = np.zeros(n, dtype=float)

    for _ in range(10):  # iterate to self-consistency
        # weighted average distance
        weights = np.where(mask, np.exp(1.0 - (dists / np.where(
            econ[:, None] > 0,
            _weighted_avg_dist(dists, mask, econ),
            np.where(mask, dists, np.inf).min(axis=1, keepdims=True)
        ))**6), 0.0)
        econ_new = weights.sum(axis=1)
        if np.allclose(econ, econ_new, atol=1e-8):
            break
        econ = econ_new

    atoms.arrays["pyscal_econ"] = econ
    return econ


def _weighted_avg_dist(dists, mask, econ):
    """Weighted average distance for ECoN iteration."""
    weights = np.where(mask, np.exp(1.0 - (dists / np.where(
        mask, dists, np.inf).min(axis=1, keepdims=True))**6), 0.0)
    wsum = weights.sum(axis=1, keepdims=True)
    wsum = np.where(wsum > 0, wsum, 1.0)
    return (weights * dists).sum(axis=1, keepdims=True) / wsum


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
    dists = atoms.arrays["pyscal_neighbordist"]
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

    neighbors = atoms.arrays["pyscal_neighbors"]
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
    dists = atoms.arrays["pyscal_neighbordist"]
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
    
    for i in range(natoms):
        neighbors_i = d["neighbors"][i]
        if not hasattr(neighbors_i, '__len__') or len(neighbors_i) == 0:
            continue
            
        dists = d["neighbordist"][i]
        diffs = d["diff"][i]
        
        for j_idx in range(len(neighbors_i)):
            rij = dists[j_idx]
            if rij < 1e-10 or rij >= cutoff:
                continue
            
            vec = np.array(diffs[j_idx])
            # Spherical coordinates
            # theta = polar angle from z axis
            # phi = azimuthal angle in xy plane
            theta = np.arccos(np.clip(vec[2] / rij, -1, 1))
            phi = np.arctan2(vec[1], vec[0])
            
            for n in range(nmax):
                R_n = _ace_radial_basis(n, rij, cutoff)
                for l in range(lmax + 1):
                    for m in range(-l, l + 1):
                        # scipy sph_harm_y(l, m, theta, phi) uses physics convention
                        Y_lm = sph_harm_y(l, m, theta, phi)
                        A[i, n, l, m + lmax] += R_n * Y_lm
    
    return A


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
                    B_desc += np.real(
                        np.conj(A[:, n1, l, m + lmax]) * A[:, n2, l, m + lmax]
                    )
                descriptors.append(B_desc)
    
    return np.column_stack(descriptors) if descriptors else np.zeros((natoms, 0))


def _ace_b_basis_nu3(A, nmax, lmax):
    """Compute nu=3 B-basis (bispectrum-like triplet correlations).
    
    B^{(3)} = sum_{m1,m2,m3} C_{l1,l2,l3}^{m1,m2,m3} * A_{n1,l1,m1} * A_{n2,l2,m2} * A_{n3,l3,m3}
    
    where the coupling coefficient ensures rotational invariance (total L=0).
    For simplicity, we use the constraint m1 + m2 + m3 = 0 with equal weights.
    
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
                                    
                                    # Product of three A-functions
                                    prod = (A[:, n1, l1, m1 + lmax] *
                                            A[:, n2, l2, m2 + lmax] *
                                            A[:, n3, l3, m3 + lmax])
                                    B_desc += np.real(prod)
                            
                            if np.any(np.abs(B_desc) > 1e-15):
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
        Neighbor cutoff radius. If None, uses the cutoff from
        find_neighbors.
    normalize : bool, default True
        If True, normalize descriptors to unit norm per atom.
        
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

