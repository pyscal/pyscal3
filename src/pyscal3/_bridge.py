"""
Bridge between ASE Atoms and the pyscal3 C++ extension.

The C++ functions expect a py::dict with specific keys.
This module converts ASE Atoms <-> that dict format,
and extracts box parameters (triclinic, rot, rotinv, boxdims)
from ASE's cell.

All pyscal-computed per-atom data is stored in atoms.arrays
(numpy arrays) or atoms.info (scalars/metadata).
"""

import contextlib
import gc
import itertools

import numpy as np
from ase import Atoms


# We prefix pyscal keys in atoms.arrays / atoms.info to avoid clashes with ASE
_PREFIX = "pyscal_"

# ---- Keys managed by pyscal C++ ----
# Everything derived from a neighbor search (all methods, incl. Voronoi).
# These are removed before a new search so that no stale data survives a
# change of neighbor method.
NEIGHBOR_DERIVED_KEYS = [
    "neighbors",
    "neighbordist",
    "neighborweight",
    "diff",
    "r",
    "phi",
    "theta",
    "cutoff",
    "temp_neighbors",
    "temp_neighbordist",
    "voronoi_volume",
    "face_vertices",
    "face_perimeters",
    "vertex_vectors",
    "vertex_numbers",
    "vertex_positions",
    "unique_vertices",
    "neighbors_found",
    "neighbor_method",
]
NEIGHBOR_KEYS = [_PREFIX + k for k in NEIGHBOR_DERIVED_KEYS]


@contextlib.contextmanager
def gc_paused():
    """Pause the cyclic garbage collector while building many small lists.

    Per-atom neighbor rows are millions of list objects that cannot form
    reference cycles. Allocating them triggers repeated collections, which
    scan every list created so far and can take most of the time of a
    neighbor search. The previous state of the collector is restored.
    """
    enabled = gc.isenabled()
    gc.disable()
    try:
        yield
    finally:
        if enabled:
            gc.enable()


def clear_neighbor_data(atoms: Atoms):
    """Remove all neighbor-derived pyscal keys from ``atoms``."""
    for key in NEIGHBOR_KEYS:
        atoms.arrays.pop(key, None)
        atoms.info.pop(key, None)


def neighbor_arrays(atoms: Atoms, *keys):
    """Flat arrays of stored per-bond neighbor data, for the C++ routines.

    Returns ``(offsets, values)``: ``offsets`` has ``len(atoms) + 1`` entries,
    and the bonds of atom ``i`` are ``offsets[i]:offsets[i + 1]`` in each
    array of ``values``, which maps every requested key (without the
    ``pyscal_`` prefix) to a flat array. ``neighbors`` gives int64 indices,
    ``diff`` an (m, 3) array, the other keys float64 arrays.

    The arrays are rebuilt from ``atoms.arrays`` (rows of equal length, which
    only need a reshape) or ``atoms.info`` (lists of lists), so they always
    match what is stored.
    """
    ensure_neighbors(atoms)
    n = len(atoms)
    counts = None
    values = {}
    for key in keys:
        store = _PREFIX + key
        if store in atoms.arrays:
            rows = atoms.arrays[store]
        elif store in atoms.info:
            rows = atoms.info[store]
        else:
            raise ValueError("No %s stored; call pyscal3.find_neighbors first." % store)
        dtype = np.int64 if key == "neighbors" else np.float64
        if isinstance(rows, np.ndarray):
            k = rows.shape[1] if rows.ndim >= 2 else 0
            c = np.full(n, k, dtype=np.int64)
            flat = rows.reshape((n * k, 3) if key == "diff" else (n * k,))
        else:
            c = np.fromiter(map(len, rows), dtype=np.int64, count=n)
            m = int(c.sum())
            if key == "diff":
                flat = np.fromiter(
                    itertools.chain.from_iterable(itertools.chain.from_iterable(rows)),
                    dtype=np.float64, count=3 * m,
                ).reshape(m, 3)
            else:
                flat = np.fromiter(itertools.chain.from_iterable(rows), dtype=dtype, count=m)
        if counts is None:
            counts = c
        elif not np.array_equal(counts, c):
            raise ValueError("Stored neighbor keys have different lengths: %s" % (keys,))
        values[key] = np.ascontiguousarray(flat, dtype=dtype)
    if counts is None:
        counts = np.zeros(n, dtype=np.int64)
    offsets = np.zeros(n + 1, dtype=np.int64)
    np.cumsum(counts, out=offsets[1:])
    return offsets, values


def rows_from_flat(offsets, flat):
    """Per-atom rows of a flat per-bond array, in the stored-key format.

    Rows of equal length give an (n, k) array (an (n, 0) float array if every
    row is empty); rows of different length give a list of lists.
    """
    offsets = np.asarray(offsets)
    n = len(offsets) - 1
    counts = np.diff(offsets)
    k = int(counts[0]) if n else 0
    if (counts == k).all():
        return np.zeros((n, 0)) if k == 0 else np.asarray(flat).reshape(n, k)
    with gc_paused():
        values = np.asarray(flat).tolist()
        return [values[a:b] for a, b in zip(offsets[:-1].tolist(), offsets[1:].tolist())]


def stored_per_atom(atoms: Atoms, key):
    """A per-atom pyscal array (``pyscal_<key>``) from atoms.arrays or atoms.info."""
    store = _PREFIX + key
    if store in atoms.arrays:
        return np.asarray(atoms.arrays[store])
    return np.asarray(atoms.info[store])


def get_box_params(atoms: Atoms):
    """
    Extract box parameters from ASE cell for the C++ functions.

    Returns
    -------
    triclinic : int
        0 for orthorhombic, 1 for triclinic
    rot : list of list of float
        Cell vectors transposed (rotation matrix)
    rotinv : list of list of float
        Inverse of rot
    boxdims : list of float
        Box side lengths [Lx, Ly, Lz]
    """
    cell = np.array(atoms.cell)

    # Check if the cell is non-orthorhombic (triclinic).
    # The C++ orthorhombic code path assumes box edges are aligned with
    # the Cartesian x/y/z axes, so the cell matrix must be diagonal.
    # A cell with mutually perpendicular but *rotated* vectors (e.g. a
    # cubic cell after rigid rotation) has zero dot-products between edges
    # but non-zero off-diagonal elements and must use the triclinic path.
    off_diag = cell - np.diag(np.diag(cell))

    triclinic = 0
    rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    rotinv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

    if np.max(np.abs(off_diag)) > 1e-6:
        triclinic = 1
        rot = cell.T.tolist()
        rotinv = np.linalg.inv(cell.T).tolist()

    boxdims = [np.linalg.norm(cell[i]) for i in range(3)]

    return triclinic, rot, rotinv, boxdims


def atoms_to_dict(atoms: Atoms) -> dict:
    """
    Convert ASE Atoms to the dict format expected by pyscal C++ functions.

    The C++ code reads: positions, ghost.
    Numpy arrays are passed directly — pybind11 converts them
    to the C++ types automatically, avoiding expensive .tolist() calls.
    """
    n = len(atoms)
    d = {
        "positions": atoms.positions,  # numpy (n,3) — pybind11 casts directly
        "ghost": [False] * n,
        "types": atoms.get_atomic_numbers(),  # numpy 1-D
    }

    # Copy any existing pyscal data from atoms.arrays (keep as numpy)
    for key in atoms.arrays:
        if key.startswith(_PREFIX):
            d[key[len(_PREFIX) :]] = atoms.arrays[key]

    # Copy pyscal info keys (may be ragged lists — keep as-is)
    for key in atoms.info:
        if key.startswith(_PREFIX):
            d[key[len(_PREFIX) :]] = atoms.info[key]

    return d


def dict_to_atoms(d: dict, atoms: Atoms, nreal=None):
    """
    Write C++ results from dict back into ASE Atoms.arrays and atoms.info.

    If nreal is given, only the first nreal entries are written
    (used when ghost atoms were added for neighbor finding).

    Handles ragged arrays (neighbors, etc.) by storing in atoms.info
    since atoms.arrays requires uniform-length arrays. Arrays with three
    or more dimensions (e.g. the neighbor vectors ``diff``) also go to
    atoms.info because ASE file writers only support per-atom scalars and
    vectors.
    """
    skip_keys = {
        "positions",
        "ghost",
        "types",
        "ids",
        "head",
        "nreal",
    }
    n = len(atoms)
    if nreal is None:
        nreal = n

    # Build index remap if ghost atoms were used
    head = d.get("head")
    _NEIGHBOR_INDEX_KEYS = {"neighbors", "temp_neighbors"}

    for key, val in d.items():
        if key in skip_keys:
            continue

        store_key = _PREFIX + key

        # Fast path: detect ragged (list-of-lists with varying sub-lengths)
        # early to skip the expensive np.asarray() probe that creates slow
        # object arrays.  Uniform list-of-lists (same sub-length) go through
        # the numpy path so they can be stored as 2-D arrays.
        if isinstance(val, list) and len(val) > 1 and isinstance(val[0], (list, tuple)):
            if len(val[0]) != len(val[1]):
                # Definitely ragged — store directly in info
                trimmed = val[:nreal]
                if key in _NEIGHBOR_INDEX_KEYS and head is not None:
                    trimmed = [
                        [head[j] if j < len(head) else j for j in row]
                        for row in trimmed
                    ]
                atoms.info[store_key] = trimmed
                continue

        # Try to store as atoms.arrays (requires same-shape numpy array)
        try:
            arr = np.asarray(val)
            if arr.ndim >= 1 and arr.dtype.kind != "O":
                # Trim to real atoms
                trimmed = arr[:nreal]
                # Remap neighbor indices that reference ghost atoms
                if key in _NEIGHBOR_INDEX_KEYS and head is not None:
                    head_arr = np.array(head)
                    trimmed = head_arr[trimmed]
                if len(trimmed) == n:
                    if arr.ndim <= 2:
                        atoms.arrays[store_key] = trimmed
                    else:
                        atoms.info[store_key] = trimmed
                    continue
        except (ValueError, TypeError, IndexError):
            pass

        # Ragged data (lists of lists) — trim and remap indices
        if isinstance(val, list) and len(val) >= nreal:
            trimmed = val[:nreal]
            # Remap neighbor indices if they reference ghost atoms
            if key in _NEIGHBOR_INDEX_KEYS and head is not None:
                trimmed = [
                    [head[j] if j < len(head) else j for j in row] for row in trimmed
                ]
            atoms.info[store_key] = trimmed
        else:
            atoms.info[store_key] = val


def ensure_neighbors(atoms: Atoms):
    """Check that neighbors have been computed for these atoms."""
    if (
        _PREFIX + "neighbors" not in atoms.info
        and _PREFIX + "neighbors" not in atoms.arrays
    ):
        raise ValueError(
            "Neighbors have not been computed. "
            "Call pyscal3.find_neighbors(atoms, ...) first."
        )


def create_attribute(d: dict, key: str, fill_with=0):
    """Create a new key in the dict filled with a default value."""
    n = len(d["positions"])
    if isinstance(fill_with, (int, float)):
        d[key] = [fill_with] * n
    else:
        d[key] = [fill_with] * n


# ---------------------------------------------------------------------------
# Ghost atom padding for small cells
# ---------------------------------------------------------------------------
_MIN_ATOMS = 200
_MIN_BOX_SIDE = 10.0  # Angstroms


def perpendicular_widths(cell):
    """Perpendicular width of the cell along each cell vector.

    The width along vector *i* is the volume divided by the area of the
    face spanned by the other two vectors.  For an orthogonal cell these
    are simply the box lengths.  The minimum-image convention used by the
    C++ routines is exact only for distances below half of the smallest
    width, so this is the quantity that decides how much padding is needed.
    """
    cell = np.asarray(cell, dtype=float)
    vol = abs(np.linalg.det(cell))
    widths = np.zeros(3)
    for i in range(3):
        cross = np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])
        area = np.linalg.norm(cross)
        widths[i] = vol / area if area > 0 else 0.0
    return widths


def guess_cutoff(atoms: Atoms, prefactor):
    """Candidate-search radius used by the adaptive/SANN/number methods.

    Mirrors the C++ estimate ``prefactor * (V / N)^(1/3)``.
    """
    return prefactor * (abs(np.linalg.det(np.asarray(atoms.cell))) / len(atoms)) ** (1.0 / 3.0)


_NONPERIODIC_PAD = 10.0  # Angstrom, vacuum added when no search radius is known


def periodic_directions(atoms: Atoms):
    """Directions that are periodic and have a non-zero cell vector."""
    lengths = np.linalg.norm(np.array(atoms.cell, dtype=float), axis=1)
    return np.array(atoms.pbc, dtype=bool) & (lengths > 0)


def effective_periodic_cell(atoms: Atoms, pad):
    """Cell to use for the (always periodic) C++ routines.

    Directions that are periodic and have a non-zero cell vector are kept.
    Every other direction (``pbc`` False, or a zero cell vector as for a
    molecule read from an XYZ file) is replaced by a vector orthogonal to
    the periodic ones whose length is the extent of the atoms along it plus
    ``2 * pad``. Periodic images along such a direction are then at least
    ``2 * pad`` apart, so no image can be found within a search radius of
    ``pad``.

    Returns
    -------
    cell : ndarray (3, 3)
    periodic : ndarray of bool
        Which directions were genuinely periodic.
    """
    cell = np.array(atoms.cell, dtype=float)
    periodic = periodic_directions(atoms)

    if periodic.all():
        if abs(np.linalg.det(cell)) <= 0:
            raise ValueError(
                "pyscal requires three non-coplanar cell vectors; got cell=%s"
                % cell.tolist()
            )
        return cell, periodic

    kept = cell[periodic]
    if len(kept) > 0 and np.linalg.matrix_rank(kept) < len(kept):
        raise ValueError(
            "The periodic cell vectors are linearly dependent: %s" % cell.tolist()
        )
    # orthonormal complement of the periodic directions
    if len(kept) == 0:
        complement = np.eye(3)
    else:
        _, _, vt = np.linalg.svd(kept)
        complement = vt[len(kept):]

    positions = atoms.positions
    new_cell = cell.copy()
    for j, i in enumerate(np.where(~periodic)[0]):
        direction = complement[j]
        proj = positions @ direction
        extent = float(proj.max() - proj.min()) if len(proj) else 0.0
        new_cell[i] = direction * (extent + 2.0 * pad)
    return new_cell, periodic


def padded_supercell(atoms: Atoms, cutoff=None):
    """
    The periodic cell used by the C++ routines that need one, replicated with
    ASE's ``repeat`` when it is too small for the requested search.

    Non-periodic directions and zero cell vectors are handled by
    :func:`effective_periodic_cell`, which adds enough vacuum that periodic
    images cannot be found within the search radius.

    The cell is replicated along periodic directions when

    * the cell has fewer than 200 atoms or a perpendicular width below
      10 Angstrom (legacy rule, keeps the adaptive estimates stable), or
    * ``cutoff`` is given and some perpendicular width is not larger than
      ``2 * cutoff`` -- the minimum-image convention would otherwise miss
      neighbors beyond half the box.

    Parameters
    ----------
    atoms : ase.Atoms
        The structure.
    cutoff : float, optional
        Largest distance the neighbor search has to resolve.

    Returns
    -------
    ase.Atoms
        Fully periodic cell whose first ``len(atoms)`` atoms are the original
        ones. It is ``atoms`` itself when nothing had to change.
    """
    n = len(atoms)
    if n == 0:
        raise ValueError("Cannot find neighbors of an empty Atoms object.")

    pad = cutoff if (cutoff is not None and cutoff > 0) else _NONPERIODIC_PAD
    work_cell, periodic = effective_periodic_cell(atoms, pad)
    if periodic.all():
        work = atoms
    else:
        work = atoms.copy()
        work.set_cell(work_cell)
        work.set_pbc(True)

    widths = perpendicular_widths(work_cell)
    reps = np.ones(3, dtype=int)

    if n < _MIN_ATOMS:
        needed = max(int(np.ceil((_MIN_ATOMS / n) ** (1.0 / 3.0))), 2)
        reps[:] = needed
        for i in range(3):
            if widths[i] * reps[i] < _MIN_BOX_SIDE:
                reps[i] = max(reps[i], int(np.ceil(_MIN_BOX_SIDE / widths[i])))

    if cutoff is not None and cutoff > 0:
        # minimum image is exact only if every width exceeds 2 * cutoff
        need = 2.0 * cutoff * (1.0 + 1e-6)
        for i in range(3):
            if widths[i] * reps[i] <= need:
                reps[i] = max(reps[i], int(np.ceil(need / widths[i])))
            if widths[i] * reps[i] <= need:
                reps[i] += 1

    # never replicate along non-periodic directions
    reps[~periodic] = 1

    if np.all(reps == 1):
        return work
    return work.repeat([int(r) for r in reps])


def pad_atoms_for_neighbor_finding(atoms: Atoms, cutoff=None):
    """
    Build the atom dict for the C++ routines that work on a periodic cell
    (Voronoi, CNA), adding ghost atoms when the cell is too small.

    The padded cell is the one of :func:`padded_supercell`. Ghost atoms are
    marked with ghost=True so results can be trimmed to the original atoms.

    Parameters
    ----------
    atoms : ase.Atoms
        The structure.
    cutoff : float, optional
        Largest distance the neighbor search has to resolve.

    Returns
    -------
    d : dict
        Atom dict (possibly with ghost atoms).
    box_params : tuple
        (triclinic, rot, rotinv, boxdims) for the (possibly padded) box.
    nreal : int
        Number of real (non-ghost) atoms.
    """
    supercell = padded_supercell(atoms, cutoff=cutoff)
    nreal = len(atoms)
    total = len(supercell)
    d = atoms_to_dict(supercell)
    if total == nreal:
        return d, get_box_params(supercell), nreal

    # Mark ghost atoms
    d["ghost"] = [False] * nreal + [True] * (total - nreal)

    # Head array: maps each atom to its real-atom index
    d["head"] = [i % nreal for i in range(total)]
    d["nreal"] = nreal

    return d, get_box_params(supercell), nreal
