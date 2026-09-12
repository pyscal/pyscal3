"""
Bridge between ASE Atoms and the pyscal3 C++ extension.

The C++ functions expect a py::dict with specific keys.
This module converts ASE Atoms <-> that dict format,
and extracts box parameters (triclinic, rot, rotinv, boxdims)
from ASE's cell.

All pyscal-computed per-atom data is stored in atoms.arrays
(numpy arrays) or atoms.info (scalars/metadata).
"""

import numpy as np
from ase import Atoms


# ---- Keys managed by pyscal C++ ----
# Neighbor keys written by C++ neighbor routines
NEIGHBOR_KEYS = [
    "pyscal_neighbors",
    "pyscal_neighbordist",
    "pyscal_neighborweight",
    "pyscal_diff",
    "pyscal_r",
    "pyscal_phi",
    "pyscal_theta",
    "pyscal_cutoff",
    "pyscal_temp_neighbors",
    "pyscal_temp_neighbordist",
]

# We prefix pyscal keys in atoms.info to avoid clashes with ASE
_PREFIX = "pyscal_"


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

    The C++ code reads: positions, mask_1, mask_2, ghost.
    Numpy arrays are passed directly — pybind11 converts them
    to the C++ types automatically, avoiding expensive .tolist() calls.
    """
    n = len(atoms)
    d = {
        "positions": atoms.positions,  # numpy (n,3) — pybind11 casts directly
        "mask_1": [False] * n,
        "mask_2": [False] * n,
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
    since atoms.arrays requires uniform-length arrays.
    """
    skip_keys = {
        "positions",
        "mask_1",
        "mask_2",
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
                    atoms.arrays[store_key] = trimmed
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


def pad_atoms_for_neighbor_finding(atoms: Atoms, cutoff=None):
    """
    Build the atom dict for the C++ neighbor search, adding ghost atoms
    (periodic images created with ASE's ``repeat``) when the cell is too
    small for the requested search.

    Padding is applied when

    * the cell has fewer than 200 atoms or a perpendicular width below
      10 Angstrom (legacy rule, keeps the adaptive estimates stable), or
    * ``cutoff`` is given and some perpendicular width is not larger than
      ``2 * cutoff`` -- the minimum-image convention would otherwise miss
      neighbors beyond half the box.

    Ghost atoms are marked with ghost=True so results can be trimmed to
    the original atoms.

    Parameters
    ----------
    atoms : ase.Atoms
        Structure with a periodic cell.
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
    n = len(atoms)
    widths = perpendicular_widths(atoms.cell)
    if n == 0:
        raise ValueError("Cannot find neighbors of an empty Atoms object.")
    if np.any(widths <= 0):
        raise ValueError(
            "pyscal requires a periodic cell with three non-zero, non-coplanar "
            "cell vectors; got cell=%s. Set atoms.cell (and atoms.pbc) first."
            % np.asarray(atoms.cell).tolist()
        )

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

    if np.all(reps == 1):
        d = atoms_to_dict(atoms)
        return d, get_box_params(atoms), n

    # Create a repeated supercell using ASE
    supercell = atoms.repeat([int(r) for r in reps])
    nreal = n
    total = len(supercell)

    # Build the dict
    d = atoms_to_dict(supercell)

    # Mark ghost atoms
    d["ghost"] = [False] * nreal + [True] * (total - nreal)

    # Head array: maps each atom to its real-atom index
    d["head"] = [i % nreal for i in range(total)]
    d["nreal"] = nreal

    return d, get_box_params(supercell), nreal
