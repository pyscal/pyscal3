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
    "pyscal_neighbors", "pyscal_neighbordist", "pyscal_neighborweight",
    "pyscal_diff", "pyscal_r", "pyscal_phi", "pyscal_theta", "pyscal_cutoff",
    "pyscal_temp_neighbors", "pyscal_temp_neighbordist",
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
    
    # Check if triclinic: compute sum of dot products of adjacent cell vector pairs
    summ = 0.0
    for i in range(3):
        box1 = cell[i - 1]
        box2 = cell[i]
        n1 = np.linalg.norm(box1)
        n2 = np.linalg.norm(box2)
        if n1 > 0 and n2 > 0:
            summ += np.dot(box1, box2) / (n1 * n2)
    
    triclinic = 0
    rot = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    rotinv = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    
    if abs(summ) > 1e-6:
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
        "positions": atoms.positions,     # numpy (n,3) — pybind11 casts directly
        "mask_1": [False] * n,
        "mask_2": [False] * n,
        "ghost": [False] * n,
        "types": atoms.get_atomic_numbers(),  # numpy 1-D
    }
    
    # Copy any existing pyscal data from atoms.arrays (keep as numpy)
    for key in atoms.arrays:
        if key.startswith(_PREFIX):
            d[key[len(_PREFIX):]] = atoms.arrays[key]
    
    # Copy pyscal info keys (may be ragged lists — keep as-is)
    for key in atoms.info:
        if key.startswith(_PREFIX):
            d[key[len(_PREFIX):]] = atoms.info[key]
    
    return d


def dict_to_atoms(d: dict, atoms: Atoms, nreal=None):
    """
    Write C++ results from dict back into ASE Atoms.arrays and atoms.info.
    
    If nreal is given, only the first nreal entries are written 
    (used when ghost atoms were added for neighbor finding).
    
    Handles ragged arrays (neighbors, etc.) by storing in atoms.info 
    since atoms.arrays requires uniform-length arrays.
    """
    skip_keys = {"positions", "mask_1", "mask_2", "ghost", "types",
                 "ids", "head", "nreal"}
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
        
        # Try to store as atoms.arrays (requires same-shape numpy array)
        try:
            arr = np.asarray(val)
            if arr.ndim >= 1 and arr.dtype.kind != 'O':
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
                    [head[j] if j < len(head) else j for j in row]
                    for row in trimmed
                ]
            atoms.info[store_key] = trimmed
        else:
            atoms.info[store_key] = val


def ensure_neighbors(atoms: Atoms):
    """Check that neighbors have been computed for these atoms."""
    if _PREFIX + "neighbors" not in atoms.info and _PREFIX + "neighbors" not in atoms.arrays:
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


def pad_atoms_for_neighbor_finding(atoms: Atoms):
    """
    If the cell has too few atoms (< 200) or is too small (< 10 Å per side),
    create a padded system using ASE's repeat() for the C++ neighbor search.
    
    Ghost atoms are marked with ghost=True so results can be trimmed to
    the original atoms.

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
    cell = np.array(atoms.cell)

    if n >= _MIN_ATOMS:
        # No padding needed
        d = atoms_to_dict(atoms)
        return d, get_box_params(atoms), n

    # Compute repetitions needed  
    needed = max(int(np.ceil((_MIN_ATOMS / n) ** (1.0 / 3.0))), 2)
    reps = [needed, needed, needed]
    
    # Also ensure box is large enough per side
    for i in range(3):
        side = np.linalg.norm(cell[i]) * reps[i]
        if side < _MIN_BOX_SIDE:
            reps[i] = max(reps[i], int(np.ceil(_MIN_BOX_SIDE / np.linalg.norm(cell[i]))))

    # Create a repeated supercell using ASE
    supercell = atoms.repeat(reps)
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
