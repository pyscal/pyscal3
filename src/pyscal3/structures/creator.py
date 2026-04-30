"""
Crystal structure creation — returns ASE Atoms objects.
"""
import copy
import os

import numpy as np
import yaml
from ase import Atoms

_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

def _load_yaml(filename):
    with open(os.path.join(_DATA_DIR, filename)) as f:
        return yaml.safe_load(f)

_structures = _load_yaml("structure_data.yaml")
_elements = _load_yaml("element_data.yaml")


def available_structures():
    """Return list of available crystal structure names."""
    return list(_structures.keys())


def available_elements():
    """Return list of available element symbols."""
    return list(_elements.keys())


def make_crystal(structure, lattice_constant=1.0, repetitions=None,
                 ca_ratio=1.633, noise=0, element=None, primitive=False):
    """
    Create a crystal structure as an ASE Atoms object.

    Parameters
    ----------
    structure : str
        Crystal type: 'sc', 'bcc', 'fcc', 'hcp', 'dhcp', 'diamond', 'a15', 'l12', 'b2'.
    lattice_constant : float
        Lattice constant in Angstroms. Default 1.0.
    repetitions : int or tuple of 3 ints, optional
        Unit cell repetitions. Default (1,1,1).
    ca_ratio : float
        c/a ratio for hcp/dhcp. Default 1.633.
    noise : float
        Standard deviation of Gaussian noise on positions. Default 0.
    element : str or list of str, optional
        Chemical element(s).
    primitive : bool
        If True, use primitive cell. Default False.

    Returns
    -------
    ase.Atoms
    """
    if repetitions is None:
        repetitions = (1, 1, 1)
    elif isinstance(repetitions, int):
        repetitions = (repetitions, repetitions, repetitions)

    if structure not in _structures:
        raise ValueError(f"Unknown structure '{structure}'. Available: {available_structures()}")

    cell_type = 'primitive' if primitive else 'conventional'
    if cell_type not in _structures[structure]:
        alt = 'conventional' if primitive else 'primitive'
        raise ValueError(f"'{cell_type}' cell not available for '{structure}', try {alt}")

    sdict = copy.deepcopy(_structures[structure][cell_type])

    # Handle elements
    unique_types = sorted(set(sdict["species"]))
    if element is not None:
        if isinstance(element, str):
            element = [element]
        if len(element) != len(unique_types):
            raise ValueError(f"Need {len(unique_types)} element(s), got {len(element)}")
        element_map = dict(zip(unique_types, element))
    else:
        element_map = {t: None for t in unique_types}

    # Fix c/a ratio for hcp
    if structure in ('hcp', 'dhcp'):
        sdict['box'][2][2] = ca_ratio

    box = lattice_constant * np.array(sdict["box"], dtype=float)
    positions = np.array([_unfold(p, box) for p in sdict["positions"]])
    if noise > 0:
        positions += np.random.normal(0, noise, positions.shape)

    types = list(sdict['species'])
    species = [element_map[t] for t in types]
    numbers = list(range(len(types)))

    # Replicate
    all_pos = [positions]
    all_types = [types]
    all_species = [species]

    for d in range(3):
        new_pos = []
        new_types = []
        new_species = []
        current_pos = np.concatenate(all_pos)
        current_types = sum(all_types, [])
        current_species = sum(all_species, [])

        for i in range(1, repetitions[d]):
            shifted = current_pos + i * box[d]
            if noise > 0:
                shifted += np.random.normal(0, noise, shifted.shape)
            new_pos.append(shifted)
            new_types.append(current_types)
            new_species.append(current_species)

        all_pos = [np.concatenate(all_pos)] + new_pos
        all_types = [sum(all_types, [])] + new_types
        all_species = [sum(all_species, [])] + new_species

    final_pos = np.concatenate(all_pos)
    final_types = sum(all_types, [])
    final_species = sum(all_species, [])

    # Build box
    cell = np.array([
        repetitions[0] * box[0],
        repetitions[1] * box[1],
        repetitions[2] * box[2],
    ])

    # Build ASE Atoms
    has_elements = any(s is not None for s in final_species)
    if has_elements:
        atoms = Atoms(
            symbols=final_species,
            positions=final_pos,
            cell=cell,
            pbc=True,
        )
    else:
        atoms = Atoms(
            positions=final_pos,
            cell=cell,
            pbc=True,
        )
        # Store integer types as tags
        atoms.set_tags(final_types)

    return atoms


def make_general_lattice(positions, types, box, lattice_constant=1.0,
                         repetitions=None, noise=0, element=None):
    """
    Create a custom lattice from fractional positions.

    Parameters
    ----------
    positions : list of [x,y,z]
        Fractional coordinates.
    types : list of int
        Per-atom type indices.
    box : list of 3 vectors
        Unit cell vectors (before lattice_constant scaling).
    lattice_constant : float
        Scaling factor. Default 1.0.
    repetitions : tuple of 3 ints
        Repetitions. Default (1,1,1).
    noise : float
        Gaussian noise std dev. Default 0.
    element : str or list of str, optional
        Chemical element(s).

    Returns
    -------
    ase.Atoms
    """
    if len(types) != len(positions):
        raise ValueError("types and positions must have same length")

    # Inject as a custom structure
    custom_structures = {
        "custom": {
            "conventional": {
                "species": types,
                "positions": positions,
                "box": box,
            }
        }
    }

    # Temporarily replace the structure database
    global _structures
    old = _structures
    _structures = custom_structures
    try:
        result = make_crystal("custom", lattice_constant=lattice_constant,
                              repetitions=repetitions, noise=noise, element=element)
    finally:
        _structures = old

    return result


def make_element(symbol, repetitions=None, primitive=False):
    """
    Create a crystal of a specific element using its known structure.

    Parameters
    ----------
    symbol : str
        Element symbol, e.g. 'Fe', 'Cu', 'Al'.
    repetitions : tuple of 3 ints, optional
        Default (1,1,1).
    primitive : bool
        Use primitive cell. Default False.

    Returns
    -------
    ase.Atoms
    """
    if symbol not in _elements:
        raise ValueError(f"Unknown element '{symbol}'. Available: {available_elements()}")

    info = _elements[symbol]
    return make_crystal(
        info['structure'],
        lattice_constant=info['lattice_constant'],
        repetitions=repetitions,
        element=symbol,
        primitive=primitive,
    )


def _unfold(p, box):
    """Convert fractional to Cartesian: p[0]*box[0] + p[1]*box[1] + p[2]*box[2]."""
    return p[0] * box[0] + p[1] * box[1] + p[2] * box[2]
