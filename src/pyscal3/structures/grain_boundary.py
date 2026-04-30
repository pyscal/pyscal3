"""
Grain boundary creation — returns ASE Atoms objects.
"""
import numpy as np
import os
import yaml

from ase import Atoms
from pyscal3 import csl as _csl

_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

def _load_yaml(filename):
    with open(os.path.join(_DATA_DIR, filename)) as f:
        return yaml.safe_load(f)

_structures = _load_yaml("structure_data.yaml")
_elements = _load_yaml("element_data.yaml")


def make_grain_boundary(axis, sigma, gb_plane, structure=None, element=None,
                        lattice_constant=1, repetitions=(1, 1, 1), overlap=0.0):
    """
    Create a symmetric tilt grain boundary as an ASE Atoms object.

    Parameters
    ----------
    axis : list of 3 int
        Tilt axis, e.g. [1, 0, 0].
    sigma : int
        CSL sigma value.
    gb_plane : list of 3 int
        Grain boundary plane Miller indices.
    structure : str, optional
        Crystal structure ('bcc', 'fcc', etc.). Required if element not given.
    element : str, optional
        Element symbol. Structure and lattice constant auto-determined.
    lattice_constant : float
        Lattice constant. Default 1.
    repetitions : tuple of 3 int
        Supercell repetitions. Default (1,1,1).
    overlap : float
        Atom overlap distance to remove. Default 0.0.

    Returns
    -------
    ase.Atoms
    """
    axis = np.array(axis)

    # Determine structure parameters
    if element is not None and element in _elements:
        structure = _elements[element]['structure']
        lattice_constant = _elements[element]['lattice_constant']
    elif structure is None:
        raise ValueError("Either 'structure' or 'element' must be provided")

    if structure not in _structures:
        raise ValueError(f"Unknown structure '{structure}'")
    if 'conventional' not in _structures[structure]:
        raise ValueError("Grain boundaries require a conventional cell")

    basis = _structures[structure]['conventional']['positions']

    # Create the grain boundary geometry
    theta, m, n = _csl.get_theta_m_n_list(axis, sigma)[0]
    R = _csl.rot(axis, theta)
    M1, M2 = _csl.create_minimal_cell_method_1(sigma, axis, R)

    validity, ortho_1, ortho_2 = _csl.find_orthogonal_cell(
        axis, sigma, theta, R, m, n, gb_plane, M1, M2, tol=0.001)

    if not validity:
        raise ValueError("Cannot create GB with the given parameters")

    # Populate with atoms
    box, atoms1, atoms2 = _csl.populate_gb(
        ortho_1, ortho_2, np.array(basis), lattice_constant,
        dim=repetitions, overlap=overlap)

    total_atoms = np.concatenate((atoms1, atoms2))

    # Build ASE Atoms
    cell = np.array(box)
    if element is not None:
        symbols = [element] * len(total_atoms)
        ase_atoms = Atoms(symbols=symbols, positions=total_atoms, cell=cell, pbc=True)
    else:
        ase_atoms = Atoms(positions=total_atoms, cell=cell, pbc=True)

    return ase_atoms
