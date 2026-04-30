"""
pyscal3.structures — Crystal structure creation utilities.

Create crystal structures, custom lattices, and grain boundaries 
as ASE Atoms objects.

Example
-------
>>> from pyscal3.structures import make_crystal, make_element
>>> atoms = make_crystal("bcc", lattice_constant=2.87, repetitions=(4,4,4))
>>> atoms = make_element("Fe", repetitions=(4,4,4))
"""
from pyscal3.structures.creator import (
    make_crystal,
    make_general_lattice, 
    make_element,
    available_structures,
    available_elements,
)
from pyscal3.structures.grain_boundary import make_grain_boundary

__all__ = [
    "make_crystal",
    "make_general_lattice",
    "make_element",
    "make_grain_boundary",
    "available_structures",
    "available_elements",
]
