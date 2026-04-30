"""Tests for short-range order parameter."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_sro_l12():
    # L12 has 2 species (e.g., Ni3Al)
    atoms = make_crystal("l12", lattice_constant=4.0, repetitions=(2, 2, 2),
                         element=["Ni", "Al"])
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)

    # Get types from atoms
    types = atoms.get_atomic_numbers()
    unique_types = sorted(np.unique(types))
    assert len(unique_types) == 2

    # SRO for same-type pair (Ni-Ni or Al-Al)
    sro = pyscal3.short_range_order(atoms,
                                     reference_type=int(unique_types[0]),
                                     compare_type=int(unique_types[0]),
                                     average=False)

    # For L12, the SRO should have characteristic non-zero values
    assert np.any(np.abs(sro) > 0.01)
