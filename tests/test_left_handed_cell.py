"""Left-handed cells (negative determinant) must work for the adaptive,
SANN and number methods (regression: signed volume gave a negative search
radius and no candidates)."""
import numpy as np
import pytest
from ase.build import bulk

import pyscal3


def _left_handed_fcc():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    cell = np.array(atoms.cell)
    cell[2] *= -1
    atoms.set_cell(cell, scale_atoms=False)
    atoms.positions[:, 2] *= -1
    assert np.linalg.det(atoms.cell) < 0
    return atoms


@pytest.mark.parametrize("kwargs", [
    dict(method="cutoff", cutoff=0),
    dict(method="cutoff", cutoff="sann"),
    dict(method="number", nmax=12),
])
def test_left_handed_cell(kwargs):
    atoms = _left_handed_fcc()
    pyscal3.find_neighbors(atoms, **kwargs)
    assert np.all(pyscal3.coordination_number(atoms) == 12)
    q6 = pyscal3.steinhardt_parameter(atoms, l=6)[0]
    assert np.allclose(q6, 0.5745242597140696, atol=1e-6)


def test_left_handed_cell_cna():
    atoms = _left_handed_fcc()
    assert pyscal3.common_neighbor_analysis(atoms)["fcc"] == len(atoms)
