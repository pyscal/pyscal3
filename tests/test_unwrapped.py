"""Unwrapped coordinates (atoms several cells away from the primary box, as
in long MD runs) must give the same neighbors as wrapped ones."""
import numpy as np
import pytest
from ase.build import bulk

import pyscal3


def _sets(atoms):
    nb = atoms.arrays.get("pyscal_neighbors")
    if nb is None:
        nb = atoms.info["pyscal_neighbors"]
    return [frozenset(int(j) for j in row) for row in nb]


def _unwrap(atoms, seed=0):
    rng = np.random.default_rng(seed)
    shifts = rng.integers(-3, 4, size=(len(atoms), 3))
    out = atoms.copy()
    out.positions = atoms.positions + shifts @ np.array(atoms.cell)
    return out


@pytest.mark.parametrize("builder", [
    lambda: bulk("Cu", "fcc", cubic=True).repeat(4),
    lambda: bulk("Cu", "fcc").repeat(6),          # 60-degree primitive cell
])
@pytest.mark.parametrize("kwargs", [
    dict(method="cutoff", cutoff=3.0, cells=False),
    dict(method="cutoff", cutoff=3.0, cells=True),
    dict(method="cutoff", cutoff=0),
    dict(method="number", nmax=12),
])
def test_unwrapped_positions_same_neighbors(builder, kwargs):
    ref = builder()
    moved = _unwrap(ref)
    pyscal3.find_neighbors(ref, **kwargs)
    pyscal3.find_neighbors(moved, **kwargs)
    assert _sets(ref) == _sets(moved)
    assert np.allclose(pyscal3.steinhardt_parameter(ref, l=6)[0],
                       pyscal3.steinhardt_parameter(moved, l=6)[0])


def test_unwrapped_positions_voronoi():
    ref = bulk("Mg", "hcp").repeat((4, 4, 3))
    moved = _unwrap(ref, seed=1)
    for a in (ref, moved):
        pyscal3.find_neighbors(a, method="voronoi")
    assert _sets(ref) == _sets(moved)
    assert np.allclose(ref.arrays["pyscal_voronoi_volume"], moved.arrays["pyscal_voronoi_volume"])


def test_get_distance_unwrapped():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    a, b = atoms.positions[0], atoms.positions[1]
    d0 = pyscal3.get_distance(atoms, a, b)
    d1 = pyscal3.get_distance(atoms, a, b + 5 * np.array(atoms.cell)[0] - 2 * np.array(atoms.cell)[2])
    assert np.isclose(d0, d1)
