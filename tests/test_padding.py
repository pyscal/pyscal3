"""Ghost-atom padding must be driven by the search radius: neighbors beyond
half the (perpendicular) box width were silently lost before."""
import numpy as np
import pytest
from ase.build import bulk
from ase.neighborlist import neighbor_list

import pyscal3
from pyscal3._bridge import perpendicular_widths, pad_atoms_for_neighbor_finding


def _ase_cn(atoms, cutoff):
    i = neighbor_list("i", atoms, cutoff)
    return np.bincount(i, minlength=len(atoms))


@pytest.mark.parametrize("cutoff", [5.0, 7.5, 9.0])
def test_cutoff_beyond_half_box(cutoff):
    atoms = bulk("Cu", "fcc", cubic=True).repeat(4)  # 256 atoms, L = 14.44
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=cutoff)
    assert np.array_equal(pyscal3.coordination_number(atoms), _ase_cn(atoms, cutoff))


def test_skewed_primitive_cell_third_shell():
    atoms = bulk("Cu", "fcc").repeat(4)  # 64 atoms, 60-degree rhombohedral cell
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.5)
    cn = pyscal3.coordination_number(atoms)
    assert np.all(cn == 42)
    assert np.array_equal(cn, _ase_cn(atoms, 4.5))


def test_shell_beyond_half_box():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=7.5, shell_thickness=1.5)
    expected = _ase_cn(atoms, 9.0) - _ase_cn(atoms, 7.5)
    assert np.array_equal(pyscal3.coordination_number(atoms), expected)


def test_elongated_box_large_cutoff():
    atoms = bulk("Cu", "fcc", cubic=True).repeat((1, 1, 30))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=7.0)
    assert np.array_equal(pyscal3.coordination_number(atoms), _ase_cn(atoms, 7.0))


def test_perpendicular_widths():
    cell = bulk("Mg", "hcp").cell
    w = perpendicular_widths(cell)
    a = np.linalg.norm(cell[0])
    assert np.allclose(w[:2], a * np.sqrt(3) / 2)
    assert np.isclose(w[2], np.linalg.norm(cell[2]))


def test_padding_reps_scale_with_cutoff():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
    d, _, nreal = pad_atoms_for_neighbor_finding(atoms, cutoff=3.0)
    assert len(d["positions"]) == 256 and nreal == 256
    d, _, nreal = pad_atoms_for_neighbor_finding(atoms, cutoff=9.0)
    assert len(d["positions"]) == 256 * 8 and nreal == 256


def test_singular_cell_raises():
    from ase.build import molecule
    with pytest.raises(ValueError, match="periodic cell"):
        pyscal3.find_neighbors(molecule("H2O"), method="cutoff", cutoff=1.2)
