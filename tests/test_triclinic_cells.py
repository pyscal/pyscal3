"""Neighbor finding must give identical results with and without cell lists,
also for triclinic / hexagonal / rotated cells (regression for cell lists that
binned raw Cartesian coordinates)."""
import numpy as np
import pytest
from ase.build import bulk, fcc111

import pyscal3

FCC_Q6 = 0.5745242597140696
HCP_Q6 = 0.4842  # ASE Mg has c/a = 1.624, slightly off the ideal 0.4848


def _neighbor_sets(atoms):
    nb = atoms.arrays.get("pyscal_neighbors")
    if nb is None:
        nb = atoms.info["pyscal_neighbors"]
    return [frozenset(int(j) for j in row) for row in nb]


def _rotated_fcc(n):
    atoms = bulk("Cu", "fcc", cubic=True).repeat(n)
    # arbitrary rotation: 30 deg about z then 20 deg about x
    cz, sz = np.cos(np.pi / 6), np.sin(np.pi / 6)
    cx, sx = np.cos(np.pi / 9), np.sin(np.pi / 9)
    rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])
    rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
    rot = rx @ rz
    atoms.set_cell(np.array(atoms.cell) @ rot.T, scale_atoms=False)
    atoms.positions = atoms.positions @ rot.T
    return atoms


CASES = {
    # name: (builder, cutoff, expected cn, expected q6)
    "hcp_hexagonal": (lambda: bulk("Mg", "hcp").repeat((6, 6, 4)), 3.5, 12, HCP_Q6),
    "fcc_primitive": (lambda: bulk("Cu", "fcc").repeat(6), 3.0, 12, FCC_Q6),
    "fcc_rotated": (lambda: _rotated_fcc(5), 3.0, 12, FCC_Q6),
    "fcc_primitive_padded": (lambda: bulk("Cu", "fcc").repeat(4), 3.0, 12, FCC_Q6),
}


@pytest.mark.parametrize("name", list(CASES))
def test_cutoff_cells_match_brute_force(name):
    builder, cutoff, cn, q6 = CASES[name]
    a = builder()
    b = builder()
    pyscal3.find_neighbors(a, method="cutoff", cutoff=cutoff, cells=True)
    pyscal3.find_neighbors(b, method="cutoff", cutoff=cutoff, cells=False)
    assert _neighbor_sets(a) == _neighbor_sets(b)
    assert np.all(pyscal3.coordination_number(a) == cn)
    _check_q6(a, b, q6)


def _check_q6(a, b, q6):
    qa = pyscal3.steinhardt_parameter(a, l=6)[0]
    qb = pyscal3.steinhardt_parameter(b, l=6)[0]
    assert np.allclose(qa, qb)
    assert np.ptp(qa) < 1e-8          # perfect crystal: identical environments
    assert np.allclose(qa, q6, atol=1e-3)


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("method_kwargs", [
    dict(method="cutoff", cutoff=0),          # adaptive
    dict(method="number", nmax=12),
])
def test_adaptive_and_number_cells_triclinic(name, method_kwargs):
    builder, cutoff, cn, q6 = CASES[name]
    a = builder()
    b = builder()
    pyscal3.find_neighbors(a, cells=True, **method_kwargs)
    pyscal3.find_neighbors(b, cells=False, **method_kwargs)
    assert _neighbor_sets(a) == _neighbor_sets(b)
    assert np.all(pyscal3.coordination_number(a) == cn)
    _check_q6(a, b, q6)


def test_shell_cells_triclinic():
    a = bulk("Cu", "fcc").repeat(6)
    b = bulk("Cu", "fcc").repeat(6)
    pyscal3.find_neighbors(a, method="cutoff", cutoff=3.0, shell_thickness=1.0, cells=True)
    pyscal3.find_neighbors(b, method="cutoff", cutoff=3.0, shell_thickness=1.0, cells=False)
    assert _neighbor_sets(a) == _neighbor_sets(b)
    # second shell of fcc has 6 atoms
    assert np.all(pyscal3.coordination_number(a) == 6)


def test_cna_on_padded_primitive_cell():
    """Used to hang: cell lists lost all candidates and the failure was ignored."""
    a = bulk("Cu", "fcc").repeat(4)
    res = pyscal3.common_neighbor_analysis(a)
    assert res["fcc"] == len(a)
    cs = pyscal3.centrosymmetry(a)
    assert np.all(cs < 1e-8)


def test_cna_hexagonal_cell():
    a = bulk("Mg", "hcp").repeat((6, 6, 4))
    res = pyscal3.common_neighbor_analysis(a)
    assert res["hcp"] == len(a)


def test_fcc111_slab_hexagonal_cell():
    slab = fcc111("Cu", size=(4, 4, 4), vacuum=10.0)
    pyscal3.find_neighbors(slab, method="cutoff", cutoff=3.0)
    cn = pyscal3.coordination_number(slab)
    values, counts = np.unique(cn, return_counts=True)
    assert values.tolist() == [9, 12]
    assert counts.tolist() == [32, 32]


def test_cells_with_cutoff_larger_than_box_side_does_not_crash():
    atoms = bulk("Cu", "fcc", cubic=True).repeat((1, 1, 30))
    a, b = atoms.copy(), atoms.copy()
    pyscal3.find_neighbors(a, method="cutoff", cutoff=7.0, cells=True)
    pyscal3.find_neighbors(b, method="cutoff", cutoff=7.0, cells=False)
    assert _neighbor_sets(a) == _neighbor_sets(b)
