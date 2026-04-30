"""Tests for Common Neighbor Analysis and Diamond Structure."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal
from ase.build import bulk


def test_cna_fcc_with_lattice_constant():
    atoms = make_crystal("fcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms, lattice_constant=4.00)
    assert res["fcc"] == 7 * 7 * 7 * 4


def test_cna_fcc_adaptive():
    atoms = make_crystal("fcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms)
    assert res["fcc"] == 7 * 7 * 7 * 4


def test_cna_bcc_adaptive():
    atoms = make_crystal("bcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms)
    assert atoms.arrays["pyscal_structure"][0] == 3


def test_cna_hcp_adaptive():
    atoms = make_crystal("hcp", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.common_neighbor_analysis(atoms)
    assert atoms.arrays["pyscal_structure"][0] == 2


def test_cna_ase_bulks():
    al = bulk("Al")
    res = pyscal3.common_neighbor_analysis(al)
    assert res["fcc"] == 1

    fe = bulk("Fe")
    res = pyscal3.common_neighbor_analysis(fe)
    assert res["bcc"] == 1

    ti = bulk("Ti")
    res = pyscal3.common_neighbor_analysis(ti)
    assert res["hcp"] == 2


def test_diamond_structure():
    atoms = make_crystal("diamond", lattice_constant=4.00, repetitions=(7, 7, 7))
    res = pyscal3.diamond_structure(atoms)
    assert atoms.arrays["pyscal_structure"][0] == 1
