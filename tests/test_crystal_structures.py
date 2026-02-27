"""Tests for crystal structure creation."""
import pyscal3
from pyscal3.structures import make_crystal, make_element, available_structures


def test_available_structures():
    structs = available_structures()
    assert "bcc" in structs
    assert "fcc" in structs
    assert "hcp" in structs


def test_create_crystal_bcc():
    atoms = make_crystal("bcc")
    assert len(atoms) == 2


def test_create_crystal_fcc():
    atoms = make_crystal("fcc")
    assert len(atoms) == 4


def test_create_crystal_hcp():
    atoms = make_crystal("hcp")
    assert len(atoms) == 4


def test_create_crystal_diamond():
    atoms = make_crystal("diamond")
    assert len(atoms) == 8


def test_create_crystal_a15():
    atoms = make_crystal("a15")
    assert len(atoms) == 8


def test_create_crystal_l12():
    atoms = make_crystal("l12")
    assert len(atoms) == 4


def test_create_crystal_b2():
    atoms = make_crystal("b2")
    assert len(atoms) == 2


def test_create_crystal_with_repetitions():
    atoms = make_crystal("bcc", repetitions=(3, 3, 3))
    assert len(atoms) == 2 * 27


def test_make_element_cu():
    atoms = make_element("Cu")
    assert len(atoms) > 0


def test_make_element_fe():
    atoms = make_element("Fe")
    assert len(atoms) > 0


def test_make_element_si():
    atoms = make_element("Si")
    assert len(atoms) > 0
