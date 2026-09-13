"""Tests for the Warren-Cowley short-range order parameter."""
import numpy as np
import pytest
import pyscal3
from pyscal3.structures import make_crystal


def test_sro_l12_ordered():
    # Cu3Au (type 1 = Au corner, type 2 = Cu faces):
    # each Cu has 4 Au + 8 Cu neighbours, c_Au = 1/4
    atoms = make_crystal("l12", lattice_constant=3.75, repetitions=(3, 3, 3),
                         element=["Au", "Cu"])
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    per_atom = pyscal3.short_range_order(atoms, "Cu", "Au", average=False)
    is_cu = np.array(atoms.get_chemical_symbols()) == "Cu"
    assert np.allclose(per_atom[is_cu], 1 - (1 / 3) / (1 / 4))
    assert np.all(np.isnan(per_atom[~is_cu]))
    assert np.isclose(pyscal3.short_range_order(atoms, "Cu", "Au"), -1 / 3)
    # Au atoms: all 12 neighbours are Cu -> p = 1, c_Cu = 3/4
    assert np.isclose(pyscal3.short_range_order(atoms, "Au", "Cu"), 1 - 1 / (3 / 4))
    # like pairs around Au: none
    assert np.isclose(pyscal3.short_range_order(atoms, "Au", "Au"), 1.0)


def test_sro_defaults_and_numbers():
    """Defaults: reference = majority species (type 2), compare = type 1."""
    atoms = make_crystal("l12", lattice_constant=3.75, repetitions=(2, 2, 2))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    assert np.isclose(pyscal3.short_range_order(atoms), -1 / 3)
    per_atom = pyscal3.short_range_order(atoms, average=False)
    assert np.all(np.isnan(per_atom[atoms.get_atomic_numbers() == 1]))
    # minority corner atoms: all 12 neighbours are type 2, c_2 = 3/4
    assert np.isclose(pyscal3.short_range_order(atoms, reference_type=1, compare_type=2), 1 - 1 / 0.75)
    # like pairs around the corner atoms: none
    assert np.isclose(pyscal3.short_range_order(atoms, reference_type=1, compare_type=1), 1.0)


def test_sro_b2():
    atoms = make_crystal("b2", lattice_constant=2.87, repetitions=(3, 3, 3),
                         element=["Ni", "Al"])
    pyscal3.find_neighbors(atoms, method="number", nmax=8)
    assert np.isclose(pyscal3.short_range_order(atoms, "Ni", "Al"), -1.0)


def test_sro_random_alloy_near_zero():
    atoms = make_crystal("l12", lattice_constant=3.75, repetitions=(6, 6, 6),
                         element=["Au", "Cu"])
    rng = np.random.default_rng(42)
    symbols = np.array(atoms.get_chemical_symbols())
    rng.shuffle(symbols)
    atoms.set_chemical_symbols(symbols.tolist())
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    assert abs(pyscal3.short_range_order(atoms, "Cu", "Au")) < 0.05


def test_sro_errors():
    atoms = make_crystal("fcc", lattice_constant=3.6, repetitions=(2, 2, 2), element="Cu")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    with pytest.raises(ValueError, match="two species"):
        pyscal3.short_range_order(atoms)
    with pytest.raises(ValueError, match="No atoms"):
        pyscal3.short_range_order(atoms, "Cu", "Au")
