"""Tests for centrosymmetry parameter."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal
from ase.build import bulk


def test_cs_bcc_nmax8():
    """Perfect BCC with nmax=8 should give zero centrosymmetry."""
    atoms = make_crystal("bcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    cs = pyscal3.centrosymmetry(atoms, nmax=8)
    assert round(np.mean(cs), 2) == 0.00


def test_cs_bcc_nmax4():
    """BCC with nmax=4 (wrong for BCC) should give non-zero centrosymmetry."""
    atoms = make_crystal("bcc", lattice_constant=4.00, repetitions=(7, 7, 7))
    cs = pyscal3.centrosymmetry(atoms, nmax=4)
    assert np.mean(cs) > 0.00


def test_cs_hcp():
    """HCP Ti with nmax=12 should give ~8.7."""
    ti = bulk("Ti", orthorhombic=True)
    cs = pyscal3.centrosymmetry(ti, nmax=12)
    assert round(np.mean(cs), 1) == 8.7
