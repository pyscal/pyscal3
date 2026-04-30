"""Tests for Steinhardt bond-order parameters."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_q4_q6_bcc():
    """q4 and q6 for BCC should match known values."""
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
    assert round(np.mean(q4), 2) == 0.51
    assert round(np.mean(q6), 2) == 0.63


def test_q4_q6_bcc_averaged():
    """Averaged q4 and q6 for BCC."""
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
    assert round(np.mean(q4), 2) == 0.51
    assert round(np.mean(q6), 2) == 0.63


def test_q3_bcc_averaged():
    """Averaged q3 for BCC should be zero."""
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    [q3] = pyscal3.steinhardt_parameter(atoms, l=[3], averaged=True)
    assert round(np.mean(q3), 2) == 0.00


def test_q12_fcc():
    """q12 for FCC (normal and averaged)."""
    atoms = make_crystal("fcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    [q12] = pyscal3.steinhardt_parameter(atoms, l=[12])
    assert round(np.mean(q12), 2) == 0.60

    [q12_avg] = pyscal3.steinhardt_parameter(atoms, l=[12], averaged=True)
    assert round(np.mean(q12_avg), 2) == 0.60


def test_q4_q6_via_ase_bulk():
    """q4/q6 for W BCC via ase.build.bulk."""
    from ase.build import bulk
    atoms = bulk("W", a=1.0, cubic=True).repeat([4, 4, 4])
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
    assert round(np.mean(q4), 2) == 0.51
    assert round(np.mean(q6), 2) == 0.63
