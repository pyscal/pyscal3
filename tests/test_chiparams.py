"""Tests for chi parameters."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_chiparams_bcc():
    atoms = make_crystal("bcc", lattice_constant=4.0, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    cp = pyscal3.chi_params(atoms)
    # BCC 14 neighbors: C(14,2)=91 pairs, 7 antiparallel (cos=-1) pairs
    assert np.sum(cp[2] - [7, 0, 0, 0, 36, 12, 0, 36, 0]) == 0


def test_chiparams_fcc():
    atoms = make_crystal("fcc", lattice_constant=4.0, repetitions=(5, 5, 5))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    cp = pyscal3.chi_params(atoms)
    assert np.sum(cp[2] - [6, 0, 0, 0, 24, 12, 0, 24, 0]) == 0


def test_chiparams_diamond():
    atoms = make_crystal("diamond", lattice_constant=4.0, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    cp = pyscal3.chi_params(atoms)
    assert np.sum(cp[2] - [0, 0, 0, 0, 6, 0, 0, 0, 0]) == 0


def test_chiparams_ase():
    from ase.build import bulk
    atoms = bulk("W", a=4, cubic=True).repeat([3, 3, 3])
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    cp = pyscal3.chi_params(atoms)
    # BCC 14 neighbors: C(14,2)=91 pairs, 7 antiparallel (cos=-1) pairs
    assert np.sum(cp[0] - [7, 0, 0, 0, 36, 12, 0, 36, 0]) == 0
