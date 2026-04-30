"""Tests for angular criteria."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_angular_diamond():
    atoms = make_crystal("diamond", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    ang = pyscal3.angular_criteria(atoms)
    assert round(np.mean(ang), 2) == 0.00
