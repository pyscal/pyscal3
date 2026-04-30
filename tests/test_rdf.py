"""Tests for radial distribution function."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_rdf_bcc():
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(10, 10, 10))
    rdf, r = pyscal3.radial_distribution_function(atoms, rmax=2)
    args = np.argsort(rdf)[::-1]
    assert r[args[0]] - 0.86 < 1e-5


def test_rdf_fcc():
    atoms = make_crystal("fcc", lattice_constant=1.0, repetitions=(10, 10, 10))
    rdf, r = pyscal3.radial_distribution_function(atoms, rmax=2)
    args = np.argsort(rdf)[::-1]
    assert r[args[0]] - 0.70 < 1e-5
