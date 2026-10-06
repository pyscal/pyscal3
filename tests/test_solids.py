"""Inputs and stored results of find_solids and find_clusters."""
import numpy as np

import pyscal3
from pyscal3.structures import make_crystal


def test_find_solids_accepts_numpy_numbers():
    # bonds may come from numpy, for example from a parameter scan
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    expected = pyscal3.find_solids(atoms, bonds=6, cluster=True)
    assert pyscal3.find_solids(atoms, bonds=np.int64(6), cluster=True) == expected
    assert pyscal3.find_solids(atoms, bonds=np.float32(0.5), cluster=True) == len(atoms)


def test_find_clusters_clears_the_largest_cluster_mask():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    assert pyscal3.find_clusters(atoms, condition=np.ones(len(atoms), bool)) == len(atoms)
    assert pyscal3.find_clusters(atoms, condition=np.zeros(len(atoms), bool)) == 0
    assert not atoms.arrays["pyscal_largest_cluster"].any()
