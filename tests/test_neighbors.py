"""Tests for neighbor finding."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal, make_element


def test_cutoff_neighbors_bcc():
    atoms = make_crystal("bcc", lattice_constant=3.127, repetitions=(10, 10, 10))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.6)

    dists = atoms.arrays["pyscal_neighbordist"]
    expected = sorted([2.708061437633939, 3.127, 3.127, 3.127, 3.127,
                3.127, 3.127, 2.708061437633939, 2.708061437633939,
                2.708061437633939, 2.708061437633939, 2.708061437633939,
                2.708061437633939, 2.708061437633939])
    assert abs(sum(np.array(sorted(dists[0])) - np.array(expected))) < 1e-5


def test_sann_neighbors_bcc():
    atoms = make_crystal("bcc", lattice_constant=3.127, repetitions=(10, 10, 10))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff="sann")
    dists = atoms.arrays["pyscal_neighbordist"]
    # SANN should find the distance 4.422245809540667 in its neighbor list
    assert any(abs(d - 4.422245809540667) < 1e-5 for d in dists[0])


def test_number_neighbors_bcc():
    atoms = make_crystal("bcc", lattice_constant=3.127, repetitions=(10, 10, 10))
    pyscal3.find_neighbors(atoms, method="number", nmax=8)
    dists = atoms.arrays["pyscal_neighbordist"]
    assert dists.shape[1] == 8
    for d in dists[0]:
        assert abs(d - 2.708061437633939) < 1e-5


def test_neighbor_shell():
    atoms = make_element("Cu", repetitions=(5, 5, 5))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0, shell_thickness=1.0, cells=False)
    dists = atoms.arrays["pyscal_neighbordist"]
    assert dists.shape[1] == 6
    assert abs(dists[0][0] - 3.61) <= 1e-2


def test_distance():
    atoms = make_crystal("bcc", lattice_constant=3.127, repetitions=(2, 2, 2))
    dist = pyscal3.get_distance(atoms, [0.0, 0.0, 0.0], [1.5635, 1.5635, 1.5635])
    assert abs(dist - 2.708061437633939) < 1e-5
