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
    # For bcc, SANN accepts the first (8) and second (6) shells:
    # R(14) = (8*2.708 + 6*3.127)/12 = 3.37 < 4.422 (third shell)
    assert dists.shape[1] == 14
    assert np.sum(np.abs(dists[0] - 2.708061437633939) < 1e-5) == 8
    assert np.sum(np.abs(dists[0] - 3.127) < 1e-5) == 6
    assert np.all(dists < 4.0)
    cutoff = atoms.arrays["pyscal_cutoff"]
    assert np.allclose(cutoff, (8 * 2.708061437633939 + 6 * 3.127) / 12)


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


def _reference_sann(atoms, rmax=6.0):
    """Independent SANN implementation using ASE's neighbor list."""
    from ase.neighborlist import neighbor_list
    i, j, d = neighbor_list("ijd", atoms, rmax)
    result = []
    for at in range(len(atoms)):
        sel = i == at
        order = np.argsort(d[sel])
        jj, dd = j[sel][order], d[sel][order]
        m = 3
        while m < len(dd) and dd[:m].sum() / (m - 2) >= dd[m]:
            m += 1
        result.append(frozenset(jj[:m].tolist()))
    return result


def test_sann_matches_reference_on_perturbed_fcc():
    atoms = make_element("Cu", repetitions=(5, 5, 5))
    rng = np.random.default_rng(3)
    atoms.positions += rng.normal(scale=0.08, size=atoms.positions.shape)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff="sann")
    nb = atoms.info.get("pyscal_neighbors", atoms.arrays.get("pyscal_neighbors"))
    got = [frozenset(int(x) for x in row) for row in nb]
    assert got == _reference_sann(atoms)


def test_sann_perfect_fcc_q6():
    atoms = make_element("Cu", repetitions=(5, 5, 5))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff="sann")
    assert np.all(pyscal3.coordination_number(atoms) == 12)
    q6 = pyscal3.steinhardt_parameter(atoms, l=6)[0]
    assert np.allclose(q6, 0.5745242597140696, atol=1e-6)
