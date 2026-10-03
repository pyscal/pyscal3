"""Tests for Ackland-Jones structure classification."""

from pathlib import Path

import numpy as np
import pytest

import pyscal3
from pyscal3.structures import make_crystal


def _classify(name, lattice_constant, noise=0.0, repetitions=(4, 4, 4)):
    """Helper to create a crystal, find neighbors, and classify.

    make_crystal draws the noise from numpy's global generator, so it is
    seeded here: unseeded, about 1 in 1000 noisy fcc samples falls below the
    90% threshold of test_ackland_fcc_noisy.
    """
    np.random.seed(12345)
    atoms = make_crystal(name, lattice_constant=lattice_constant,
                         repetitions=repetitions, noise=noise)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    labels, names = pyscal3.identify_ackland_jones(atoms)
    return atoms, labels, names


# ------------------------------------------------------------------
# Perfect crystals
# ------------------------------------------------------------------

def test_ackland_fcc():
    """Perfect FCC should be classified as FCC (label 1)."""
    atoms, labels, names = _classify("fcc", 4.05)
    assert np.all(labels == 1)
    assert all(n == "fcc" for n in names)


def test_ackland_bcc():
    """Perfect BCC should be classified as BCC (label 3)."""
    atoms, labels, names = _classify("bcc", 2.87)
    assert np.all(labels == 3)
    assert all(n == "bcc" for n in names)


def test_ackland_hcp():
    """Perfect HCP should be classified as HCP (label 2)."""
    atoms, labels, names = _classify("hcp", 3.21)
    assert np.all(labels == 2)
    assert all(n == "hcp" for n in names)


# ------------------------------------------------------------------
# With thermal noise
# ------------------------------------------------------------------

def test_ackland_fcc_noisy():
    """Most atoms in a noisy FCC should still be classified as FCC."""
    _, labels, _ = _classify("fcc", 4.05, noise=0.05)
    fcc_frac = (labels == 1).sum() / len(labels)
    assert fcc_frac > 0.9, f"Only {fcc_frac:.0%} classified as FCC"


def test_ackland_bcc_noisy():
    """Noisy BCC should have at least some atoms classified as BCC.

    Note: BCC is more sensitive to noise than FCC/HCP because the
    antiparallel-pair count (χ₀) easily drops from 7 to 6 when atoms
    are displaced, causing misclassification.  This is a known
    limitation of chi-parameter based methods.
    """
    _, labels, _ = _classify("bcc", 2.87, noise=0.05)
    bcc_frac = (labels == 3).sum() / len(labels)
    assert bcc_frac > 0.1, f"Only {bcc_frac:.0%} classified as BCC"


def test_ackland_hcp_noisy():
    """Most atoms in a noisy HCP should still be classified as HCP."""
    _, labels, _ = _classify("hcp", 3.21, noise=0.05)
    hcp_frac = (labels == 2).sum() / len(labels)
    assert hcp_frac > 0.9, f"Only {hcp_frac:.0%} classified as HCP"


# ------------------------------------------------------------------
# Return format and storage
# ------------------------------------------------------------------

def test_ackland_returns_labels_and_names():
    """Should return (labels_array, names_list)."""
    atoms, labels, names = _classify("fcc", 4.05)
    assert isinstance(labels, np.ndarray)
    assert labels.dtype in (np.int32, np.int64, int)
    assert isinstance(names, list)
    assert len(labels) == len(atoms)
    assert len(names) == len(atoms)


def test_ackland_stores_in_arrays():
    """Results should be stored in atoms.arrays."""
    atoms, labels, names = _classify("fcc", 4.05)
    assert "pyscal_ackland_label" in atoms.arrays
    assert "pyscal_structure" in atoms.arrays
    np.testing.assert_array_equal(atoms.arrays["pyscal_ackland_label"], labels)


# ------------------------------------------------------------------
# Discrimination
# ------------------------------------------------------------------

def test_ackland_fcc_bcc_discrimination():
    """FCC and BCC should get different labels."""
    _, fcc_labels, _ = _classify("fcc", 4.05)
    _, bcc_labels, _ = _classify("bcc", 2.87)
    assert fcc_labels[0] != bcc_labels[0]


def test_ackland_fcc_hcp_discrimination():
    """FCC and HCP should get different labels."""
    _, fcc_labels, _ = _classify("fcc", 4.05)
    _, hcp_labels, _ = _classify("hcp", 3.21)
    assert fcc_labels[0] != hcp_labels[0]


# ------------------------------------------------------------------
# The method of Ackland and Jones (2006)
# ------------------------------------------------------------------

def _reference(atoms, rcut=6.0):
    """Independent implementation on ASE neighbor lists, following the paper
    and the original LAMMPS compute ackland/atom."""
    from ase.neighborlist import neighbor_list
    i, D = neighbor_list("iD", atoms, rcut)
    edges = [-0.945, -0.915, -0.755, -0.195, 0.195, 0.245, 0.795]
    out = np.zeros(len(atoms), int)
    for a in range(len(atoms)):
        v = D[i == a]
        d2 = (v * v).sum(1)
        order = np.argsort(d2, kind="stable")
        v, d2 = v[order], d2[order]
        r0 = d2[:6].mean()
        n0, n1 = (d2 < 1.45 * r0).sum(), (d2 < 1.55 * r0).sum()
        u = v[:n0] / np.sqrt(d2[:n0])[:, None]
        chi = np.zeros(8, int)
        for x in range(n0):
            for y in range(x + 1, n0):
                chi[np.searchsorted(edges, u[x] @ u[y], side="right")] += 1
        with np.errstate(all="ignore"):
            d_bcc = 0.35 * chi[4] / (chi[5] + chi[6] - chi[4])
        d_cp = abs(1 - chi[6] / 24)
        d_fcc = 0.61 * (abs(chi[0] + chi[1] - 6) + chi[2]) / 6
        d_hcp = (abs(chi[0] - 3) + abs(chi[:4].sum() - 9)) / 12
        if chi[0] == 7:
            d_bcc = 0
        elif chi[0] == 6:
            d_fcc = 0
        elif chi[0] <= 3:
            d_hcp = 0
        if chi[7] > 0:
            out[a] = 0
        elif chi[4] < 3:
            out[a] = 0 if (n1 > 13 or n1 < 11) else 4
        elif d_bcc <= d_cp:
            out[a] = 0 if n1 < 11 else 3
        elif n1 > 12 or n1 < 11:
            out[a] = 0
        else:
            out[a] = 1 if d_fcc < d_hcp else 2
    return out


@pytest.mark.parametrize("name,a", [("fcc", 3.61), ("bcc", 2.87)])
def test_ackland_matches_reference(name, a):
    from ase.build import bulk
    atoms = bulk("Cu", name, a=a, cubic=True).repeat(4)
    atoms.rattle(0.06 * a / np.sqrt(2), seed=3)
    labels, _ = pyscal3.identify_ackland_jones(atoms)
    assert np.array_equal(labels, _reference(atoms))
    # the noise is large enough that several classes occur
    assert len(np.unique(labels)) >= 3


def test_ackland_liquid_is_not_a_crystal():
    from ase.io import read
    path = Path(__file__).resolve().parent.parent / "examples" / "conf.lqd.Al.dump"
    liquid = read(path, format="lammps-dump-text")
    labels, _ = pyscal3.identify_ackland_jones(liquid)
    assert np.mean(labels == 0) > 0.7
    assert np.mean(labels == 2) < 0.2


def test_ackland_stores_integer_labels_and_keeps_neighbors():
    atoms = make_crystal("bcc", lattice_constant=2.87, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.6)
    before = atoms.info["pyscal_bond_neighbors"].copy()
    labels, _ = pyscal3.identify_ackland_jones(atoms)
    assert atoms.arrays["pyscal_structure"].dtype.kind == "i"
    assert np.array_equal(atoms.arrays["pyscal_structure"], labels)
    assert atoms.arrays["pyscal_ackland_chi"].shape == (len(atoms), 8)
    assert np.array_equal(atoms.info["pyscal_bond_neighbors"], before)
    assert np.all(labels == 3)


def test_ackland_needs_no_neighbor_list():
    atoms = make_crystal("hcp", lattice_constant=3.21, repetitions=(3, 3, 3))
    labels, _ = pyscal3.identify_ackland_jones(atoms)
    assert np.all(labels == 2)
