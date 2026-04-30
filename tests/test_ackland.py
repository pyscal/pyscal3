"""Tests for Ackland-Jones structure classification."""

import numpy as np
import pytest

import pyscal3
from pyscal3.structures import make_crystal


def _classify(name, lattice_constant, noise=0.0, repetitions=(4, 4, 4)):
    """Helper to create a crystal, find neighbors, and classify."""
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
