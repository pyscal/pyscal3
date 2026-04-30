"""Tests for coordination number variants."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pyscal3
from pyscal3.structures import make_crystal


# ------------------------------------------------------------------
# Coordination Number
# ------------------------------------------------------------------

def test_cn_fcc():
    """FCC should have CN = 12."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    cn = pyscal3.coordination_number(atoms)
    assert np.all(cn == 12)


def test_cn_bcc():
    """BCC should have CN = 14 (8+6) with adaptive cutoff."""
    atoms = make_crystal("bcc", lattice_constant=2.87, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    cn = pyscal3.coordination_number(atoms)
    assert np.all(cn == 14)


def test_cn_stores_in_arrays():
    """CN should be stored in atoms.arrays."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.coordination_number(atoms)
    assert "pyscal_cn" in atoms.arrays


# ------------------------------------------------------------------
# Effective Coordination Number (ECoN)
# ------------------------------------------------------------------

def test_econ_fcc_perfect():
    """ECoN for perfect FCC should be exactly 12 (all neighbors equidistant)."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    econ = pyscal3.effective_coordination_number(atoms)
    assert_allclose(econ, 12.0, atol=1e-6)


def test_econ_bcc():
    """ECoN for BCC should be less than 14 (second shell neighbors weighted less)."""
    atoms = make_crystal("bcc", lattice_constant=2.87, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    econ = pyscal3.effective_coordination_number(atoms)
    # ECoN < CN because second-shell neighbors have reduced weight
    assert econ[0] < 14.0
    assert econ[0] > 8.0  # but well above 8 (nearest-only)


def test_econ_noise_decreases():
    """Adding noise should change ECoN from the perfect value."""
    perfect = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(perfect, method="cutoff", cutoff=0)
    econ_perf = pyscal3.effective_coordination_number(perfect)
    assert econ_perf.std() < 1e-10  # all identical for perfect

    noisy = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3), noise=0.05)
    # Use explicit cutoff for noisy structures (adaptive can be unreliable)
    pyscal3.find_neighbors(noisy, method="cutoff", cutoff=3.5)
    econ_noisy = pyscal3.effective_coordination_number(noisy)
    assert econ_noisy.std() > 0.01


def test_econ_stores_in_arrays():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.effective_coordination_number(atoms)
    assert "pyscal_econ" in atoms.arrays


# ------------------------------------------------------------------
# Generalized Coordination Number (GCN)
# ------------------------------------------------------------------

def test_gcn_fcc_bulk():
    """GCN for bulk FCC = CN = 12 (all neighbors fully coordinated)."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    gcn = pyscal3.generalized_coordination_number(atoms)
    assert_allclose(gcn, 12.0, atol=1e-6)


def test_gcn_custom_cnmax():
    """GCN with explicit cn_max should scale properly."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    gcn = pyscal3.generalized_coordination_number(atoms, cn_max=12)
    assert_allclose(gcn, 12.0, atol=1e-6)


def test_gcn_stores_in_arrays():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.generalized_coordination_number(atoms)
    assert "pyscal_gcn" in atoms.arrays


# ------------------------------------------------------------------
# Local Density
# ------------------------------------------------------------------

def test_local_density_fcc():
    """Local density for FCC should be uniform in a perfect crystal."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    rho = pyscal3.local_density(atoms)
    assert rho.std() / rho.mean() < 1e-10  # uniform


def test_local_density_fcc_vs_bcc():
    """BCC with smaller lattice should have higher local density."""
    fcc = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    bcc = make_crystal("bcc", lattice_constant=2.87, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(fcc, method="cutoff", cutoff=0)
    pyscal3.find_neighbors(bcc, method="cutoff", cutoff=0)

    rho_fcc = pyscal3.local_density(fcc)
    rho_bcc = pyscal3.local_density(bcc)
    # BCC has smaller distances -> higher density estimate
    assert rho_bcc.mean() > rho_fcc.mean()


def test_local_density_stores_in_arrays():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.local_density(atoms)
    assert "pyscal_local_density" in atoms.arrays
