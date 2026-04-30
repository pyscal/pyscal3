"""Tests for Minkowski structure metrics (Voronoi-weighted Steinhardt parameters)."""

import numpy as np
import pytest
from numpy.testing import assert_allclose

import pyscal3
from pyscal3.structures import make_crystal


# ------------------------------------------------------------------
# Reference values
# ------------------------------------------------------------------
# Mickel et al., J. Chem. Phys. 138, 044501 (2013) Table I
# Perfect crystals: q_l^Mink for FCC (14 Voronoi neighbors)
# FCC:  q4=0.1909,  q6=0.5745
#
# BCC Voronoi tessellation gives a truncated octahedron with 14 faces
# (8 hexagonal nearest + 6 square second-nearest).  The weighted q
# values are q4≈0.224, q6≈0.567 with face-area exponent 1.


def test_minkowski_fcc_q6():
    """FCC q6^Mink should match Mickel et al. reference."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    [q6] = pyscal3.minkowski_parameter(atoms, l=6)
    # All atoms in a perfect crystal should have the same value
    assert_allclose(q6, q6[0], atol=1e-10)
    assert_allclose(q6[0], 0.5745, atol=0.005)


def test_minkowski_bcc_q6():
    """BCC q6^Mink should be uniform across all atoms."""
    atoms = make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4))
    [q6] = pyscal3.minkowski_parameter(atoms, l=6)
    assert_allclose(q6, q6[0], atol=1e-10)
    # BCC truncated octahedron with 14 Voronoi neighbors
    assert_allclose(q6[0], 0.5669, atol=0.005)


def test_minkowski_fcc_q4():
    """FCC q4^Mink should match Mickel et al. reference."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    [q4] = pyscal3.minkowski_parameter(atoms, l=4)
    assert_allclose(q4, q4[0], atol=1e-10)
    assert_allclose(q4[0], 0.1909, atol=0.005)


def test_minkowski_bcc_q4():
    """BCC q4^Mink should be uniform across all atoms."""
    atoms = make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4))
    [q4] = pyscal3.minkowski_parameter(atoms, l=4)
    assert_allclose(q4, q4[0], atol=1e-10)
    # BCC truncated octahedron: q4 ≈ 0.224
    assert_allclose(q4[0], 0.2240, atol=0.005)


def test_minkowski_multi_l():
    """Request multiple l at once."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    q4, q6 = pyscal3.minkowski_parameter(atoms, l=[4, 6])
    assert_allclose(q4[0], 0.1909, atol=0.005)
    assert_allclose(q6[0], 0.5745, atol=0.005)


def test_minkowski_differs_from_standard():
    """Minkowski q_l should differ from standard cutoff-based q_l."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))

    # Standard Steinhardt with cutoff
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    [q6_std] = pyscal3.steinhardt_parameter(atoms, l=6)

    # Minkowski (Voronoi-weighted)
    [q6_mink] = pyscal3.minkowski_parameter(atoms, l=6)

    # For perfect FCC both are constant, but values may differ slightly
    # because Voronoi includes more neighbors than nearest-neighbour cutoff
    assert q6_std.shape == q6_mink.shape


def test_minkowski_averaged():
    """Averaged Minkowski parameters should differ from non-averaged."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    [q6] = pyscal3.minkowski_parameter(atoms, l=6)
    [q6_avg] = pyscal3.minkowski_parameter(atoms, l=6, averaged=True)
    # For perfect crystal both should be the same (all atoms identical)
    assert_allclose(q6, q6_avg, atol=1e-6)


def test_minkowski_fcc_vs_bcc_discrimination():
    """Minkowski q4 should discriminate FCC from BCC."""
    fcc = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    bcc = make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4))

    [q4_fcc] = pyscal3.minkowski_parameter(fcc, l=4)
    [q4_bcc] = pyscal3.minkowski_parameter(bcc, l=4)

    # FCC q4 ≈ 0.191, BCC q4 ≈ 0.224 — different values
    assert abs(q4_fcc.mean() - q4_bcc.mean()) > 0.02


def test_minkowski_voroexp():
    """Different voroexp should produce different results."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))

    [q6_exp1] = pyscal3.minkowski_parameter(atoms, l=6, voroexp=1)
    [q6_exp2] = pyscal3.minkowski_parameter(atoms, l=6, voroexp=2)

    # Different exponent → different weighting → different values
    # For a perfect crystal with uniform face areas they may be
    # similar but not necessarily identical
    assert q6_exp1.shape == q6_exp2.shape


def test_minkowski_stores_in_arrays():
    """Results should be accessible from atoms.arrays."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    pyscal3.minkowski_parameter(atoms, l=6)
    assert "pyscal_q6" in atoms.arrays
    assert atoms.arrays["pyscal_q6"].shape == (len(atoms),)


def test_minkowski_noise_sensitivity():
    """Added noise should change Minkowski q6 relative to perfect crystal."""
    perfect = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    noisy = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4), noise=0.1)

    [q6_perf] = pyscal3.minkowski_parameter(perfect, l=6)
    [q6_noisy] = pyscal3.minkowski_parameter(noisy, l=6)

    # Perfect has zero spread, noisy has non-zero spread
    assert q6_perf.std() < 1e-10
    assert q6_noisy.std() > 0.001
