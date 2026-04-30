"""Tests for angular_distribution_function and bond_length_distribution."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3 as pc


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _setup(structure, cubic=True, cutoff=0, method="cutoff"):
    atoms = bulk(structure, cubic=cubic)
    atoms = atoms.repeat(4)
    pc.find_neighbors(atoms, method=method, cutoff=cutoff)
    return atoms


# ---------------------------------------------------------------------------
# ADF tests
# ---------------------------------------------------------------------------


class TestAngularDistribution:
    """Angular Distribution Function tests."""

    def test_adf_fcc_peaks(self):
        """FCC should show peaks near 60, 90, 120 and 180 degrees."""
        atoms = _setup("Al", cutoff=0)
        adf, angles = pc.angular_distribution_function(atoms, bins=180)
        assert adf.shape == (180,)
        assert angles.shape == (180,)
        # Find peak positions (bin centres with highest density)
        peak_bins = np.where(adf > adf.max() * 0.3)[0]
        peak_angles = angles[peak_bins] + 0.5  # bin centres
        # FCC expected: 60, 90, 120, 180
        for expected in [60, 90, 120, 180]:
            assert np.any(np.abs(peak_angles - expected) < 5), (
                f"FCC ADF missing peak near {expected}°; peaks at {peak_angles}"
            )

    def test_adf_bcc_peaks(self):
        """BCC with 8 nearest neighbors shows peaks near 70.5, 109.5, 180."""
        atoms = bulk("Fe", cubic=True).repeat(4)
        # Use number method to get exactly 8 nearest neighbors
        pc.find_neighbors(atoms, method="number", nmax=8)
        adf, angles = pc.angular_distribution_function(atoms, bins=180)
        peak_bins = np.where(adf > adf.max() * 0.3)[0]
        peak_angles = angles[peak_bins] + 0.5
        # BCC 8-nn expected: ~70.5° and ~109.5°
        assert np.any(np.abs(peak_angles - 70.5) < 5), (
            f"BCC ADF missing peak near 70.5°; peaks at {peak_angles}"
        )

    def test_adf_normalised(self):
        """ADF should integrate approximately to 1 (probability density)."""
        atoms = _setup("Al", cutoff=0)
        adf, angles = pc.angular_distribution_function(atoms, bins=180)
        dtheta = 180.0 / 180  # bin width in degrees
        integral = np.sum(adf) * dtheta
        assert abs(integral - 1.0) < 0.05, f"ADF integral = {integral}"

    def test_adf_stored_in_info(self):
        """Results should be stored in atoms.info."""
        atoms = _setup("Al", cutoff=0)
        pc.angular_distribution_function(atoms, bins=90)
        assert "pyscal_adf" in atoms.info
        assert "pyscal_adf_angles" in atoms.info
        assert len(atoms.info["pyscal_adf"]) == 90

    def test_adf_custom_bins(self):
        """Custom bin count should work."""
        atoms = _setup("Al", cutoff=0)
        adf, angles = pc.angular_distribution_function(atoms, bins=36)
        assert adf.shape == (36,)
        assert angles.shape == (36,)


# ---------------------------------------------------------------------------
# BLDF tests
# ---------------------------------------------------------------------------


class TestBondLengthDistribution:
    """Bond-Length Distribution Function tests."""

    def test_bldf_fcc_single_peak(self):
        """Perfect FCC should have a single sharp peak at the NN distance."""
        atoms = _setup("Al", cutoff=0)
        bldf, r = pc.bond_length_distribution(atoms, bins=200)
        # Peak should be near 2.86 Å (Al FCC nearest-neighbor distance)
        peak_r = r[np.argmax(bldf)]
        assert abs(peak_r - 2.86) < 0.15, f"FCC BLDF peak at {peak_r}"

    def test_bldf_bcc_two_peaks(self):
        """BCC with 14 neighbors should show two peaks (1st + 2nd shell)."""
        atoms = bulk("Fe", cubic=True).repeat(4)
        pc.find_neighbors(atoms, method="cutoff", cutoff=0)
        bldf, r = pc.bond_length_distribution(atoms, bins=200)
        # Find peaks: distinct maxima
        # Fe NN ≈ 2.48, 2nd shell ≈ 2.87
        peak_idx = np.argmax(bldf)
        peak_r = r[peak_idx]
        assert 2.3 < peak_r < 2.6, f"BCC BLDF first peak at {peak_r}"

    def test_bldf_normalised(self):
        """BLDF should integrate approximately to 1."""
        atoms = _setup("Al", cutoff=0)
        bldf, r = pc.bond_length_distribution(atoms, bins=200)
        dr = r[1] - r[0]
        integral = np.sum(bldf) * dr
        assert abs(integral - 1.0) < 0.1, f"BLDF integral = {integral}"

    def test_bldf_stored_in_info(self):
        """Results should be stored in atoms.info."""
        atoms = _setup("Al", cutoff=0)
        pc.bond_length_distribution(atoms, bins=50)
        assert "pyscal_bldf" in atoms.info
        assert "pyscal_bldf_r" in atoms.info
        assert len(atoms.info["pyscal_bldf"]) == 50

    def test_bldf_custom_range(self):
        """Custom rmin/rmax should work."""
        atoms = _setup("Al", cutoff=0)
        bldf, r = pc.bond_length_distribution(atoms, bins=100,
                                               rmin=2.0, rmax=4.0)
        assert r[0] >= 2.0
        assert r[-1] < 4.0

    def test_bldf_voronoi_neighbors(self):
        """BLDF should work with Voronoi neighbors too."""
        atoms = bulk("Cu", cubic=True).repeat(3)
        pc.find_neighbors(atoms, method="voronoi")
        bldf, r = pc.bond_length_distribution(atoms, bins=100)
        assert bldf.shape == (100,)
        peak_r = r[np.argmax(bldf)]
        # Cu FCC NN ~ 2.55
        assert abs(peak_r - 2.55) < 0.2, f"Cu BLDF peak at {peak_r}"


# ---------------------------------------------------------------------------
# Combined tests
# ---------------------------------------------------------------------------


class TestADFBLDFCombined:
    """Test that ADF and BLDF can be computed together."""

    def test_both_on_same_atoms(self):
        """Both functions should work on the same atoms object."""
        atoms = _setup("Al", cutoff=0)
        adf, angles = pc.angular_distribution_function(atoms, bins=90)
        bldf, r = pc.bond_length_distribution(atoms, bins=100)
        assert adf.shape == (90,)
        assert bldf.shape == (100,)

    def test_hcp_adf_distinct_from_fcc(self):
        """HCP ADF should differ from FCC ADF."""
        fcc = bulk("Al", cubic=True).repeat(4)
        pc.find_neighbors(fcc, method="cutoff", cutoff=0)
        adf_fcc, _ = pc.angular_distribution_function(fcc, bins=180)

        hcp = bulk("Mg", "hcp", a=3.21, c=5.21).repeat((4, 4, 3))
        pc.find_neighbors(hcp, method="cutoff", cutoff=0)
        adf_hcp, _ = pc.angular_distribution_function(hcp, bins=180)

        # They should not be identical
        assert not np.allclose(adf_fcc, adf_hcp, atol=0.01)
