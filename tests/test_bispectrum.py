"""Tests for bispectrum components (SNAP-style descriptors)."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3


class TestBispectrumBasic:
    """Basic bispectrum functionality tests."""
    
    def test_bispectrum_fcc(self):
        """Bispectrum computes for FCC structure."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        assert "bispectrum" in result
        assert result["bispectrum"].shape[0] == len(atoms)
        assert result["descriptor_size"] == result["bispectrum"].shape[1]
    
    def test_bispectrum_bcc(self):
        """Bispectrum computes for BCC structure."""
        atoms = bulk("Fe", "bcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        assert result["bispectrum"].shape[0] == len(atoms)
        assert np.all(np.isfinite(result["bispectrum"]))
    
    def test_bispectrum_hcp(self):
        """Bispectrum computes for HCP structure."""
        atoms = bulk("Mg", "hcp").repeat((2, 2, 2))
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        assert result["bispectrum"].shape[0] == len(atoms)


class TestBispectrumDimensions:
    """Test descriptor dimensions."""
    
    def test_descriptor_size_increases_with_jmax(self):
        """More components with higher j_max."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result1 = pyscal3.bispectrum(atoms, j_max=1, r_cut=5.0)
        result2 = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        assert result2["descriptor_size"] > result1["descriptor_size"]
    
    def test_indices_returned(self):
        """Check that indices are returned."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        assert "indices" in result
        assert len(result["indices"]) == result["descriptor_size"]


class TestBispectrumNormalization:
    """Test normalization options."""
    
    def test_normalized_unit_norm(self):
        """Normalized bispectrum vectors have unit L2 norm."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0, normalized=True)
        
        norms = np.linalg.norm(result["bispectrum"], axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-10)
    
    def test_unnormalized(self):
        """Unnormalized bispectrum has varying norms."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0, normalized=False)
        
        norms = np.linalg.norm(result["bispectrum"], axis=1)
        # For uniform crystal, all norms should be similar but not 1
        assert not np.allclose(norms, 1.0)


class TestBispectrumRotationalInvariance:
    """Test rotational invariance property."""
    
    def test_rotation_invariance(self):
        """Bispectrum should be invariant under rotation."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        # Rotation matrix (90 degrees around z)
        R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        
        # Original
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        result1 = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        # Rotated
        atoms_rot = atoms.copy()
        atoms_rot.positions = atoms_rot.positions @ R.T
        atoms_rot.cell = atoms_rot.cell @ R.T
        pyscal3.find_neighbors(atoms_rot, method="cutoff", cutoff=6.0)
        result2 = pyscal3.bispectrum(atoms_rot, j_max=2, r_cut=5.0)
        
        # Bispectrum vectors should be similar (allowing for numerical precision)
        np.testing.assert_allclose(
            result1["bispectrum"], 
            result2["bispectrum"], 
            atol=0.1  # Looser tolerance for bispectrum
        )


class TestBispectrumUniformCrystal:
    """Test behavior in uniform crystals."""
    
    def test_uniform_fcc(self):
        """All atoms in perfect FCC should have same bispectrum."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        bs = result["bispectrum"]
        
        # All rows should be similar
        mean_row = np.mean(bs, axis=0)
        for i in range(len(atoms)):
            np.testing.assert_allclose(bs[i], mean_row, atol=1e-6)
    
    def test_fcc_vs_bcc_different(self):
        """FCC and BCC should have different bispectrum signatures."""
        fcc = bulk("Cu", "fcc", cubic=True).repeat(2)
        bcc = bulk("Fe", "bcc", cubic=True).repeat(2)
        
        pyscal3.find_neighbors(fcc, method="cutoff", cutoff=6.0)
        pyscal3.find_neighbors(bcc, method="cutoff", cutoff=6.0)
        
        # Use unnormalized to see the actual differences
        result_fcc = pyscal3.bispectrum(fcc, j_max=2, r_cut=5.0, normalized=False)
        result_bcc = pyscal3.bispectrum(bcc, j_max=2, r_cut=5.0, normalized=False)
        
        # Mean bispectrum vectors should have different sums (magnitudes)
        sum_fcc = np.sum(np.abs(result_fcc["bispectrum"]))
        sum_bcc = np.sum(np.abs(result_bcc["bispectrum"]))
        
        # They should have different magnitudes due to different coordination
        # FCC has 12 nearest neighbors, BCC has 8
        assert not np.isclose(sum_fcc, sum_bcc, rtol=0.1)


class TestBispectrumStorage:
    """Test that results are stored in atoms.arrays."""
    
    def test_stored_in_arrays(self):
        """Results stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        assert "pyscal_bispectrum" in atoms.arrays
        assert atoms.arrays["pyscal_bispectrum"].shape[0] == len(atoms)


class TestBispectrumEdgeCases:
    """Test edge cases."""
    
    def test_single_atom(self):
        """Single atom with no neighbors."""
        from ase import Atoms
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bispectrum(atoms, j_max=2, r_cut=5.0)
        
        # Should return zeros when no neighbors
        np.testing.assert_allclose(result["bispectrum"][0], 0.0)
    
    def test_j_max_1(self):
        """Small j_max value."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.bispectrum(atoms, j_max=1, r_cut=5.0)
        
        assert result["descriptor_size"] > 0
        assert result["bispectrum"].shape == (len(atoms), result["descriptor_size"])
