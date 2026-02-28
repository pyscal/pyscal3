"""Tests for SOAP (Smooth Overlap of Atomic Positions) descriptors."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3


class TestSOAPBasic:
    """Basic SOAP functionality tests."""
    
    def test_soap_fcc(self):
        """SOAP computes for FCC structure."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        
        assert "power_spectrum" in result
        assert result["power_spectrum"].shape[0] == len(atoms)
        assert result["descriptor_size"] == result["power_spectrum"].shape[1]
    
    def test_soap_bcc(self):
        """SOAP computes for BCC structure."""
        atoms = bulk("Fe", "bcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        
        assert result["power_spectrum"].shape[0] == len(atoms)
        assert np.all(np.isfinite(result["power_spectrum"]))
    
    def test_soap_hcp(self):
        """SOAP computes for HCP structure."""
        atoms = bulk("Mg", "hcp").repeat((2, 2, 2))
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        
        assert result["power_spectrum"].shape[0] == len(atoms)


class TestSOAPDimensions:
    """Test descriptor dimensions."""
    
    def test_descriptor_size(self):
        """Check descriptor size formula: n_max*(n_max+1)/2 * (l_max+1)."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        # Test different n_max, l_max combinations
        for n_max, l_max in [(4, 2), (6, 4), (8, 6)]:
            result = pyscal3.soap(atoms, r_cut=5.0, n_max=n_max, l_max=l_max)
            expected_size = (n_max * (n_max + 1) // 2) * (l_max + 1)
            assert result["descriptor_size"] == expected_size
            assert result["power_spectrum"].shape[1] == expected_size
    
    def test_parameters_stored(self):
        """Check that parameters are returned."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=4.5, n_max=5, l_max=3, sigma=0.4)
        
        assert result["r_cut"] == 4.5
        assert result["n_max"] == 5
        assert result["l_max"] == 3
        assert result["sigma"] == 0.4


class TestSOAPNormalization:
    """Test normalization options."""
    
    def test_normalized_unit_norm(self):
        """Normalized SOAP vectors have unit L2 norm."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3, normalized=True)
        
        norms = np.linalg.norm(result["power_spectrum"], axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-10)
    
    def test_unnormalized(self):
        """Unnormalized SOAP has varying norms."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3, normalized=False)
        
        norms = np.linalg.norm(result["power_spectrum"], axis=1)
        # For uniform crystal, all norms should be similar but not 1
        assert not np.allclose(norms, 1.0)


class TestSOAPRotationalInvariance:
    """Test rotational invariance property."""
    
    def test_rotation_invariance(self):
        """SOAP should be invariant under rotation of the entire structure."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        # Rotation matrix (90 degrees around z)
        R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        
        # Original
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        result1 = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        
        # Rotated
        atoms_rot = atoms.copy()
        atoms_rot.positions = atoms_rot.positions @ R.T
        atoms_rot.cell = atoms_rot.cell @ R.T
        pyscal3.find_neighbors(atoms_rot, method="cutoff", cutoff=6.0)
        result2 = pyscal3.soap(atoms_rot, r_cut=5.0, n_max=4, l_max=3)
        
        # SOAP vectors should be identical (up to numerical precision)
        np.testing.assert_allclose(
            result1["power_spectrum"], 
            result2["power_spectrum"], 
            atol=1e-6
        )


class TestSOAPUniformCrystal:
    """Test behavior in uniform crystals."""
    
    def test_uniform_fcc(self):
        """All atoms in perfect FCC should have same SOAP."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        ps = result["power_spectrum"]
        
        # All rows should be identical
        for i in range(1, len(atoms)):
            np.testing.assert_allclose(ps[0], ps[i], atol=1e-10)
    
    def test_fcc_vs_bcc_different(self):
        """FCC and BCC should have different SOAP signatures."""
        fcc = bulk("Cu", "fcc", cubic=True).repeat(2)
        bcc = bulk("Fe", "bcc", cubic=True).repeat(2)
        
        pyscal3.find_neighbors(fcc, method="cutoff", cutoff=6.0)
        pyscal3.find_neighbors(bcc, method="cutoff", cutoff=6.0)
        
        result_fcc = pyscal3.soap(fcc, r_cut=5.0, n_max=4, l_max=3)
        result_bcc = pyscal3.soap(bcc, r_cut=5.0, n_max=4, l_max=3)
        
        # Mean SOAP vectors should be different
        mean_fcc = np.mean(result_fcc["power_spectrum"], axis=0)
        mean_bcc = np.mean(result_bcc["power_spectrum"], axis=0)
        
        assert not np.allclose(mean_fcc, mean_bcc)


class TestSOAPStorage:
    """Test that results are stored in atoms.arrays."""
    
    def test_stored_in_arrays(self):
        """Results stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        
        assert "pyscal_soap" in atoms.arrays
        assert atoms.arrays["pyscal_soap"].shape[0] == len(atoms)


class TestSOAPParameters:
    """Test different parameter settings."""
    
    def test_different_sigma(self):
        """Different sigma values produce different results."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result1 = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3, sigma=0.3)
        result2 = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3, sigma=0.6)
        
        assert not np.allclose(result1["power_spectrum"], result2["power_spectrum"])
    
    def test_different_cutoff(self):
        """Different cutoff radii produce different radial basis computations."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        # With different r_cut, the radial basis normalization changes
        result1 = pyscal3.soap(atoms, r_cut=3.5, n_max=4, l_max=3, normalized=False)
        result2 = pyscal3.soap(atoms, r_cut=5.5, n_max=4, l_max=3, normalized=False)
        
        # Raw (unnormalized) power spectra should differ due to different radial bases
        # At minimum, the magnitudes should be different
        sum1 = np.sum(result1["power_spectrum"])
        sum2 = np.sum(result2["power_spectrum"])
        assert abs(sum1 - sum2) / max(abs(sum1), abs(sum2)) > 0.01  # At least 1% different


class TestSOAPEdgeCases:
    """Test edge cases."""
    
    def test_single_atom(self):
        """Single atom with no neighbors."""
        from ase import Atoms
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=4, l_max=3)
        
        # Should return zeros when no neighbors
        np.testing.assert_allclose(result["power_spectrum"][0], 0.0)
    
    def test_small_n_max_l_max(self):
        """Small n_max and l_max."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.soap(atoms, r_cut=5.0, n_max=2, l_max=1)
        
        # n_max*(n_max+1)/2 * (l_max+1) = 2*3/2 * 2 = 6
        assert result["descriptor_size"] == 6
        assert result["power_spectrum"].shape == (len(atoms), 6)
