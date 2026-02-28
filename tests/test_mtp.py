"""Tests for MTP (Moment Tensor Potentials) descriptors."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3


class TestMTPBasic:
    """Basic MTP functionality tests."""
    
    def test_mtp_fcc(self):
        """MTP computes for FCC structure."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        assert "descriptors" in result
        assert result["descriptors"].shape[0] == len(atoms)
        assert result["descriptor_size"] == result["descriptors"].shape[1]
    
    def test_mtp_bcc(self):
        """MTP computes for BCC structure."""
        atoms = bulk("Fe", "bcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        assert result["descriptors"].shape[0] == len(atoms)
        assert np.all(np.isfinite(result["descriptors"]))
    
    def test_mtp_hcp(self):
        """MTP computes for HCP structure."""
        atoms = bulk("Mg", "hcp").repeat((2, 2, 2))
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        assert result["descriptors"].shape[0] == len(atoms)


class TestMTPLevel:
    """Test level parameter."""
    
    def test_level_increases_descriptors(self):
        """Higher level means more descriptors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result6 = pyscal3.moment_tensor_descriptors(atoms, level=6)
        result8 = pyscal3.moment_tensor_descriptors(atoms, level=8)
        result10 = pyscal3.moment_tensor_descriptors(atoms, level=10)
        
        assert result8["descriptor_size"] >= result6["descriptor_size"]
        assert result10["descriptor_size"] >= result8["descriptor_size"]
    
    def test_basis_specs_returned(self):
        """Basis function specifications are returned."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        assert "basis_specs" in result
        assert len(result["basis_specs"]) == result["descriptor_size"]


class TestMTPNormalization:
    """Test normalization options."""
    
    def test_normalized_unit_norm(self):
        """Normalized MTP vectors have unit L2 norm."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8, normalized=True)
        
        norms = np.linalg.norm(result["descriptors"], axis=1)
        np.testing.assert_allclose(norms, 1.0, atol=1e-10)
    
    def test_unnormalized(self):
        """Unnormalized MTP has varying norms."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8, normalized=False)
        
        norms = np.linalg.norm(result["descriptors"], axis=1)
        assert not np.allclose(norms, 1.0)


class TestMTPRotationalInvariance:
    """Test rotational invariance property."""
    
    def test_rotation_invariance(self):
        """MTP should be invariant under rotation."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        # Rotation matrix (90 degrees around z)
        R = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
        
        # Original
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        result1 = pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        # Rotated
        atoms_rot = atoms.copy()
        atoms_rot.positions = atoms_rot.positions @ R.T
        atoms_rot.cell = atoms_rot.cell @ R.T
        pyscal3.find_neighbors(atoms_rot, method="cutoff", cutoff=6.0)
        result2 = pyscal3.moment_tensor_descriptors(atoms_rot, level=8)
        
        # MTP descriptors should be similar (allowing for some numerical tolerance)
        np.testing.assert_allclose(
            result1["descriptors"], 
            result2["descriptors"], 
            atol=0.05
        )


class TestMTPUniformCrystal:
    """Test behavior in uniform crystals."""
    
    def test_uniform_fcc(self):
        """All atoms in perfect FCC should have same MTP descriptors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8)
        desc = result["descriptors"]
        
        # All rows should be similar
        mean_row = np.mean(desc, axis=0)
        for i in range(len(atoms)):
            np.testing.assert_allclose(desc[i], mean_row, atol=1e-6)
    
    def test_fcc_vs_bcc_different(self):
        """FCC and BCC should have different MTP signatures."""
        fcc = bulk("Cu", "fcc", cubic=True).repeat(2)
        bcc = bulk("Fe", "bcc", cubic=True).repeat(2)
        
        pyscal3.find_neighbors(fcc, method="cutoff", cutoff=6.0)
        pyscal3.find_neighbors(bcc, method="cutoff", cutoff=6.0)
        
        result_fcc = pyscal3.moment_tensor_descriptors(fcc, level=8, normalized=False)
        result_bcc = pyscal3.moment_tensor_descriptors(bcc, level=8, normalized=False)
        
        # Mean descriptors should be different
        sum_fcc = np.sum(np.abs(result_fcc["descriptors"]))
        sum_bcc = np.sum(np.abs(result_bcc["descriptors"]))
        
        assert not np.isclose(sum_fcc, sum_bcc, rtol=0.1)


class TestMTPStorage:
    """Test that results are stored in atoms.arrays."""
    
    def test_stored_in_arrays(self):
        """Results stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        assert "pyscal_mtp" in atoms.arrays
        assert atoms.arrays["pyscal_mtp"].shape[0] == len(atoms)


class TestMTPEdgeCases:
    """Test edge cases."""
    
    def test_single_atom(self):
        """Single atom with no neighbors."""
        from ase import Atoms
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.moment_tensor_descriptors(atoms, level=8)
        
        # Should return zeros when no neighbors
        np.testing.assert_allclose(result["descriptors"][0], 0.0)
    
    def test_n_radial_affects_descriptors(self):
        """More radial functions = more descriptors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result2 = pyscal3.moment_tensor_descriptors(atoms, level=8, n_radial=2)
        result4 = pyscal3.moment_tensor_descriptors(atoms, level=8, n_radial=4)
        
        assert result4["descriptor_size"] >= result2["descriptor_size"]
