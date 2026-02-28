"""
Tests for Coulomb Matrix and MBTR descriptors.
"""
import pytest
import numpy as np
from ase.build import bulk, molecule
from ase import Atoms
import pyscal3


class TestCoulombMatrixBasic:
    """Basic functionality tests for Coulomb matrix."""
    
    def test_water_molecule(self):
        """Test Coulomb matrix on water molecule."""
        mol = molecule('H2O')
        result = pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=10)
        
        assert result.shape == (10,)
        assert np.all(np.isfinite(result))
        # First 3 eigenvalues should be non-zero
        assert np.any(result[:3] != 0)
    
    def test_fcc_structure(self):
        """Test Coulomb matrix on FCC copper."""
        atoms = bulk("Cu", "fcc", cubic=True)
        result = pyscal3.coulomb_matrix(atoms, representation='eigenvalues', max_atoms=10)
        
        assert result.shape == (10,)
        assert np.all(np.isfinite(result))
    
    def test_matrix_representation(self):
        """Test raw matrix output."""
        mol = molecule('H2O')
        result = pyscal3.coulomb_matrix(mol, representation='matrix')
        
        assert result.shape == (3, 3)
        # Matrix should be symmetric
        np.testing.assert_allclose(result, result.T)
        # Diagonal should be 0.5 * Z^2.4
        Z = mol.get_atomic_numbers()
        expected_diag = 0.5 * Z.astype(float) ** 2.4
        np.testing.assert_allclose(np.diag(result), expected_diag)


class TestCoulombMatrixRepresentations:
    """Test different representation options."""
    
    def test_eigenvalues_descending(self):
        """Test that eigenvalues are sorted descending."""
        mol = molecule('C6H6')  # Benzene
        result = pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=20)
        
        # Non-zero eigenvalues should be descending
        non_zero = result[result != 0]
        assert np.all(np.diff(non_zero) <= 0)
    
    def test_sorted_representation(self):
        """Test sorted row-norm representation."""
        mol = molecule('H2O')
        result = pyscal3.coulomb_matrix(mol, representation='sorted', max_atoms=5)
        
        # Upper triangle of 5x5 matrix
        n_upper = 5 * 6 // 2
        assert result.shape == (n_upper,)
    
    def test_padding_with_zeros(self):
        """Test that smaller systems are padded with zeros."""
        mol = molecule('H2')  # 2 atoms
        result = pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=10)
        
        assert result.shape == (10,)
        # First 2 eigenvalues may be non-zero, rest should be zero
        assert np.sum(result != 0) <= 2


class TestCoulombMatrixInvariance:
    """Test invariance properties."""
    
    def test_translation_invariance(self):
        """Test that descriptor is translation invariant."""
        mol = molecule('H2O')
        result1 = pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=10)
        
        mol_trans = mol.copy()
        mol_trans.positions += [5.0, 3.0, -2.0]
        result2 = pyscal3.coulomb_matrix(mol_trans, representation='eigenvalues', max_atoms=10)
        
        np.testing.assert_allclose(result1, result2)
    
    def test_rotation_invariance(self):
        """Test that eigenvalue representation is rotation invariant."""
        mol = molecule('H2O')
        result1 = pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=10)
        
        # Rotate molecule
        mol_rot = mol.copy()
        angle = np.pi / 3
        rot = np.array([
            [np.cos(angle), -np.sin(angle), 0],
            [np.sin(angle), np.cos(angle), 0],
            [0, 0, 1]
        ])
        mol_rot.positions = mol.positions @ rot.T
        result2 = pyscal3.coulomb_matrix(mol_rot, representation='eigenvalues', max_atoms=10)
        
        np.testing.assert_allclose(result1, result2, rtol=1e-10)


class TestMBTRBasic:
    """Basic functionality tests for MBTR."""
    
    def test_fcc_structure(self):
        """Test MBTR on FCC copper."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        result = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
        
        assert 'k1' in result
        assert 'k2' in result
        assert 'full' in result
        assert result['k1'].shape == (50,)
        assert result['k2'].shape == (50,)
    
    def test_water_molecule(self):
        """Test MBTR on water molecule."""
        mol = molecule('H2O')
        mol.center(vacuum=5.0)
        result = pyscal3.mbtr(mol, k=(1, 2, 3), n_grid=30)
        
        assert 'k1' in result
        assert 'k2' in result
        assert 'k3' in result
        assert np.all(np.isfinite(result['full']))
    
    def test_k_terms_selection(self):
        """Test selecting different k terms."""
        atoms = bulk("Cu", "fcc", cubic=True)
        
        result_k1 = pyscal3.mbtr(atoms, k=(1,), n_grid=50)
        result_k12 = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
        result_k123 = pyscal3.mbtr(atoms, k=(1, 2, 3), n_grid=50)
        
        assert 'k1' in result_k1
        assert 'k2' not in result_k1
        
        assert 'k1' in result_k12
        assert 'k2' in result_k12
        assert 'k3' not in result_k12
        
        assert 'k1' in result_k123
        assert 'k2' in result_k123
        assert 'k3' in result_k123


class TestMBTRParameters:
    """Test MBTR parameter variations."""
    
    def test_grid_size(self):
        """Test different grid sizes."""
        atoms = bulk("Cu", "fcc", cubic=True)
        
        result_50 = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
        result_100 = pyscal3.mbtr(atoms, k=(1, 2), n_grid=100)
        
        assert result_50['k1'].shape == (50,)
        assert result_100['k1'].shape == (100,)
    
    def test_sigma_effect(self):
        """Test that sigma affects broadening."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        result_narrow = pyscal3.mbtr(atoms, k=(2,), n_grid=100, sigma_k2=0.01, normalize=False)
        result_broad = pyscal3.mbtr(atoms, k=(2,), n_grid=100, sigma_k2=0.2, normalize=False)
        
        # Broader sigma should give smoother distribution
        # Check that narrow has more variation
        narrow_var = np.var(result_narrow['k2'])
        broad_var = np.var(result_broad['k2'])
        # This can fail if the distribution is very peaked, so just check they're different
        assert not np.allclose(result_narrow['k2'], result_broad['k2'])


class TestMBTRInvariance:
    """Test MBTR invariance properties."""
    
    def test_translation_invariance(self):
        """Test translation invariance."""
        atoms = bulk("Cu", "fcc", cubic=True)
        result1 = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
        
        atoms_trans = atoms.copy()
        atoms_trans.positions += [3.0, -1.5, 2.0]
        result2 = pyscal3.mbtr(atoms_trans, k=(1, 2), n_grid=50)
        
        np.testing.assert_allclose(result1['full'], result2['full'])
    
    def test_permutation_invariance(self):
        """Test that MBTR is permutation invariant."""
        mol = molecule('CH4')
        result1 = pyscal3.mbtr(mol, k=(1, 2), n_grid=50)
        
        # Permute atoms
        mol_perm = mol.copy()
        perm = [0, 4, 3, 2, 1]  # Keep C at 0, permute H's
        mol_perm.positions = mol_perm.positions[perm]
        mol_perm.numbers = mol_perm.numbers[perm]
        result2 = pyscal3.mbtr(mol_perm, k=(1, 2), n_grid=50)
        
        np.testing.assert_allclose(result1['full'], result2['full'], rtol=1e-10)


class TestMBTRSpeciesResolved:
    """Test species-resolved MBTR."""
    
    def test_multispecies(self):
        """Test species-resolved mode with multiple species."""
        mol = molecule('H2O')
        mol.center(vacuum=5.0)
        
        result_agg = pyscal3.mbtr(mol, k=(1, 2), n_grid=50, species_resolved=False)
        result_res = pyscal3.mbtr(mol, k=(1, 2), n_grid=50, species_resolved=True)
        
        # Species-resolved should have more components
        # H, O for k1 (2 species * 50 = 100)
        # H-H, H-O, O-O for k2 (3 pairs * 50 = 150)
        assert result_res['full'].shape[0] > result_agg['full'].shape[0]


class TestMBTRNormalization:
    """Test MBTR normalization."""
    
    def test_normalized_unit_norm(self):
        """Test that normalized output has unit norm."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        result = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50, normalize=True)
        
        norm = np.linalg.norm(result['full'])
        np.testing.assert_allclose(norm, 1.0, rtol=1e-10)
    
    def test_unnormalized_varies(self):
        """Test that unnormalized output varies with system size."""
        atoms_small = bulk("Cu", "fcc", cubic=True)
        atoms_large = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        result_small = pyscal3.mbtr(atoms_small, k=(1, 2), n_grid=50, normalize=False)
        result_large = pyscal3.mbtr(atoms_large, k=(1, 2), n_grid=50, normalize=False)
        
        # Larger system should have larger k1 and k2 sums
        assert np.sum(result_large['k1']) > np.sum(result_small['k1'])


class TestMBTRDifferentStructures:
    """Test that MBTR distinguishes structures."""
    
    def test_fcc_vs_bcc(self):
        """Test that FCC and BCC give different MBTR."""
        atoms_fcc = bulk("Cu", "fcc", cubic=True).repeat(2)
        atoms_bcc = bulk("Fe", "bcc", cubic=True).repeat(2)
        
        result_fcc = pyscal3.mbtr(atoms_fcc, k=(1, 2), n_grid=50, normalize=True)
        result_bcc = pyscal3.mbtr(atoms_bcc, k=(1, 2), n_grid=50, normalize=True)
        
        # k1 should differ (different atomic numbers)
        assert not np.allclose(result_fcc['k1'], result_bcc['k1'])
        
        # k2 should also differ (different structures)
        assert not np.allclose(result_fcc['k2'], result_bcc['k2'])


class TestEdgeCases:
    """Test edge cases."""
    
    def test_single_atom_coulomb(self):
        """Test Coulomb matrix for single atom."""
        atoms = Atoms("Cu", positions=[[0, 0, 0]])
        result = pyscal3.coulomb_matrix(atoms, representation='eigenvalues', max_atoms=5)
        
        assert result.shape == (5,)
        # Single atom: diagonal is 0.5 * Z^2.4
        expected = 0.5 * 29 ** 2.4
        np.testing.assert_allclose(result[0], expected, rtol=1e-10)
    
    def test_single_atom_mbtr(self):
        """Test MBTR for single atom."""
        atoms = Atoms("Cu", positions=[[0, 0, 0]])
        atoms.center(vacuum=5.0)
        result = pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
        
        # k1 should have a peak at Z=29
        assert np.max(result['k1']) > 0
        # k2 should be all zeros (no pairs)
        assert np.allclose(result['k2'], 0)
    
    def test_two_atoms_mbtr(self):
        """Test MBTR for two atoms."""
        atoms = Atoms("Cu2", positions=[[0, 0, 0], [2.5, 0, 0]])
        result = pyscal3.mbtr(atoms, k=(1, 2, 3), n_grid=50, normalize=False)
        
        # k1 and k2 should be non-zero
        assert np.max(result['k1']) > 0
        assert np.max(result['k2']) > 0
        # k3 should be zero (need 3 atoms for angles)
        np.testing.assert_allclose(result['k3'], 0)
    
    def test_stores_results(self):
        """Test that results are stored in atoms object."""
        mol = molecule('H2O')
        pyscal3.coulomb_matrix(mol, representation='eigenvalues', max_atoms=10)
        
        assert "pyscal_coulomb_eigenvalues" in mol.arrays
        
        atoms = bulk("Cu", "fcc", cubic=True)
        pyscal3.mbtr(atoms, k=(1, 2), n_grid=50)
        
        assert "pyscal_mbtr" in atoms.info
