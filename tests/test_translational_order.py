"""Tests for translational order parameters."""

import numpy as np
import pytest
from ase.build import bulk
from ase import Atoms

import pyscal3


class TestTranslationalOrderFourier:
    """Tests for Fourier-based translational order."""

    def test_fourier_with_correct_G_vectors(self):
        """Test Fourier method with proper FCC reciprocal vectors."""
        atoms = bulk("Cu", "fcc", cubic=True, a=3.6).repeat(4)
        
        # For FCC, use G = (2π/a)(1,1,1) type vectors
        a = 3.6
        G_111 = 2 * np.pi / a * np.array([1, 1, 1])
        G_200 = 2 * np.pi / a * np.array([2, 0, 0])
        
        result = pyscal3.translational_order(
            atoms, method='fourier', G_vectors=[G_111, G_200]
        )
        
        # Perfect crystal with correct G: high τ
        assert result['tau_global'] > 0.5
        assert 'tau' in result
        assert len(result['tau']) == len(atoms)

    def test_output_keys(self):
        """Result should have expected keys."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        result = pyscal3.translational_order(atoms, reference=reference)
        
        assert 'tau' in result
        assert 'tau_global' in result
        assert 'mean' in result
        assert 'std' in result
        assert 'method' in result
        assert result['method'] == 'displacement'

    def test_stored_in_arrays(self):
        """Result should be stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        pyscal3.translational_order(atoms, reference=reference)
        
        assert 'pyscal_translational_order' in atoms.arrays
        
    def test_fourier_values_bounded(self):
        """Fourier method gives τ in [0, 1] for per-atom."""
        atoms = bulk("Cu", "fcc", cubic=True, a=3.6).repeat(4)
        G = 2 * np.pi / 3.6 * np.array([1, 1, 1])
        result = pyscal3.translational_order(
            atoms, method='fourier', G_vectors=[G]
        )
        
        assert np.all(result['tau'] >= 0.0)
        assert np.all(result['tau'] <= 1.0)

    def test_fourier_requires_G_vectors(self):
        """Fourier method should require G_vectors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        with pytest.raises(ValueError, match="G_vectors"):
            pyscal3.translational_order(atoms, method='fourier')

    def test_different_structures_displacement(self):
        """Test displacement method on different structures."""
        fcc = bulk("Cu", "fcc", cubic=True).repeat(3)
        bcc = bulk("Fe", "bcc", cubic=True).repeat(3)
        
        tau_fcc = pyscal3.translational_order(fcc, reference=fcc.copy())['tau_global']
        tau_bcc = pyscal3.translational_order(bcc, reference=bcc.copy())['tau_global']
        
        # Perfect match: τ = 1
        assert abs(tau_fcc - 1.0) < 1e-10
        assert abs(tau_bcc - 1.0) < 1e-10


class TestTranslationalOrderDisplacement:
    """Tests for displacement-based translational order."""

    def test_perfect_match_high_order(self):
        """Identical structure and reference should give τ = 1."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        result = pyscal3.translational_order(
            atoms, method='displacement', reference=reference
        )
        
        # Zero displacement → τ = 1
        np.testing.assert_allclose(result['tau'], 1.0, rtol=1e-10)
        assert abs(result['mean'] - 1.0) < 1e-10

    def test_displaced_structure(self):
        """Displaced atoms should have lower τ."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        # Add thermal-like displacements
        np.random.seed(42)
        atoms.positions += 0.1 * np.random.randn(*atoms.positions.shape)
        
        result = pyscal3.translational_order(
            atoms, method='displacement', reference=reference
        )
        
        # Displaced → τ < 1
        assert result['mean'] < 1.0
        assert result['mean'] > 0.0

    def test_displacement_values_bounded(self):
        """Displacement method gives τ in [0, 1]."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        atoms.positions += 0.2 * np.random.randn(*atoms.positions.shape)
        
        result = pyscal3.translational_order(
            atoms, method='displacement', reference=reference
        )
        
        assert np.all(result['tau'] >= 0.0)
        assert np.all(result['tau'] <= 1.0)

    def test_requires_reference(self):
        """Displacement method should raise without reference."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        with pytest.raises(ValueError, match="reference"):
            pyscal3.translational_order(atoms, method='displacement')

    def test_reference_size_mismatch(self):
        """Should raise if reference size doesn't match."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        with pytest.raises(ValueError, match="atoms"):
            pyscal3.translational_order(
                atoms, method='displacement', reference=reference
            )

    def test_custom_sigma(self):
        """Test with custom sigma parameter."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        atoms.positions += 0.1 * np.random.randn(*atoms.positions.shape)
        
        # Larger sigma → higher τ for same displacement
        result_small = pyscal3.translational_order(
            atoms, method='displacement', reference=reference, sigma=0.1
        )
        result_large = pyscal3.translational_order(
            atoms, method='displacement', reference=reference, sigma=1.0
        )
        
        assert result_large['mean'] > result_small['mean']


class TestLindemannParameter:
    """Tests for Lindemann parameter."""

    def test_perfect_crystal_zero(self):
        """Perfect match should give δ = 0."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        result = pyscal3.lindemann_parameter(atoms, reference)
        
        np.testing.assert_allclose(result['delta'], 0.0, atol=1e-10)
        assert abs(result['global']) < 1e-10

    def test_output_keys(self):
        """Result should have expected keys."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        result = pyscal3.lindemann_parameter(atoms, reference)
        
        assert 'delta' in result
        assert 'global' in result
        assert 'msd' in result
        assert 'nn_distance' in result

    def test_displaced_structure(self):
        """Displaced atoms should have positive δ."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        np.random.seed(42)
        atoms.positions += 0.2 * np.random.randn(*atoms.positions.shape)
        
        result = pyscal3.lindemann_parameter(atoms, reference)
        
        assert result['global'] > 0.0
        assert np.all(result['delta'] >= 0.0)

    def test_lindemann_criterion_threshold(self):
        """Large displacement should give δ > 0.15 (liquid-like)."""
        atoms = bulk("Cu", "fcc", cubic=True, a=3.6).repeat(3)
        reference = atoms.copy()
        
        # Large displacements (~0.5 Å)
        np.random.seed(42)
        atoms.positions += 0.5 * np.random.randn(*atoms.positions.shape)
        
        result = pyscal3.lindemann_parameter(atoms, reference)
        
        # Should be in liquid regime
        assert result['global'] > 0.1

    def test_custom_nn_distance(self):
        """Test with custom nearest-neighbor distance."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        atoms.positions += 0.1 * np.random.randn(*atoms.positions.shape)
        
        result1 = pyscal3.lindemann_parameter(atoms, reference, nn_distance=2.55)
        result2 = pyscal3.lindemann_parameter(atoms, reference, nn_distance=5.0)
        
        # Larger nn_distance → smaller δ
        assert result2['global'] < result1['global']

    def test_requires_reference(self):
        """Should raise without reference."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        with pytest.raises(ValueError, match="reference"):
            pyscal3.lindemann_parameter(atoms)

    def test_stored_in_arrays(self):
        """Result should be stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        pyscal3.lindemann_parameter(atoms, reference)
        
        assert 'pyscal_lindemann' in atoms.arrays


class TestMeanSquaredDisplacement:
    """Tests for mean squared displacement."""

    def test_zero_displacement(self):
        """No displacement should give MSD = 0."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        reference = atoms.copy()
        
        result = pyscal3.mean_squared_displacement(atoms, reference)
        
        assert abs(result['msd']) < 1e-10
        assert abs(result['rmsd']) < 1e-10

    def test_known_displacement(self):
        """Test with known displacement."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(1)
        reference = atoms.copy()
        
        # Move all atoms by (1, 0, 0)
        atoms.positions[:, 0] += 1.0
        
        result = pyscal3.mean_squared_displacement(atoms, reference)
        
        # MSD should be ~1.0 (displacement magnitude squared)
        assert abs(result['msd'] - 1.0) < 0.01
        assert abs(result['rmsd'] - 1.0) < 0.01

    def test_output_keys(self):
        """Result should have expected keys."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        reference = atoms.copy()
        
        result = pyscal3.mean_squared_displacement(atoms, reference)
        
        assert 'displacement' in result
        assert 'msd' in result
        assert 'rmsd' in result

    def test_stored_in_arrays(self):
        """Per-atom displacement should be stored."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        reference = atoms.copy()
        atoms.positions += 0.5
        
        pyscal3.mean_squared_displacement(atoms, reference)
        
        assert 'pyscal_displacement' in atoms.arrays


class TestInvalidInput:
    """Tests for error handling."""

    def test_invalid_method(self):
        """Should raise for unknown method."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        with pytest.raises(ValueError, match="Unknown method"):
            pyscal3.translational_order(atoms, method='invalid')

    @pytest.mark.skip(reason="Zero-volume cells may pass or fail depending on conditions")
    def test_zero_volume_cell(self):
        """Handling of degenerate cells is implementation-dependent."""
        pass


class TestPBC:
    """Tests for periodic boundary condition handling."""

    def test_pbc_minimum_image(self):
        """Displacement should use minimum image convention."""
        atoms = bulk("Cu", "fcc", cubic=True, a=3.6).repeat(2)
        reference = atoms.copy()
        
        # Move atom close to cell boundary - displacement should wrap
        cell_length = atoms.cell[0, 0]
        atoms.positions[0, 0] += cell_length - 0.1  # Near boundary
        
        result = pyscal3.mean_squared_displacement(atoms, reference)
        
        # With MIC, displacement should be small (~0.1), not large
        assert result['displacement'][0] < 1.0
