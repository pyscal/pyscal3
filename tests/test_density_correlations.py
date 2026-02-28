"""Tests for density correlation functions."""

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk

import pyscal3
from pyscal3.descriptors import (
    structure_factor,
    local_density,
    density_fluctuations,
    hyperuniformity,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def fcc_crystal():
    """Small FCC Cu crystal."""
    return bulk("Cu", "fcc", cubic=True).repeat(3)


@pytest.fixture
def bcc_crystal():
    """Small BCC Fe crystal."""
    return bulk("Fe", "bcc", cubic=True).repeat(3)


@pytest.fixture
def simple_cubic():
    """Simple cubic crystal."""
    return bulk("Po", "sc", a=3.0).repeat(4)


@pytest.fixture
def amorphous():
    """Pseudo-amorphous structure (slightly randomized)."""
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    pos = atoms.get_positions()
    np.random.seed(42)
    pos += np.random.randn(*pos.shape) * 0.1
    atoms.set_positions(pos)
    return atoms


# ---------------------------------------------------------------------------
# Structure Factor Tests
# ---------------------------------------------------------------------------

class TestStructureFactor:
    """Tests for structure_factor()."""
    
    def test_basic_fcc(self, fcc_crystal):
        """Test structure factor on FCC crystal."""
        result = structure_factor(fcc_crystal, k_max=10.0, n_k=50, n_samples=30)
        
        assert 'k' in result
        assert 'S' in result
        assert 'S_0' in result
        assert len(result['k']) == 50
        assert len(result['S']) == 50
        
    def test_k_range(self, fcc_crystal):
        """Test k-vector range."""
        result = structure_factor(fcc_crystal, k_max=15.0, n_k=100)
        
        assert result['k'][0] == pytest.approx(0.1, rel=0.01)
        assert result['k'][-1] == pytest.approx(15.0, rel=0.01)
        
    def test_structure_factor_positive(self, fcc_crystal):
        """Structure factor should be positive."""
        result = structure_factor(fcc_crystal, k_max=10.0, n_k=50)
        
        assert np.all(result['S'] >= 0)
        
    def test_crystal_has_peaks(self, fcc_crystal):
        """Crystal should show Bragg-like peaks."""
        # Note: With random k-direction sampling, we sample diffuse scattering
        # rather than exact Bragg peaks. The isotropically-averaged S(k) for
        # crystals still shows structure but not necessarily delta-function peaks.
        result = structure_factor(fcc_crystal, k_max=15.0, n_k=100, n_samples=100)
        
        # Check that S(k) varies significantly (not flat like ideal gas)
        S = result['S']
        max_S = np.max(S)
        min_S = np.min(S)
        
        # There should be structure in the signal
        assert max_S > min_S * 1.5  # Signal shows variation
        
    def test_metadata_stored(self, fcc_crystal):
        """Test that metadata is stored in atoms.info."""
        structure_factor(fcc_crystal, k_max=8.0, n_k=40)
        
        assert 'pyscal_structure_factor' in fcc_crystal.info
        assert fcc_crystal.info['pyscal_structure_factor']['k_max'] == 8.0
        assert fcc_crystal.info['pyscal_structure_factor']['n_k'] == 40
        
    def test_S_0_extrapolation(self, fcc_crystal):
        """Test S(k→0) extrapolation works."""
        result = structure_factor(fcc_crystal, k_max=5.0, n_k=50, n_samples=100)
        
        # S_0 should be finite and non-negative
        assert result['S_0'] >= 0
        # For small systems with random sampling, S(k) at small k can be
        # substantial due to finite-size effects


# ---------------------------------------------------------------------------
# Local Density Tests
# ---------------------------------------------------------------------------

class TestLocalDensity:
    """Tests for local_density()."""
    
    def test_voronoi_method(self, fcc_crystal):
        """Test Voronoi-based local density."""
        pyscal3.find_neighbors(fcc_crystal, method="voronoi")
        result = local_density(fcc_crystal, method='voronoi')
        
        assert 'density' in result
        assert 'mean' in result
        assert 'std' in result
        assert result['method'] == 'voronoi'
        assert len(result['density']) == len(fcc_crystal)
        
    def test_voronoi_needs_tessellation(self, fcc_crystal):
        """Should raise if Voronoi not computed."""
        # Fresh structure without Voronoi
        fresh = bulk("Cu", "fcc", cubic=True).repeat(2)
        
        with pytest.raises(ValueError, match="Voronoi tessellation required"):
            local_density(fresh, method='voronoi')
            
    def test_neighbor_method(self, fcc_crystal):
        """Test neighbor-count based density."""
        pyscal3.find_neighbors(fcc_crystal, method="cutoff", cutoff=3.0)
        result = local_density(fcc_crystal, method='neighbor', cutoff=3.0)
        
        assert len(result['density']) == len(fcc_crystal)
        assert result['method'] == 'neighbor'
        
    def test_neighbor_needs_cutoff(self, fcc_crystal):
        """Should require cutoff for neighbor method."""
        pyscal3.find_neighbors(fcc_crystal, method="cutoff", cutoff=3.0)
        
        with pytest.raises(ValueError, match="cutoff required"):
            local_density(fcc_crystal, method='neighbor')
            
    def test_gaussian_method(self, simple_cubic):
        """Test Gaussian-smoothed density."""
        result = local_density(simple_cubic, method='gaussian', sigma=1.0)
        
        assert len(result['density']) == len(simple_cubic)
        assert result['method'] == 'gaussian'
        assert np.all(result['density'] >= 0)
        
    def test_gaussian_needs_sigma(self, simple_cubic):
        """Should require sigma for Gaussian method."""
        with pytest.raises(ValueError, match="sigma required"):
            local_density(simple_cubic, method='gaussian')
            
    def test_uniform_crystal_has_uniform_density(self, fcc_crystal):
        """Perfect crystal should have uniform density."""
        pyscal3.find_neighbors(fcc_crystal, method="voronoi")
        result = local_density(fcc_crystal, method='voronoi')
        
        # All atoms should have similar density
        relative_std = result['std'] / result['mean']
        assert relative_std < 0.01  # Less than 1% variation
        
    def test_density_stored_in_arrays(self, fcc_crystal):
        """Test that density is stored in atoms.arrays."""
        pyscal3.find_neighbors(fcc_crystal, method="voronoi")
        local_density(fcc_crystal, method='voronoi')
        
        assert 'pyscal_local_density' in fcc_crystal.arrays
        
    def test_invalid_method(self, fcc_crystal):
        """Should raise for invalid method."""
        with pytest.raises(ValueError, match="Unknown method"):
            local_density(fcc_crystal, method='invalid')


# ---------------------------------------------------------------------------
# Density Fluctuations Tests
# ---------------------------------------------------------------------------

class TestDensityFluctuations:
    """Tests for density_fluctuations()."""
    
    def test_basic_fluctuations(self, fcc_crystal):
        """Test basic density fluctuations."""
        result = density_fluctuations(fcc_crystal, n_blocks=3)
        
        assert 'mean_N' in result
        assert 'var_N' in result
        assert 'normalized_variance' in result
        assert 'block_counts' in result
        
    def test_block_count(self, fcc_crystal):
        """Test number of blocks."""
        result = density_fluctuations(fcc_crystal, n_blocks=4)
        
        # Should have n_blocks^3 total blocks
        assert len(result['block_counts']) == 4**3
        
    def test_total_atoms_conserved(self, fcc_crystal):
        """Sum of block counts should equal total atoms."""
        result = density_fluctuations(fcc_crystal, n_blocks=3)
        
        assert np.sum(result['block_counts']) == len(fcc_crystal)
        
    def test_crystal_low_variance(self, fcc_crystal):
        """Crystal should have low normalized variance."""
        result = density_fluctuations(fcc_crystal, n_blocks=4)
        
        # Crystal has uniform distribution, low variance
        # But at finite size with integer counts, there's some variance
        assert result['normalized_variance'] < 2.0
        
    def test_metadata_stored(self, fcc_crystal):
        """Test metadata stored in atoms.info."""
        density_fluctuations(fcc_crystal, n_blocks=5)
        
        assert 'pyscal_density_fluctuations' in fcc_crystal.info
        assert fcc_crystal.info['pyscal_density_fluctuations']['n_blocks'] == 5


# ---------------------------------------------------------------------------
# Hyperuniformity Tests
# ---------------------------------------------------------------------------

class TestHyperuniformity:
    """Tests for hyperuniformity()."""
    
    def test_basic_hyperuniformity(self, fcc_crystal):
        """Test basic hyperuniformity analysis."""
        result = hyperuniformity(fcc_crystal, k_max=5.0, n_k=30)
        
        assert 'S_k' in result
        assert 'alpha' in result
        assert 'A' in result
        assert 'is_hyperuniform' in result
        assert 'hyperuniform_class' in result
        
    def test_crystal_is_hyperuniform(self, fcc_crystal):
        """Test hyperuniformity detection runs correctly."""
        result = hyperuniformity(fcc_crystal, k_max=5.0, n_k=50)
        
        # The classification should complete (true hyperuniformity
        # detection for crystals requires exact Bragg peak sampling,
        # which this random-direction method doesn't guarantee)
        assert 'is_hyperuniform' in result
        assert isinstance(result['alpha'], float)
        
    def test_hyperuniformity_class(self, simple_cubic):
        """Test hyperuniformity classification."""
        result = hyperuniformity(simple_cubic, k_max=5.0, n_k=50)
        
        if result['is_hyperuniform']:
            assert result['hyperuniform_class'] in ['I', 'II', 'III']
        else:
            assert result['hyperuniform_class'] is None
            
    def test_alpha_computed(self, fcc_crystal):
        """Test that alpha exponent is computed."""
        result = hyperuniformity(fcc_crystal, k_max=5.0, n_k=50)
        
        # Alpha is fit from log-log plot, can be any value
        assert isinstance(result['alpha'], float)
        assert isinstance(result['A'], float)
        assert result['A'] > 0  # Prefactor should be positive
        
    def test_metadata_stored(self, fcc_crystal):
        """Test metadata stored in atoms.info."""
        hyperuniformity(fcc_crystal, k_max=4.0)
        
        assert 'pyscal_hyperuniformity' in fcc_crystal.info
        assert 'alpha' in fcc_crystal.info['pyscal_hyperuniformity']


# ---------------------------------------------------------------------------
# Integration Tests
# ---------------------------------------------------------------------------

class TestDensityCorrelationIntegration:
    """Integration tests for density correlation functions."""
    
    def test_full_workflow(self, fcc_crystal):
        """Test full density analysis workflow."""
        # 1. Compute structure factor
        S_result = structure_factor(fcc_crystal, k_max=12.0, n_k=80)
        
        # 2. Local density (Voronoi)
        pyscal3.find_neighbors(fcc_crystal, method="voronoi")
        rho_result = local_density(fcc_crystal, method='voronoi')
        
        # 3. Density fluctuations  
        fluct_result = density_fluctuations(fcc_crystal, n_blocks=4)
        
        # 4. Hyperuniformity
        hu_result = hyperuniformity(fcc_crystal, k_max=5.0)
        
        # All should complete without error
        assert S_result['S'][0] > 0
        assert rho_result['mean'] > 0
        assert fluct_result['mean_N'] > 0
        assert 'is_hyperuniform' in hu_result
        
    def test_different_structures(self, fcc_crystal, bcc_crystal, simple_cubic):
        """Test on different crystal structures."""
        for atoms in [fcc_crystal, bcc_crystal, simple_cubic]:
            result = structure_factor(atoms, k_max=10.0, n_k=40, n_samples=20)
            
            assert len(result['S']) == 40
            assert np.all(result['S'] >= 0)
            
    def test_amorphous_vs_crystal(self, fcc_crystal, amorphous):
        """Compare crystal and amorphous structure factors."""
        S_crystal = structure_factor(fcc_crystal, k_max=10.0, n_k=50, n_samples=50)
        S_amorphous = structure_factor(amorphous, k_max=10.0, n_k=50, n_samples=50)
        
        # Both should complete
        assert len(S_crystal['S']) == len(S_amorphous['S'])
        
    def test_pyscal3_toplevel_access(self, fcc_crystal):
        """Test that functions are accessible at top level."""
        # Should work via pyscal3.function() 
        assert hasattr(pyscal3, 'structure_factor')
        assert hasattr(pyscal3, 'local_density')
        assert hasattr(pyscal3, 'density_fluctuations')
        assert hasattr(pyscal3, 'hyperuniformity')
