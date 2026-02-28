"""Tests for bond orientational order W_l parameters."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3


class TestWigner3j:
    """Test Wigner 3j symbol implementation."""
    
    def test_selection_rules(self):
        """Verify selection rules give zero."""
        from pyscal3.descriptors import _wigner_3j
        
        # m1 + m2 + m3 != 0 should give 0
        assert _wigner_3j(2, 2, 2, 1, 1, 1) == 0.0
        
        # Triangle inequality violation
        assert _wigner_3j(1, 1, 10, 0, 0, 0) == 0.0
    
    def test_known_values(self):
        """Test against known Wigner 3j values."""
        from pyscal3.descriptors import _wigner_3j
        
        # Known value: (2 2 2; 0 0 0) = ±sqrt(2/35) (sign depends on convention)
        w = _wigner_3j(2, 2, 2, 0, 0, 0)
        expected = np.sqrt(2/35)
        assert abs(abs(w) - expected) < 1e-10
        
        # Known value: (4 4 4; 0 0 0) should be non-zero and consistent
        w440 = _wigner_3j(4, 4, 4, 0, 0, 0)
        assert abs(w440) > 0.01  # Should be non-trivial
        
        # Test symmetry: (l l l; m1 m2 m3) should be related by permutations
        w1 = _wigner_3j(4, 4, 4, 2, 1, -3)
        w2 = _wigner_3j(4, 4, 4, 1, 2, -3)
        # Under column permutation, 3j symbol picks up a phase
        assert abs(abs(w1) - abs(w2)) < 1e-10


class TestW4W6:
    """Tests for W_4 and W_6 parameters."""
    
    def test_w4_w6_fcc(self):
        """W_4 and W_6 for FCC should be close to reference values."""
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=[4, 6])
        
        # FCC reference values (normalized):
        # Ŵ_4 ≈ -0.159, Ŵ_6 ≈ -0.013
        w4_mean = np.mean(result["w4"])
        w6_mean = np.mean(result["w6"])
        
        # Check approximate values (may vary with cutoff)
        assert -0.2 < w4_mean < -0.1  # Roughly -0.159
        assert -0.05 < w6_mean < 0.01  # Roughly -0.013
    
    def test_w4_w6_bcc(self):
        """W_4 and W_6 for BCC."""
        atoms = bulk("Fe", "bcc", a=2.87, cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=[4, 6])
        
        # BCC reference: Ŵ_4 ≈ +0.159, Ŵ_6 ≈ +0.013
        w4_mean = np.mean(result["w4"])
        w6_mean = np.mean(result["w6"])
        
        # BCC W_4 is positive (opposite sign from FCC)
        assert 0.1 < w4_mean < 0.2
    
    def test_fcc_vs_bcc_distinguishable(self):
        """FCC and BCC should have distinguishable W_4 values."""
        fcc = bulk("Cu", "fcc", a=3.6, cubic=True).repeat(3)
        bcc = bulk("Fe", "bcc", a=2.87, cubic=True).repeat(3)
        
        pyscal3.find_neighbors(fcc, method="cutoff", cutoff=3.0)
        pyscal3.find_neighbors(bcc, method="cutoff", cutoff=3.0)
        
        w_fcc = pyscal3.bond_orientational_order_w(fcc, l=4)
        w_bcc = pyscal3.bond_orientational_order_w(bcc, l=4)
        
        # FCC has negative W_4, BCC has positive W_4
        assert np.mean(w_fcc["w4"]) < 0
        assert np.mean(w_bcc["w4"]) > 0


class TestNormalization:
    """Tests for normalized vs unnormalized W_l."""
    
    def test_raw_values_stored(self):
        """Raw (unnormalized) values should be stored."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=6, normalized=True)
        
        assert "w6_raw" in result
        assert "w6" in result
        
        # Raw and normalized should differ (normalization changes values)
        assert not np.allclose(result["w6"], result["w6_raw"])
    
    def test_unnormalized_mode(self):
        """Test with normalized=False."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result_norm = pyscal3.bond_orientational_order_w(atoms, l=6, normalized=True)
        result_raw = pyscal3.bond_orientational_order_w(atoms, l=6, normalized=False)
        
        # Raw mode should give same as raw values from normalized mode
        assert np.allclose(result_raw["w6"], result_norm["w6_raw"])


class TestStoredResults:
    """Tests for results stored on atoms."""
    
    def test_stored_in_arrays(self):
        """Results should be stored in atoms.arrays."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        pyscal3.bond_orientational_order_w(atoms, l=[4, 6])
        
        assert "pyscal_w4" in atoms.arrays
        assert "pyscal_w6" in atoms.arrays
        assert "pyscal_w4_raw" in atoms.arrays
        assert "pyscal_w6_raw" in atoms.arrays


class TestUniformCrystals:
    """Test that perfect crystals have uniform W values."""
    
    def test_uniform_w_values(self):
        """All atoms in perfect crystal should have same W value."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=6)
        
        # Standard deviation should be very small
        std = np.std(result["w6"])
        assert std < 1e-6


class TestMultipleL:
    """Tests for computing multiple l values at once."""
    
    def test_multiple_l_values(self):
        """Should compute multiple l values in one call."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=[4, 6, 8])
        
        assert "w4" in result
        assert "w6" in result
        assert "w8" in result
    
    def test_single_l_as_int(self):
        """Should accept single l as int."""
        atoms = bulk("Cu", "fcc", cubic=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=6)
        
        assert "w6" in result


class TestHCP:
    """Test HCP structure."""
    
    def test_hcp_w4_w6(self):
        """HCP should have characteristic W values."""
        atoms = bulk("Ti", "hcp", a=2.95, c=4.68).repeat([4, 4, 2])
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.2)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=[4, 6])
        
        w4_mean = np.mean(result["w4"])
        w6_mean = np.mean(result["w6"])
        
        # HCP reference: Ŵ_4 ≈ +0.134, Ŵ_6 ≈ -0.012
        # Note: may vary with cutoff and c/a ratio
        assert 0.05 < w4_mean < 0.2  # Positive
        assert -0.05 < w6_mean < 0.02


class TestEdgeCases:
    """Edge cases and robustness."""
    
    def test_single_atom(self):
        """Single atom with no neighbors."""
        from ase import Atoms
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
        
        result = pyscal3.bond_orientational_order_w(atoms, l=6)
        
        # With no neighbors, q_lm are all zero, so W_l should be zero
        # (or nan-like if normalization divides by zero, but we handle that)
        assert np.isfinite(result["w6"][0])  # Should not be nan/inf
