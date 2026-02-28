"""Tests for Behler-Parrinello symmetry functions."""

import numpy as np
import pytest
from ase.build import bulk, molecule
from ase import Atoms

import pyscal3


class TestG2RadialFunctions:
    """Tests for G2 radial symmetry functions."""
    
    def test_g2_basic(self):
        """G2 should be non-zero for system with neighbors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0, g4_params=[])
        
        assert "G2" in result
        assert result["G2"].shape[0] == len(atoms)
        assert result["G2"].shape[1] > 0
        # All atoms equivalent in perfect crystal, so each row should be identical
        # Check std across atoms (axis=0) is small
        std_per_feature = np.std(result["G2"], axis=0)
        assert np.all(std_per_feature < 1e-6)
    
    def test_g2_cutoff_respected(self):
        """G2 should only count neighbors within cutoff."""
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat(3)
        
        # Find neighbors with large cutoff
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=8.0)
        
        # Compute with smaller ACSF cutoff
        result_small = pyscal3.symmetry_functions(atoms, cutoff=4.0, g4_params=[])
        result_large = pyscal3.symmetry_functions(atoms, cutoff=6.0, g4_params=[])
        
        # Larger cutoff should give larger G2 values (more neighbors)
        assert np.mean(result_large["G2"]) > np.mean(result_small["G2"])
    
    def test_g2_rs_shift(self):
        """G2 with different Rs should detect different distance shells."""
        atoms = bulk("Cu", "fcc", a=3.6, cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        # G2 with Rs at 1NN distance (≈2.55 Å for FCC Cu)
        result = pyscal3.symmetry_functions(
            atoms, cutoff=6.0,
            g2_params=[(0.1, 2.55)],  # Narrow Gaussian at 1NN
            g4_params=[]
        )
        
        assert result["G2"].shape == (len(atoms), 1)
        assert np.all(result["G2"] > 0)


class TestG4AngularFunctions:
    """Tests for G4 angular symmetry functions."""
    
    def test_g4_basic(self):
        """G4 should be non-zero for system with >= 2 neighbors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0)
        
        assert "G4" in result
        assert result["G4"].shape[0] == len(atoms)
        assert np.all(result["G4"] >= 0)  # G4 is always >= 0
    
    def test_g4_lambda_plus_minus(self):
        """λ=+1 prefers 0° angles, λ=-1 prefers 180°."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        
        result_plus = pyscal3.symmetry_functions(
            atoms, cutoff=4.0,
            g2_params=[],
            g4_params=[(0.01, 2, 1)]  # λ=+1
        )
        
        result_minus = pyscal3.symmetry_functions(
            atoms, cutoff=4.0,
            g2_params=[],
            g4_params=[(0.01, 2, -1)]  # λ=-1
        )
        
        # FCC has 60°/90°/120° angles; λ=-1 should give different values
        # Just check they're different
        assert not np.allclose(result_plus["G4"], result_minus["G4"], atol=1e-3)
    
    def test_g4_zeta_angular_resolution(self):
        """Higher ζ gives sharper angular peaks."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        
        result_low = pyscal3.symmetry_functions(
            atoms, cutoff=4.0, g2_params=[],
            g4_params=[(0.01, 1, 1)]  # ζ=1 (broader)
        )
        
        result_high = pyscal3.symmetry_functions(
            atoms, cutoff=4.0, g2_params=[],
            g4_params=[(0.01, 4, 1)]  # ζ=4 (sharper)
        )
        
        # Higher zeta gives smaller values due to 2^(1-zeta) prefactor
        assert np.mean(result_high["G4"]) < np.mean(result_low["G4"])


class TestG5AngularFunctions:
    """Tests for G5 angular functions (without j-k distance)."""
    
    def test_g5_basic(self):
        """G5 should work when specified."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(
            atoms, cutoff=6.0,
            g5_params=[(0.01, 1, 1), (0.01, 2, -1)]
        )
        
        assert "G5" in result
        assert result["G5"].shape[0] == len(atoms)
        assert result["G5"].shape[1] == 2  # 2 params, 1 element
    
    def test_g5_larger_than_g4(self):
        """G5 should be larger than G4 (no j-k cutoff penalty)."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=5.0)
        
        params = [(0.01, 2, 1)]
        result = pyscal3.symmetry_functions(
            atoms, cutoff=5.0,
            g2_params=[],
            g4_params=params,
            g5_params=params
        )
        
        # G5 has no f_c(r_jk) factor, so should be >= G4
        assert np.mean(result["G5"]) >= np.mean(result["G4"])


class TestMultiElement:
    """Tests for multi-element systems."""
    
    def test_binary_alloy(self):
        """Test with binary alloy - should have more features."""
        # Create ordered B2 structure (like CsCl)
        atoms = bulk("NaCl", crystalstructure="rocksalt", a=5.64, cubic=True)
        atoms = atoms.repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0)
        
        # Should have 2 element types
        assert len(result["element_types"]) == 2
        
        # G2: n_params * n_elements = 8 * 2 = 16
        assert result["G2"].shape[1] == 16
        
        # G4: n_params * n_element_pairs = 8 * 3 = 24
        # Element pairs: (Na,Na), (Na,Cl), (Cl,Cl)
        assert result["G4"].shape[1] == 24
    
    def test_element_filter(self):
        """Test filtering to specific elements."""
        atoms = bulk("NaCl", crystalstructure="rocksalt", a=5.64, cubic=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        # Only consider Na-Na interactions
        result = pyscal3.symmetry_functions(
            atoms, cutoff=6.0,
            element_filter=[11],  # Na only
            g4_params=[]
        )
        
        assert result["element_types"] == [11]
        assert result["G2"].shape[1] == 8  # 8 params * 1 element


class TestOutputFormat:
    """Tests for output format and storage."""
    
    def test_features_concatenation(self):
        """Features should be concatenation of G2, G4, G5."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(
            atoms, cutoff=6.0,
            g2_params=[(0.01, 0.0)],
            g4_params=[(0.01, 1, 1)],
            g5_params=[(0.01, 2, 1)]
        )
        
        n_g2 = result["G2"].shape[1]
        n_g4 = result["G4"].shape[1]
        n_g5 = result["G5"].shape[1]
        
        assert result["features"].shape[1] == n_g2 + n_g4 + n_g5
        assert np.allclose(result["features"][:, :n_g2], result["G2"])
    
    def test_stored_on_atoms(self):
        """Results should be stored in atoms.arrays and atoms.info."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        pyscal3.symmetry_functions(atoms, cutoff=6.0)
        
        assert "pyscal_acsf" in atoms.arrays
        assert "pyscal_acsf_params" in atoms.info
        
        params = atoms.info["pyscal_acsf_params"]
        assert params["cutoff"] == 6.0
        assert "g2_params" in params
        assert "n_features" in params
    
    def test_feature_names(self):
        """Feature names should match features."""
        atoms = bulk("Cu", "fcc", cubic=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0)
        
        assert len(result["feature_names"]) == result["features"].shape[1]
        assert all("G2" in name or "G4" in name for name in result["feature_names"])


class TestEdgeCases:
    """Edge cases and error handling."""
    
    def test_isolated_atom(self):
        """Isolated atom with no neighbors should get zero features."""
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0)
        
        assert result["G2"].shape == (1, 8)
        assert np.all(result["G2"] == 0)
        assert np.all(result["G4"] == 0)
    
    def test_dimer(self):
        """Dimer should have non-zero G2 but zero G4."""
        atoms = Atoms("Cu2", positions=[[0, 0, 0], [2.5, 0, 0]], 
                      cell=[10, 10, 10], pbc=False)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0)
        
        assert np.all(result["G2"] > 0)  # Each atom sees one neighbor
        assert np.all(result["G4"] == 0)  # Need >= 2 neighbors for G4
    
    def test_no_g4_params(self):
        """Should work with empty G4 params."""
        atoms = bulk("Cu", "fcc", cubic=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=6.0)
        
        result = pyscal3.symmetry_functions(atoms, cutoff=6.0, g4_params=[])
        
        assert "G2" in result
        assert "G4" not in result  # Empty params = no key
        assert result["features"].shape[1] == result["G2"].shape[1]


class TestStructureDiscrimination:
    """Test that ACSFs can distinguish different structures."""
    
    def test_fcc_vs_bcc(self):
        """FCC and BCC should have different ACSF fingerprints."""
        fcc = bulk("Cu", "fcc", a=3.6, cubic=True).repeat(2)
        bcc = bulk("Cu", "bcc", a=2.87, cubic=True).repeat(2)
        
        pyscal3.find_neighbors(fcc, method="cutoff", cutoff=6.0)
        pyscal3.find_neighbors(bcc, method="cutoff", cutoff=6.0)
        
        result_fcc = pyscal3.symmetry_functions(fcc, cutoff=6.0)
        result_bcc = pyscal3.symmetry_functions(bcc, cutoff=6.0)
        
        # Mean fingerprints should be different
        mean_fcc = np.mean(result_fcc["features"], axis=0)
        mean_bcc = np.mean(result_bcc["features"], axis=0)
        
        # Should be measurably different (not identical)
        relative_diff = np.linalg.norm(mean_fcc - mean_bcc) / np.linalg.norm(mean_fcc)
        assert relative_diff > 0.01  # At least 1% different
