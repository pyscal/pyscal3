"""Tests for Wigner-Seitz defect analysis."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3


class TestWignerSeitzBasics:
    """Basic functionality tests."""
    
    def test_perfect_crystal_no_defects(self):
        """Perfect crystal should show no defects."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 0
        assert result["interstitial_count"] == 0
        assert len(result["vacancy_indices"]) == 0
        assert len(result["interstitial_sites"]) == 0
        # All sites occupied exactly once
        assert np.all(result["occupancy"] == 1)
    
    def test_single_vacancy(self):
        """Removing one atom should create one vacancy."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        del current[0]  # Remove first atom
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 1
        assert result["interstitial_count"] == 0
        assert len(result["vacancy_indices"]) == 1
        assert 0 in result["vacancy_indices"]  # Site 0 is vacant
    
    def test_multiple_vacancies(self):
        """Multiple deletions create multiple vacancies."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        # Delete atoms 0, 5, 10
        del current[10]
        del current[5]
        del current[0]
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 3
        assert result["interstitial_count"] == 0
    
    def test_single_interstitial(self):
        """Adding an atom at a site creates an interstitial."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        # Add interstitial near site 0
        pos0 = current[0].position
        current.append("Cu")
        current.positions[-1] = pos0 + [0.1, 0.1, 0.1]  # Offset slightly
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        # Site 0 now has 2 atoms
        assert result["vacancy_count"] == 0
        assert result["interstitial_count"] == 1
        assert 0 in result["interstitial_sites"]
    
    def test_frenkel_pair(self):
        """Move atom from one site to another - one vacancy, one interstitial."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        # Move atom 0 to near atom 10
        current.positions[0] = current.positions[10] + [0.1, 0.0, 0.0]
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 1
        assert result["interstitial_count"] == 1


class TestPerTypeOccupancies:
    """Tests for per-type occupancy tracking."""
    
    def test_binary_alloy_perfect(self):
        """Binary alloy with no defects."""
        # Create NaCl structure
        ref = bulk("NaCl", crystalstructure="rocksalt", a=5.64, cubic=True)
        current = ref.copy()
        
        result = pyscal3.wigner_seitz_analysis(current, ref, per_type_occupancies=True)
        
        assert "occupancy_by_type" in result
        assert "Na" in result["occupancy_by_type"]
        assert "Cl" in result["occupancy_by_type"]
        assert result["vacancy_count"] == 0
    
    def test_antisite_detection(self):
        """Track type occupancies for antisite detection."""
        ref = bulk("NaCl", crystalstructure="rocksalt", a=5.64, cubic=True)
        current = ref.copy()
        
        # Swap Na and Cl at adjacent sites (create antisite pair)
        na_idx = [i for i, s in enumerate(ref.get_chemical_symbols()) if s == "Na"][0]
        cl_idx = [i for i, s in enumerate(ref.get_chemical_symbols()) if s == "Cl"][0]
        
        current.symbols[na_idx] = "Cl"
        current.symbols[cl_idx] = "Na"
        
        result = pyscal3.wigner_seitz_analysis(current, ref, per_type_occupancies=True)
        
        # Check type occupancies changed
        na_occ = result["occupancy_by_type"]["Na"]
        cl_occ = result["occupancy_by_type"]["Cl"]
        
        # Original Na site now has Cl
        assert cl_occ[na_idx] == 1
        assert na_occ[na_idx] == 0
        
        # Original Cl site now has Na
        assert na_occ[cl_idx] == 1
        assert cl_occ[cl_idx] == 0


class TestAffineMapping:
    """Tests for affine mapping with cell distortion."""
    
    def test_expanded_cell_no_mapping(self):
        """Expanded cell without mapping may detect false defects."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        
        # Expand cell by 5%
        current.set_cell(current.get_cell() * 1.05, scale_atoms=True)
        
        # Without mapping, atoms drift from reference sites
        result_no_map = pyscal3.wigner_seitz_analysis(current, ref, affine_mapping="none")
        
        # With mapping, atoms should map back perfectly
        result_mapped = pyscal3.wigner_seitz_analysis(current, ref, affine_mapping="to_reference")
        
        # Mapped should show no defects
        assert result_mapped["vacancy_count"] == 0
        assert result_mapped["interstitial_count"] == 0
    
    def test_invalid_affine_mapping(self):
        """Invalid affine_mapping should raise."""
        ref = bulk("Cu", "fcc", cubic=True)
        current = ref.copy()
        
        with pytest.raises(ValueError, match="affine_mapping"):
            pyscal3.wigner_seitz_analysis(current, ref, affine_mapping="invalid")


class TestStoredResults:
    """Tests for results stored on atoms object."""
    
    def test_results_stored_on_atoms(self):
        """Results should be stored in atoms.arrays and atoms.info."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(2)
        current = ref.copy()
        del current[0]
        
        pyscal3.wigner_seitz_analysis(current, ref)
        
        assert "pyscal_ws_site_index" in current.arrays
        assert "pyscal_ws_occupancy" in current.arrays
        assert "pyscal_ws_vacancy_count" in current.info
        assert "pyscal_ws_interstitial_count" in current.info
        
        assert current.info["pyscal_ws_vacancy_count"] == 1
        assert current.info["pyscal_ws_interstitial_count"] == 0


class TestIdentifyDefectAtoms:
    """Tests for identify_defect_atoms convenience function."""
    
    def test_masks_correct(self):
        """Check that atom masks are correct."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(2)
        current = ref.copy()
        
        # Move atom 0 to near atom 1 (creates vacancy + interstitial)
        current.positions[0] = current.positions[1] + [0.1, 0.0, 0.0]
        
        result = pyscal3.identify_defect_atoms(current, ref)
        
        assert "perfect_mask" in result
        assert "interstitial_mask" in result
        assert "vacancy_positions" in result
        
        # Should have some perfect atoms and some at shared sites
        assert np.sum(result["perfect_mask"]) > 0
        assert np.sum(result["interstitial_mask"]) >= 2  # At least atoms 0 and 1
        
        # Should have one vacancy position
        assert len(result["vacancy_positions"]) == 1
    
    def test_defect_summary(self):
        """Test defect summary string."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(2)
        current = ref.copy()
        del current[0]
        
        result = pyscal3.identify_defect_atoms(current, ref)
        
        assert "1 vacancies" in result["defect_summary"]
    
    def test_no_defects_summary(self):
        """Perfect crystal shows no defects in summary."""
        ref = bulk("Cu", "fcc", cubic=True)
        current = ref.copy()
        
        result = pyscal3.identify_defect_atoms(current, ref)
        
        assert result["defect_summary"] == "No defects"


class TestPBC:
    """Tests for periodic boundary condition handling."""
    
    def test_atom_wrapped_across_boundary(self):
        """Atom wrapped to other side of cell should still be assigned correctly."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        
        # Move atom 0 across periodic boundary
        cell = current.get_cell()
        # Shift by nearly a full cell length (should wrap)
        current.positions[0] = current.positions[0] + cell[0] * 0.99
        
        # Wrap positions back into cell
        current.wrap()
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        # Should still detect no defects (atom near same site via PBC)
        assert result["vacancy_count"] == 0
        assert result["interstitial_count"] == 0
    
    def test_non_cubic_cell(self):
        """Test with HCP cell."""
        ref = bulk("Ti", "hcp", a=2.95, c=4.68).repeat([3, 3, 2])
        current = ref.copy()
        del current[0]
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 1
        assert result["interstitial_count"] == 0


class TestEdgeCases:
    """Edge cases and error handling."""
    
    def test_different_atom_counts(self):
        """Reference and current can have different atom counts."""
        ref = bulk("Cu", "fcc", cubic=True).repeat(3)
        current = ref.copy()
        
        # Remove several atoms
        for i in [50, 30, 10, 0]:
            del current[i]
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 4
    
    def test_small_system(self):
        """Test with minimal atom count."""
        ref = bulk("Cu", "fcc", a=3.6)  # Just 1 atom
        current = ref.copy()
        
        result = pyscal3.wigner_seitz_analysis(current, ref)
        
        assert result["vacancy_count"] == 0
        assert result["interstitial_count"] == 0
        assert result["occupancy"][0] == 1
