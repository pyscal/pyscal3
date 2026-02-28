"""Tests for polyhedral_template_matching."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3 as pc


# ---------------------------------------------------------------------------
# Perfect crystal tests
# ---------------------------------------------------------------------------

class TestPTMPerfectCrystals:
    """PTM on ideal crystal structures."""

    def test_fcc(self):
        atoms = bulk("Al", cubic=True).repeat(3)
        types = pc.polyhedral_template_matching(atoms)
        assert np.all(types == 1), "FCC atoms should all be type 1"
        assert np.all(atoms.arrays["pyscal_ptm_rmsd"] < 1e-4)

    def test_bcc(self):
        atoms = bulk("Fe", cubic=True).repeat(3)
        types = pc.polyhedral_template_matching(atoms)
        assert np.all(types == 3), "BCC atoms should all be type 3"
        assert np.all(atoms.arrays["pyscal_ptm_rmsd"] < 1e-4)

    def test_hcp(self):
        atoms = bulk("Mg", "hcp", a=3.21, c=5.21).repeat((3, 3, 3))
        types = pc.polyhedral_template_matching(atoms)
        assert np.all(types == 2), "HCP atoms should all be type 2"
        # Non-ideal c/a gives small but non-zero RMSD
        assert np.all(atoms.arrays["pyscal_ptm_rmsd"] < 0.01)

    def test_diamond_cubic(self):
        atoms = bulk("Si", "diamond", cubic=True).repeat(2)
        types = pc.polyhedral_template_matching(atoms, structures="all")
        assert np.all(types == 6), "Diamond cubic should be type 6"

    def test_sc(self):
        """Simple cubic Po should be identified as SC."""
        atoms = bulk("Po", "sc", a=3.35).repeat(3)
        types = pc.polyhedral_template_matching(atoms, structures="all")
        assert np.all(types == 5), "SC atoms should be type 5"


# ---------------------------------------------------------------------------
# Output / storage tests
# ---------------------------------------------------------------------------

class TestPTMOutputs:
    """Check output arrays are stored correctly."""

    def test_arrays_stored(self):
        atoms = bulk("Al", cubic=True).repeat(2)
        pc.polyhedral_template_matching(atoms)
        assert "pyscal_ptm_type" in atoms.arrays
        assert "pyscal_ptm_rmsd" in atoms.arrays
        assert "pyscal_ptm_orientation" in atoms.arrays
        assert "pyscal_ptm_interatomic_distance" in atoms.arrays

    def test_orientation_quaternion_shape(self):
        atoms = bulk("Al", cubic=True).repeat(2)
        pc.polyhedral_template_matching(atoms)
        q = atoms.arrays["pyscal_ptm_orientation"]
        assert q.shape == (len(atoms), 4)

    def test_interatomic_distance_fcc(self):
        """FCC Al nearest-neighbor distance ~ 2.86 A."""
        atoms = bulk("Al", cubic=True).repeat(2)
        pc.polyhedral_template_matching(atoms)
        iad = atoms.arrays["pyscal_ptm_interatomic_distance"]
        assert np.allclose(iad, 2.86, atol=0.05)


# ---------------------------------------------------------------------------
# Parameter / API tests
# ---------------------------------------------------------------------------

class TestPTMParameters:
    """Test parameter handling."""

    def test_structure_string(self):
        """Single string 'fcc' should work."""
        atoms = bulk("Al", cubic=True).repeat(2)
        types = pc.polyhedral_template_matching(atoms, structures="fcc")
        assert np.all(types == 1)

    def test_structure_list(self):
        """List of structure names should work."""
        atoms = bulk("Fe", cubic=True).repeat(2)
        types = pc.polyhedral_template_matching(
            atoms, structures=["fcc", "bcc"]
        )
        assert np.all(types == 3)

    def test_invalid_structure_raises(self):
        atoms = bulk("Al", cubic=True)
        with pytest.raises(ValueError, match="Unknown PTM structure"):
            pc.polyhedral_template_matching(atoms, structures="nonsense")

    def test_rmsd_cutoff_strict(self):
        """With an extremely strict RMSD cutoff, noisy atoms become 'other'."""
        atoms = bulk("Al", cubic=True).repeat(3)
        rng = np.random.default_rng(42)
        atoms.positions += rng.normal(0, 0.1, atoms.positions.shape)
        types = pc.polyhedral_template_matching(
            atoms, rmsd_cutoff=0.001
        )
        # Many should be "other" = 0 due to strict cutoff
        assert np.sum(types == 0) > 0

    def test_ptm_types_dict(self):
        """PTM_TYPES should map integer codes to names."""
        assert pc.PTM_TYPES[0] == "other"
        assert pc.PTM_TYPES[1] == "fcc"
        assert pc.PTM_TYPES[2] == "hcp"
        assert pc.PTM_TYPES[3] == "bcc"


# ---------------------------------------------------------------------------
# Noisy structure test
# ---------------------------------------------------------------------------

class TestPTMNoisy:
    """PTM on thermally perturbed structures."""

    def test_fcc_with_noise(self):
        """PTM should still identify FCC with moderate noise."""
        atoms = bulk("Al", cubic=True).repeat(3)
        rng = np.random.default_rng(42)
        atoms.positions += rng.normal(0, 0.05, atoms.positions.shape)
        types = pc.polyhedral_template_matching(atoms, rmsd_cutoff=0.15)
        frac_fcc = np.mean(types == 1)
        assert frac_fcc > 0.9, f"Only {frac_fcc:.1%} FCC with noise"

    def test_bcc_with_noise(self):
        """PTM should still identify BCC with moderate noise."""
        atoms = bulk("Fe", cubic=True).repeat(3)
        rng = np.random.default_rng(42)
        atoms.positions += rng.normal(0, 0.05, atoms.positions.shape)
        types = pc.polyhedral_template_matching(atoms, rmsd_cutoff=0.15)
        frac_bcc = np.mean(types == 3)
        assert frac_bcc > 0.9, f"Only {frac_bcc:.1%} BCC with noise"
