"""
Tests for ACE (Atomic Cluster Expansion) descriptors.
"""

import pytest
import numpy as np
from ase.build import bulk, molecule
from ase import Atoms
import pyscal3


class TestACEBasic:
    """Basic functionality tests for ACE descriptors."""

    def test_fcc_structure(self):
        """Test ACE on FCC copper."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)

        assert "nu1" in result
        assert "nu2" in result
        assert "full" in result
        assert result["full"].shape[0] == len(atoms)
        assert result["nu1"].shape[0] == len(atoms)
        assert result["nu2"].shape[0] == len(atoms)

    def test_bcc_structure(self):
        """Test ACE on BCC iron."""
        atoms = bulk("Fe", "bcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)

        assert result["full"].shape[0] == len(atoms)
        assert np.all(np.isfinite(result["full"]))

    def test_hcp_structure(self):
        """Test ACE on HCP magnesium."""
        atoms = bulk("Mg", "hcp").repeat((3, 3, 2))
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)

        assert result["full"].shape[0] == len(atoms)
        assert np.all(np.isfinite(result["full"]))


class TestACEParameters:
    """Test different parameter combinations."""

    def test_nmax_variation(self):
        """Test different nmax values."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)

        result2 = pyscal3.ace(atoms, nmax=2, lmax=2, nu_max=2)
        result4 = pyscal3.ace(atoms, nmax=4, lmax=2, nu_max=2)

        # More radial basis = more descriptors
        assert result4["nu1"].shape[1] > result2["nu1"].shape[1]
        assert result4["nu2"].shape[1] > result2["nu2"].shape[1]

    def test_lmax_variation(self):
        """Test different lmax values."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)

        result2 = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)
        result4 = pyscal3.ace(atoms, nmax=3, lmax=4, nu_max=2)

        # More angular momentum = more nu2 descriptors
        assert result4["nu2"].shape[1] > result2["nu2"].shape[1]

    def test_nu_max_levels(self):
        """Test different correlation orders."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)

        result1 = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=1)
        result2 = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)
        result3 = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=3)

        # nu_max=1 should only have nu1
        assert "nu1" in result1
        assert "nu2" not in result1

        # nu_max=2 should have nu1 and nu2
        assert "nu1" in result2
        assert "nu2" in result2
        assert "nu3" not in result2

        # nu_max=3 should have nu1, nu2, nu3
        assert "nu1" in result3
        assert "nu2" in result3
        assert "nu3" in result3

        # More correlation = more total descriptors
        assert result2["full"].shape[1] > result1["full"].shape[1]
        assert result3["full"].shape[1] > result2["full"].shape[1]

    def test_cutoff_parameter(self):
        """Test explicit cutoff specification."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=5.0)

        # Use a different cutoff for ACE
        result = pyscal3.ace(atoms, nmax=3, lmax=2, cutoff=4.0)

        assert "pyscal_ace_params" in atoms.info
        assert atoms.info["pyscal_ace_params"]["cutoff"] == 4.0


class TestACEInvariance:
    """Test invariance properties of ACE descriptors."""

    def test_rotation_invariance_cubic_symmetry(self):
        """Test rotation invariance using a cubic symmetry operation (90° z).

        This keeps the cell orthogonal so pyscal3's neighbor finder works.
        """
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result_orig = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2, normalize=False)

        # 90° rotation about z axis (keeps cubic cell shape)
        rot90 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

        atoms_rot = atoms.copy()
        atoms_rot.positions = atoms.positions @ rot90.T
        atoms_rot.cell = np.array(atoms.cell) @ rot90.T
        atoms_rot.wrap()
        pyscal3.find_neighbors(atoms_rot, method="cutoff", cutoff=4.0)
        result_rot = pyscal3.ace(atoms_rot, nmax=3, lmax=2, nu_max=2, normalize=False)

        np.testing.assert_allclose(
            result_orig["full"].mean(axis=0),
            result_rot["full"].mean(axis=0),
            atol=1e-10,
        )

    def test_rotation_invariance_arbitrary(self):
        """Test rotation invariance with arbitrary angle using cluster in large box.

        Uses an isolated cluster inside a large cubic box to avoid
        non-orthogonal cell issues in the neighbor finder.
        """
        from ase import Atoms

        cutoff = 4.0
        nmax, lmax = 3, 2

        # Build cluster centered in a large cubic box
        fcc_big = bulk("Cu", "fcc", cubic=True).repeat(6)
        center = fcc_big.positions.mean(axis=0)
        dists = np.linalg.norm(fcc_big.positions - center, axis=1)
        keep = dists < 8.0
        cluster_pos = fcc_big.positions[keep] - center

        L = 100.0
        cluster = Atoms(
            "Cu" * int(keep.sum()),
            positions=cluster_pos + L / 2,
            cell=[L, L, L],
            pbc=True,
        )
        pyscal3.find_neighbors(cluster, method="cutoff", cutoff=cutoff)
        desc1 = pyscal3.ace(
            cluster, nmax=nmax, lmax=lmax, nu_max=2, cutoff=cutoff, normalize=False
        )

        # 37° rotation around [1,1,1]
        ax = np.array([1.0, 1.0, 1.0])
        ax /= np.linalg.norm(ax)
        ang = 37 * np.pi / 180
        K = np.array([[0, -ax[2], ax[1]], [ax[2], 0, -ax[0]], [-ax[1], ax[0], 0]])
        rot = np.eye(3) + np.sin(ang) * K + (1 - np.cos(ang)) * K @ K

        cluster_r = Atoms(
            "Cu" * int(keep.sum()),
            positions=cluster_pos @ rot.T + L / 2,
            cell=[L, L, L],
            pbc=True,
        )
        pyscal3.find_neighbors(cluster_r, method="cutoff", cutoff=cutoff)
        desc2 = pyscal3.ace(
            cluster_r, nmax=nmax, lmax=lmax, nu_max=2, cutoff=cutoff, normalize=False
        )

        # Compare only interior atoms (all neighbors present)
        interior = np.linalg.norm(cluster.positions - L / 2, axis=1) < 4.0
        interior_idx = np.where(interior)[0]
        assert len(interior_idx) > 0, "No interior atoms found"

        np.testing.assert_allclose(
            desc1["full"][interior_idx].mean(axis=0),
            desc2["full"][interior_idx].mean(axis=0),
            atol=1e-10,
        )

    def test_translation_invariance(self):
        """Test that ACE descriptors are translation invariant."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result_orig = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2, normalize=False)

        # Translate structure - wrap to keep within cell
        atoms_trans = atoms.copy()
        atoms_trans.positions += np.array([1.5, 2.3, -0.7])
        atoms_trans.wrap()  # Wrap positions into cell
        pyscal3.find_neighbors(atoms_trans, method="cutoff", cutoff=4.0)
        result_trans = pyscal3.ace(
            atoms_trans, nmax=3, lmax=2, nu_max=2, normalize=False
        )

        # Mean descriptors should match
        mean_orig = result_orig["full"].mean(axis=0)
        mean_trans = result_trans["full"].mean(axis=0)
        np.testing.assert_allclose(mean_orig, mean_trans, rtol=1e-10, atol=1e-10)

    def test_equivalent_atoms_same_descriptors(self):
        """In perfect crystal, all atoms should have same descriptors."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)

        # All rows should be nearly identical
        mean_desc = result["full"].mean(axis=0)
        deviations = np.abs(result["full"] - mean_desc)
        assert np.all(deviations < 0.01)


class TestACENormalization:
    """Test normalization options."""

    def test_normalized_unit_norm(self):
        """Test that normalized descriptors have unit norm."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2, normalize=True)

        norms = np.linalg.norm(result["full"], axis=1)
        np.testing.assert_allclose(norms, 1.0, rtol=1e-10)

    def test_unnormalized_varies(self):
        """Test that unnormalized descriptors vary in magnitude."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2, normalize=False)

        norms = np.linalg.norm(result["full"], axis=1)
        # All atoms should have same environment, so norms should be similar
        assert np.std(norms) / np.mean(norms) < 0.01


class TestACEDifferentStructures:
    """Test that ACE distinguishes different structures."""

    def test_fcc_vs_bcc(self):
        """Test that FCC and BCC have different descriptors."""
        atoms_fcc = bulk("Cu", "fcc", cubic=True).repeat(2)
        atoms_bcc = bulk("Fe", "bcc", cubic=True).repeat(2)

        pyscal3.find_neighbors(atoms_fcc, method="cutoff", cutoff=4.0)
        pyscal3.find_neighbors(atoms_bcc, method="cutoff", cutoff=4.0)

        result_fcc = pyscal3.ace(atoms_fcc, nmax=4, lmax=3, nu_max=2, normalize=False)
        result_bcc = pyscal3.ace(atoms_bcc, nmax=4, lmax=3, nu_max=2, normalize=False)

        # Mean descriptors should differ
        mean_fcc = result_fcc["full"].mean(axis=0)
        mean_bcc = result_bcc["full"].mean(axis=0)

        # They should not be identical
        assert not np.allclose(mean_fcc, mean_bcc, rtol=0.1)


class TestACEEdgeCases:
    """Test edge cases and special scenarios."""

    def test_single_atom(self):
        """Test ACE on single atom with no neighbors."""
        atoms = Atoms("Cu", positions=[[0, 0, 0]], cell=[10, 10, 10], pbc=True)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)

        # Should return zeros for single atom with no neighbors
        assert result["full"].shape[0] == 1
        assert np.all(np.isfinite(result["full"]))

    def test_minimal_parameters(self):
        """Test with minimal parameter values."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        result = pyscal3.ace(atoms, nmax=1, lmax=0, nu_max=1)

        # Should still work
        assert result["full"].shape[0] == len(atoms)
        assert result["nu1"].shape[1] == 1  # nmax=1, only l=0

    def test_stores_in_atoms(self):
        """Test that results are stored in atoms object."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)
        pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)

        assert "pyscal_ace" in atoms.arrays
        assert "pyscal_ace_params" in atoms.info
        assert atoms.info["pyscal_ace_params"]["nmax"] == 3
        assert atoms.info["pyscal_ace_params"]["lmax"] == 2
        assert atoms.info["pyscal_ace_params"]["nu_max"] == 2


class TestACEDescriptorCount:
    """Test that descriptor counts are as expected."""

    def test_nu1_count(self):
        """Test nu=1 descriptor count equals nmax."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)

        for nmax in [2, 3, 4, 5]:
            result = pyscal3.ace(atoms, nmax=nmax, lmax=2, nu_max=1)
            assert result["nu1"].shape[1] == nmax

    def test_nu2_count(self):
        """Test nu=2 descriptor count scales with nmax and lmax."""
        atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.0)

        # For nu=2: count = sum over n1<=n2 and l
        # = (nmax * (nmax+1) / 2) * (lmax+1)
        result = pyscal3.ace(atoms, nmax=3, lmax=2, nu_max=2)
        expected_nu2 = (3 * 4 // 2) * 3  # 6 * 3 = 18
        assert result["nu2"].shape[1] == expected_nu2
