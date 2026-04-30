"""Tests for deformation descriptors."""

import numpy as np
import pytest
from ase.build import bulk

import pyscal3 as pc


def _make_deformed_pair(structure="Al", strain=0.01, cutoff=4.0):
    """Create reference and deformed atom configurations."""
    ref = bulk(structure, cubic=True).repeat(3)
    pc.find_neighbors(ref, method="cutoff", cutoff=cutoff)

    # Apply uniaxial strain in x
    deformed = ref.copy()
    cell = deformed.get_cell()
    cell[0, 0] *= (1 + strain)
    deformed.set_cell(cell, scale_atoms=True)
    pc.find_neighbors(deformed, method="cutoff", cutoff=cutoff)

    return ref, deformed


# ---------------------------------------------------------------------------
# Atomic strain
# ---------------------------------------------------------------------------

class TestAtomicStrain:

    def test_zero_strain_reference(self):
        """Identical configs should give zero strain."""
        ref = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(ref, method="cutoff", cutoff=4.0)

        strain = pc.atomic_strain(ref, ref)
        # Strain should be ~0 for identical configs
        assert np.allclose(strain, 0, atol=1e-10)

    def test_uniaxial_strain_direction(self):
        """Uniaxial strain should appear predominantly in xx component."""
        ref, deformed = _make_deformed_pair("Al", strain=0.05)
        E = pc.atomic_strain(deformed, ref)

        # All Exx should be approximately the same and ~ half the engineering strain
        # (Green-Lagrange: E_xx ≈ ε + ε²/2 ≈ ε for small strain)
        mean_exx = np.nanmean(E[:, 0, 0])
        assert mean_exx > 0.04, f"Exx = {mean_exx} too small"

        # Other diagonal components should be smaller (Poisson contraction),
        # but not necessarily zero
        mean_eyy = np.nanmean(E[:, 1, 1])
        mean_ezz = np.nanmean(E[:, 2, 2])
        assert abs(mean_eyy) < abs(mean_exx), "Eyy should be smaller than Exx"

    def test_shape(self):
        """Strain tensor should be (N, 3, 3)."""
        ref, deformed = _make_deformed_pair("Al", strain=0.01)
        E = pc.atomic_strain(deformed, ref)
        assert E.shape == (len(deformed), 3, 3)

    def test_stored_in_arrays(self):
        ref, deformed = _make_deformed_pair("Al", strain=0.01)
        pc.atomic_strain(deformed, ref)
        assert "pyscal_strain" in deformed.arrays


# ---------------------------------------------------------------------------
# Von Mises strain
# ---------------------------------------------------------------------------

class TestVonMisesStrain:

    def test_zero_for_identical(self):
        """Von Mises should be zero for identical configs."""
        ref = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(ref, method="cutoff", cutoff=4.0)
        vm = pc.von_mises_strain(ref, ref)
        assert np.allclose(vm, 0, atol=1e-10)

    def test_positive_for_shear(self):
        """Shear deformation should give positive von Mises strain."""
        ref = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(ref, method="cutoff", cutoff=4.0)

        # Apply simple shear
        deformed = ref.copy()
        pos = deformed.get_positions()
        pos[:, 0] += 0.05 * pos[:, 1]  # shear in xy
        deformed.set_positions(pos)
        pc.find_neighbors(deformed, method="cutoff", cutoff=4.0)

        vm = pc.von_mises_strain(deformed, ref)
        assert np.all(np.isnan(vm) | (vm > 0))

    def test_shape(self):
        ref, deformed = _make_deformed_pair("Al", strain=0.01)
        vm = pc.von_mises_strain(deformed, ref)
        assert vm.shape == (len(deformed),)


# ---------------------------------------------------------------------------
# D²min
# ---------------------------------------------------------------------------

class TestD2min:

    def test_zero_for_affine(self):
        """D²min should be near-zero for purely affine deformation."""
        ref, deformed = _make_deformed_pair("Al", strain=0.02)
        d2 = pc.d2min(deformed, ref)
        # Perfectly affine: D²min should be small but may not be exactly
        # zero due to boundary effects
        assert np.nanmean(d2) < 0.01, f"D²min = {np.nanmean(d2)} too large"

    def test_positive_for_nonaffine(self):
        """Random displacement should give positive D²min."""
        ref = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(ref, method="cutoff", cutoff=4.0)

        deformed = ref.copy()
        rng = np.random.default_rng(42)
        deformed.positions += rng.normal(0, 0.1, deformed.positions.shape)
        pc.find_neighbors(deformed, method="cutoff", cutoff=4.0)

        d2 = pc.d2min(deformed, ref)
        assert np.nanmean(d2) > 0

    def test_stored_in_arrays(self):
        ref, deformed = _make_deformed_pair("Al", strain=0.01)
        pc.d2min(deformed, ref)
        assert "pyscal_d2min" in deformed.arrays


# ---------------------------------------------------------------------------
# Slip vector
# ---------------------------------------------------------------------------

class TestSlipVector:

    def test_zero_for_identical(self):
        """Slip vector should be zero for identical configs."""
        ref = bulk("Al", cubic=True).repeat(3)
        pc.find_neighbors(ref, method="cutoff", cutoff=4.0)
        sv = pc.slip_vector(ref, ref)
        assert np.allclose(sv, 0, atol=1e-10)

    def test_shape(self):
        ref, deformed = _make_deformed_pair("Al", strain=0.01)
        sv = pc.slip_vector(deformed, ref)
        assert sv.shape == (len(deformed), 3)

    def test_stored_in_arrays(self):
        ref, deformed = _make_deformed_pair("Al", strain=0.01)
        pc.slip_vector(deformed, ref)
        assert "pyscal_slip_vector" in deformed.arrays


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------

class TestDeformationIntegration:

    def test_all_functions_on_same_pair(self):
        """All deformation descriptors should work together."""
        ref, deformed = _make_deformed_pair("Fe", strain=0.03, cutoff=3.5)

        E = pc.atomic_strain(deformed, ref)
        vm = pc.von_mises_strain(deformed, ref)
        d2 = pc.d2min(deformed, ref)
        sv = pc.slip_vector(deformed, ref)

        assert E.shape[0] == len(deformed)
        assert vm.shape[0] == len(deformed)
        assert d2.shape[0] == len(deformed)
        assert sv.shape[0] == len(deformed)
