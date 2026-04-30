"""Tests for Wigner W_l (third-order rotational invariant).

Reference values from:
  Steinhardt, Nelson & Ronchetti, Phys. Rev. B 28, 784 (1983), Table I.
  Lechner & Dellago, J. Chem. Phys. 129, 114707 (2008).
"""

import numpy as np
import pytest
import pyscal3
from pyscal3.structures import make_crystal


# --------------- perfect-crystal fixtures ---------------


@pytest.fixture
def fcc_atoms():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    return atoms


@pytest.fixture
def bcc_atoms():
    atoms = make_crystal("bcc", lattice_constant=3.16, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    return atoms


@pytest.fixture
def hcp_atoms():
    atoms = make_crystal("hcp", lattice_constant=3.20, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    return atoms


# --------------- FCC tests ---------------


def test_what6_fcc(fcc_atoms):
    """W-hat_6 for FCC = -0.01316 (Steinhardt 1983, Table I)."""
    [w6] = pyscal3.wigner_w_parameter(fcc_atoms, l=[6])
    assert np.allclose(w6, -0.01316, atol=1e-4)


def test_what4_fcc(fcc_atoms):
    """W-hat_4 for FCC = -0.15932 (Steinhardt 1983, Table I)."""
    [w4] = pyscal3.wigner_w_parameter(fcc_atoms, l=[4])
    assert np.allclose(w4, -0.15932, atol=1e-4)


# --------------- BCC tests ---------------


def test_what6_bcc(bcc_atoms):
    """W-hat_6 for BCC = +0.01316 (Steinhardt 1983, Table I)."""
    [w6] = pyscal3.wigner_w_parameter(bcc_atoms, l=[6])
    assert np.allclose(w6, 0.01316, atol=1e-4)


def test_what4_bcc(bcc_atoms):
    """W-hat_4 for BCC = +0.15932 (Steinhardt 1983, Table I)."""
    [w4] = pyscal3.wigner_w_parameter(bcc_atoms, l=[4])
    assert np.allclose(w4, 0.15932, atol=1e-4)


# --------------- HCP tests ---------------


def test_what6_hcp(hcp_atoms):
    """W-hat_6 for HCP = -0.01244 (Steinhardt 1983, Table I)."""
    [w6] = pyscal3.wigner_w_parameter(hcp_atoms, l=[6])
    assert np.allclose(w6, -0.01244, atol=1e-4)


# --------------- odd l → zero ---------------


def test_odd_l_zero(fcc_atoms):
    """W_l = 0 for odd l (3j symbol vanishes when 3l is odd)."""
    [w5] = pyscal3.wigner_w_parameter(fcc_atoms, l=[5])
    assert np.allclose(w5, 0.0, atol=1e-15)


# --------------- unnormalised W_l ---------------


def test_raw_w6_fcc(fcc_atoms):
    """Raw (unnormalised) W_6 for FCC should be finite and non-zero."""
    [w6_raw] = pyscal3.wigner_w_parameter(fcc_atoms, l=[6], normalized=False)
    # exact value depends on neighbour-averaged q_lm magnitudes;
    # just verify it is non-zero and sign matches W-hat
    assert w6_raw.mean() < 0  # same sign as W-hat_6 for FCC


# --------------- multi-l call ---------------


def test_multi_l(fcc_atoms):
    """Passing l=[4,6] returns two arrays with correct values."""
    w4, w6 = pyscal3.wigner_w_parameter(fcc_atoms, l=[4, 6])
    assert np.allclose(w4, -0.15932, atol=1e-4)
    assert np.allclose(w6, -0.01316, atol=1e-4)


# --------------- single int l ---------------


def test_single_int_l(bcc_atoms):
    """Passing l as int (not list) should work."""
    [w6] = pyscal3.wigner_w_parameter(bcc_atoms, l=6)
    assert np.allclose(w6, 0.01316, atol=1e-4)


# --------------- averaged (Lechner-Dellago) ---------------


def test_averaged_what6_fcc(fcc_atoms):
    """Averaged W-hat_6 for perfect FCC equals the per-atom value."""
    [w6_avg] = pyscal3.wigner_w_parameter(fcc_atoms, l=[6], averaged=True)
    assert np.allclose(w6_avg, -0.01316, atol=1e-4)


def test_averaged_what6_bcc(bcc_atoms):
    """Averaged W-hat_6 for perfect BCC equals the per-atom value."""
    [w6_avg] = pyscal3.wigner_w_parameter(bcc_atoms, l=[6], averaged=True)
    assert np.allclose(w6_avg, 0.01316, atol=1e-4)


# --------------- FCC vs BCC sign discrimination ---------------


def test_fcc_bcc_sign_discriminate():
    """W-hat_6 has opposite sign for FCC and BCC — useful for classification."""
    fcc = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    bcc = make_crystal("bcc", lattice_constant=3.16, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(fcc, method="cutoff", cutoff=0)
    pyscal3.find_neighbors(bcc, method="cutoff", cutoff=0)

    [w6_fcc] = pyscal3.wigner_w_parameter(fcc, l=[6])
    [w6_bcc] = pyscal3.wigner_w_parameter(bcc, l=[6])

    assert w6_fcc.mean() < 0
    assert w6_bcc.mean() > 0
