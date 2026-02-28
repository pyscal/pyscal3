"""Tests for Common Neighbor Parameter (CNP).

Reference:
  Tsuzuki, Branicio & Rino, Comput. Phys. Commun. 177, 518 (2007).
  https://doi.org/10.1016/j.cpc.2007.05.018
"""
import numpy as np
import pytest
import pyscal3
from pyscal3.structures import make_crystal


@pytest.fixture
def fcc_atoms():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    return atoms


@pytest.fixture
def bcc_atoms():
    atoms = make_crystal("bcc", lattice_constant=3.16, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    return atoms


def test_cnp_perfect_fcc(fcc_atoms):
    """Perfect FCC has CNP = 0 for all atoms."""
    cnp = pyscal3.common_neighbor_parameter(fcc_atoms)
    assert np.allclose(cnp, 0.0, atol=1e-10)


def test_cnp_perfect_bcc(bcc_atoms):
    """Perfect BCC has CNP = 0 for all atoms."""
    cnp = pyscal3.common_neighbor_parameter(bcc_atoms)
    assert np.allclose(cnp, 0.0, atol=1e-10)


def test_cnp_hcp_nonzero():
    """Perfect HCP has CNP > 0 (distinguishes from FCC)."""
    hcp = make_crystal("hcp", lattice_constant=3.20, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(hcp, method="cutoff", cutoff=0)
    cnp = pyscal3.common_neighbor_parameter(hcp)
    # HCP has uniform but non-zero CNP
    assert cnp.mean() > 0.1
    assert np.std(cnp) < 1e-6  # all atoms should have same value


def test_cnp_noisy_nonzero(fcc_atoms):
    """Noisy FCC has CNP > 0 (detects disorder)."""
    noisy = make_crystal("fcc", lattice_constant=4.05, repetitions=(4, 4, 4), noise=0.1)
    pyscal3.find_neighbors(noisy, method="cutoff", cutoff=0)
    cnp = pyscal3.common_neighbor_parameter(noisy)
    assert cnp.mean() > 0.0


def test_cnp_stored_in_arrays(fcc_atoms):
    """CNP results are stored in atoms.arrays with pyscal_ prefix."""
    pyscal3.common_neighbor_parameter(fcc_atoms)
    assert "pyscal_cnp" in fcc_atoms.arrays
    assert fcc_atoms.arrays["pyscal_cnp"].shape == (len(fcc_atoms),)


def test_cnp_returns_array(fcc_atoms):
    """Return type is numpy array with correct length."""
    cnp = pyscal3.common_neighbor_parameter(fcc_atoms)
    assert isinstance(cnp, np.ndarray)
    assert len(cnp) == len(fcc_atoms)


def test_cnp_fcc_vs_hcp_discrimination():
    """CNP discriminates FCC (CNP=0) from HCP (CNP>0)."""
    fcc = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    hcp = make_crystal("hcp", lattice_constant=3.20, repetitions=(3, 3, 3))

    pyscal3.find_neighbors(fcc, method="cutoff", cutoff=0)
    pyscal3.find_neighbors(hcp, method="cutoff", cutoff=0)

    cnp_fcc = pyscal3.common_neighbor_parameter(fcc)
    cnp_hcp = pyscal3.common_neighbor_parameter(hcp)

    assert cnp_fcc.mean() < cnp_hcp.mean()
