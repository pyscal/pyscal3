import pytest
import os,sys,inspect
import numpy as np
from pyscal3.ase import calculate_centrosymmetry, calculate_cna, calculate_steinhardt_parameter
from ase.build import bulk


def test_calculate_centrosymmetry():
    ti_hcp = bulk("Ti", orthorhombic=True)
    q = calculate_centrosymmetry(ti_hcp, nmax=12)
    assert np.round(np.mean(np.array(q)), decimals=2) == 8.7


def test_calculate_cna():
    al_fcc = bulk("Al")
    fe_bcc = bulk("Fe")
    ti_hcp = bulk("Ti")

    cna = calculate_cna(al_fcc)
    assert cna["fcc"] == 1

    cna = calculate_cna(fe_bcc)
    assert cna["bcc"] == 1

    cna = calculate_cna(ti_hcp)
    assert cna["hcp"] == 2


def test_calculate_steinhardt_parameter():
    w_bcc = bulk("W", a=1.0, cubic=True).repeat([4, 4, 4])
    q = calculate_steinhardt_parameter(
        ase_atoms=w_bcc,
        q=[4, 6],
        nmax=12,
        method='cutoff',
        cutoff=0.9,
    )
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51, "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63, "Calculated q4 value is wrong!"

    w_bcc = bulk("W", a=1.0, cubic=True).repeat([4, 4, 4])
    q = calculate_steinhardt_parameter(
        ase_atoms=w_bcc,
        q=[4, 6],
        nmax=12,
        method='cutoff',
        cutoff=0.9,
        averaged=True,
    )
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51, "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63, "Calculated q4 value is wrong!"
