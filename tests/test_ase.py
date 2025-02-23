import pytest
import os,sys,inspect
import numpy as np
from pyscal3.ase import calculate_centrosymmetry, calculate_cna
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