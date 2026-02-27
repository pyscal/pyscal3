import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc
from ase.build import bulk

def test_cs_ges():
    sys = pc.System.create.lattice.bcc(repetitions = [7, 7, 7], lattice_constant=4.00)
    q = sys.calculate.centrosymmetry(nmax=8)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00

    q = sys.calculate.centrosymmetry(nmax=4)
    assert np.round(np.mean(np.array(q)), decimals=2) > 0.00

def test_hcp():
    ti_hcp = bulk("Ti", orthorhombic=True)
    sys = pc.System(ti_hcp, format='ase')
    q = sys.calculate.centrosymmetry(nmax=12)
    assert np.round(np.mean(np.array(q)), decimals=2) == 8.7


