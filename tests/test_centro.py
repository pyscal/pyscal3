import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc

def test_cs_ges():
    sys = pc.System.create.lattice.bcc(repetitions = [7, 7, 7], lattice_constant=4.00)
    q = sys.calculate.centrosymmetry(nmax=8)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00

    q = sys.calculate.centrosymmetry(nmax=4)
    assert np.round(np.mean(np.array(q)), decimals=2) > 0.00
