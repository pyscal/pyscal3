import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc

def test_angular():
    sys = pc.System.create.lattice.diamond(repetitions = [4, 4, 4])
    sys.find.neighbors(method = 'cutoff', cutoff=0)
    sys.calculate.angular_criteria()

    assert np.round(np.mean(np.array(sys.atoms.angular_parameters.diamond_angle)), decimals=2) == 0.00
