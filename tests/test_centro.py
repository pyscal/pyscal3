import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc
import pyscal3.crystal_structures as pcs

def test_cs_ges():
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [7, 7, 7], lattice_constant=4.00)
    sys = pc.System()
    sys.box = boxdims
    sys.atoms = atoms

    q = sys.calculate_centrosymmetry(nmax=8)
    #assert np.round(np.mean(np.array(q)), decimals=2) == 0.00

    #q = sys.calculate_centrosymmetry(nmax=12)
    #assert np.round(np.mean(np.array(q)), decimals=2) > 0.00
