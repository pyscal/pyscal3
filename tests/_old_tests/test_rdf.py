import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc


def test_rdf_bcc():
    sys = pc.System.create.lattice.bcc(repetitions = [10, 10, 10])
    rdf, r = sys.calculate.radial_distribution_function(rmax=2)

    args = np.argsort(rdf)[::-1]
    assert(r[args[0]]-0.86 < 1E-5)

def test_rdf_fcc():
    sys = pc.System.create.lattice.fcc(repetitions = [10, 10, 10])
    rdf, r = sys.calculate.radial_distribution_function(rmax=2)

    args = np.argsort(rdf)[::-1]
    assert(r[args[0]]-0.70 < 1E-5)
