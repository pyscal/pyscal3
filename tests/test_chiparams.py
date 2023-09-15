import pytest
import numpy as np
import pyscal3.core as pc

def test_chiparamsbcc():
    sys = pc.System.create.lattice.bcc(repetitions = [3, 3, 3], lattice_constant=4)
    sys.find.neighbors(method='cutoff', cutoff=0)
    sys.calculate.chi_params()
    chip2 = [3, 0, 0, 0, 36, 12, 0, 36, 0]
    assert np.sum(np.array(sys.atoms.angular_parameters.chi_params[2])-np.array(chip2)) == 0


def test_chiparamsfcc():
    sys = pc.System.create.lattice.fcc(repetitions = [5, 5, 5], lattice_constant=4)
    sys.find.neighbors(method='cutoff', cutoff=0)
    sys.calculate.chi_params()
    chip2 =  [6, 0, 0, 0, 24, 12, 0, 24, 0]
    assert np.sum(np.array(sys.atoms.angular_parameters.chi_params[2])-np.array(chip2)) == 0

def test_chiparamsdia():
    sys = pc.System.create.lattice.diamond(repetitions = [3, 3, 3], lattice_constant=4)
    sys.find.neighbors(method='cutoff', cutoff=0)
    sys.calculate.chi_params()
    chip2 =  [0, 0, 0, 0, 6, 0, 0, 0, 0]
    assert np.sum(np.array(sys.atoms.angular_parameters.chi_params[2])-np.array(chip2)) == 0
