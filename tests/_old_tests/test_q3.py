import pytest
import os
import numpy as np
import pyscal3.core as pc


def test_q_3():
    sys = pc.System.create.lattice.bcc(repetitions = [4, 4, 4])

    #sys.get_neighbors(method = 'voronoi')
    sys.find.neighbors(method = 'cutoff', cutoff=0.9)

    q = sys.calculate.steinhardt_parameter(3, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00
