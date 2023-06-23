import pytest
import os
import numpy as np
import pyscal3.core as pc
from pyscal3.crystal_structures import Structure


def test_q_3():
    sys = Structure().lattice.bcc(repetitions = [4, 4, 4])

    #sys.get_neighbors(method = 'voronoi')
    sys.find_neighbors(method = 'cutoff', cutoff=0.9)

    q = sys.calculate_q(3, averaged=True)
    assert np.round(np.mean(np.array(q)), decimals=2) == 0.00
