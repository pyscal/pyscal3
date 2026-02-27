from pyscal3.core import System
import numpy as np


def test_sro():
	sys = System.create.lattice.l12(lattice_constant=4, repetitions=(2,2,2))
	sys.find.neighbors(method='cutoff')
	sys.chemical.short_range_order()
	assert sys.atoms.chemical.short_range_order[0] == 1
	assert np.isclose(sys.atoms.chemical.short_range_order[1], -0.33333333)