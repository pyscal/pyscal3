import pytest
import os
import numpy as np
import pyscal3.core as pc


def test_entropy():
	sys = pc.System("tests/files/conf.fcc.Al.dump")
	sys.find.neighbors(method="cutoff", cutoff=0)
	
	lat = (sys.box[0][0]/5)

	assert np.mean(sys.calculate.entropy(1.4*lat, average=True, local=True))+4.1782399110084985 < 0.001