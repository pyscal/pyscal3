import pyscal3.core as pc
import os
import numpy as np
from ase.build import bulk

def test_system_init():
	sys = pc.System("tests/files/conf.dump", customkeys=["vx", "vy"])
	assert sys.natoms == 500
	assert sys.atoms.vx[0] == '0.0394436'

	sys.read.file("tests/files/conf.dump.gz")
	assert sys.natoms == 500

	sys.read.file("tests/files/conf.bcc.scaled.dump")
	assert sys.natoms == 2

	sys.read.file("tests/files/POSCAR", format="poscar")
	assert sys.natoms == 42

	sys = pc.System("tests/files/conf.dump", customkeys=["vx", "vy"])
	sys.write.file("test.dump", customkeys=["vx"])

	sys = pc.System("test.dump", customkeys=["vx"])
	assert sys.natoms == 500
	assert sys.atoms.vx[0] == '0.0394436'

	sys = pc.System("tests/files/conf.dump", customkeys=["vx", "vy"])
	aseobj = sys.write.ase(species=["Au"])
	assert np.sum(sys.atoms.positions[0]-aseobj.positions[0]) < 1E-5