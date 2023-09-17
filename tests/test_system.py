import pyscal3.core as pc
import os
import numpy as np
from ase.build import bulk
from pyscal3.atoms import Atoms


def test_system_init():
	sys = pc.System.create.lattice.bcc(repetitions = [10,10,10], lattice_constant=3.127)
	assert len(sys.atoms["positions"]) == 10*10*10*2
	assert sys.triclinic == 0
	assert np.abs(sys.boxdims[0] - sys.box[0][0]) < 1E-5


def test_system_triclinic():
	struct = bulk('Cu').repeat(10)
	sys = pc.System()
	sys.box = np.array(struct.cell)
	atoms = {}
	atoms["positions"] = struct.positions
	sys.atoms = atoms
	assert sys.triclinic == 1
	tb = (sys.rot== np.array(struct.cell).T)
	assert np.prod(tb) == 1
	tb = (sys.rotinv==np.linalg.inv(np.array(struct.cell).T))
	assert np.prod(tb) == 1

def test_nop():
	sys = pc.System.create.lattice.bcc(repetitions = [2, 2, 2], lattice_constant=3.127)

	assert sys.natoms == 16
	assert len(sys.atoms['positions']) == 128

	for a in sys.iter_atoms():
		assert np.sum(a["positions"]) == 0
		break

	natoms = {'positions':[[0,0,0]]}
	sys.atoms += natoms
	assert sys.natoms == 17

def test_embed():
	cu = bulk('Cu')
	sys = pc.System()
	sys.box = np.array(cu.cell)
	sys.atoms = Atoms({"positions": cu.positions})
	sys.modify.embed_in_cubic_box()
	assert np.abs(sys.box[0][0] - 5.105310960166873) < 1E-5

def test_distance():
	sys = pc.System.create.lattice.bcc(repetitions = [2, 2, 2], lattice_constant=3.127)
	dist = sys.calculate.distance([0.0, 0.0, 0.0], [1.5635, 1.5635, 1.5635])
	assert np.abs(dist - 2.708061437633939) < 1E-5	

def test_composition():
	sys = pc.System.create.lattice.l12(repetitions = [2, 2, 2], lattice_constant=3.127)
	c = sys.concentration
	assert c[1] == 0.25
	assert c[2] == 0.75

	c = sys.composition
	assert c[1] == 0.25
	assert c[2] == 0.75

def test_volume():
	sys = pc.System.create.lattice.fcc(repetitions = [10, 10, 10])
	assert sys.volume == 1000

def test_system_init():
	sys = pc.System.create.lattice.custom([[0, 0, 0], [0.5, 0.5, 0.5]], [1, 2], 
		[[1,0,0], [0,1,0], [0,0,1]])
	assert sys.natoms == 2