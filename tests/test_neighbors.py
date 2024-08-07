import pyscal3.core as pc
import os
import numpy as np
from ase.build import bulk

def test_system_init():
    sys = pc.System.create.lattice.bcc(repetitions = [10,10,10], lattice_constant=3.127)
    sys.find.neighbors(method="cutoff", cutoff=3.6)
    a1 = np.array(sys.atoms.neighbors.distance)
    a2 = np.array([2.708061437633939,
 3.127,
 3.126999999999999,
 2.7080614376339383,
 3.127,
 3.126999999999999,
 2.7080614376339383,
 2.708061437633937,
 3.127,
 3.126999999999999,
 2.7080614376339383,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339356])
    assert np.sum(a1-a2) < 1E-5
    assert np.sum(sys.atoms.neighbors.distance[0]-a2) < 1E-5

    sys.find.neighbors(method="cutoff", cutoff=3.6)
    a1 = np.array(sys.atoms.neighbors.distance[0])
    a2 = np.array([2.7080614376339356,
 2.708061437633937,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339383,
 2.7080614376339383,
 2.7080614376339383,
 2.708061437633939,
 3.126999999999999,
 3.126999999999999,
 3.126999999999999,
 3.127,
 3.127,
 3.127])
    assert np.sum(a1-a2) < 1E-5

    sys.find.neighbors(method="cutoff", cutoff='sann')
    a1 = np.array(sys.atoms.neighbors.distance[0])
    a2 = np.array([2.7080614376339356,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339383,
 2.7080614376339383,
 2.7080614376339383,
 2.708061437633939,
 3.126999999999999,
 3.126999999999999,
 3.126999999999999,
 3.127,
 3.127,
 3.127,
 4.422245809540667])
    assert np.sum(a1-a2) < 1E-5

    sys.find.neighbors(method="number", nmax=8)
    a1 = np.array(sys.atoms.neighbors.distance[0])
    a2 = np.array([2.7080614376339356,
 2.708061437633937,
 2.708061437633937,
 2.708061437633937,
 2.7080614376339383,
 2.7080614376339383,
 2.7080614376339383,
 2.708061437633939])
    assert np.sum(a1-a2) < 1E-5


def test_neighbor_shell():
   sys = pc.System.create.element.Cu(repetitions=(5,5,5))
   sys.find.neighbors(method='cutoff', cutoff=3, shell_thickness=1, cells=False)
   assert len(sys.atoms.neighbors.distance[0]) == 6
   assert np.abs(sys.atoms.neighbors.distance[0][0]-3.61) <= 1E-3

   sys = pc.System.create.element.Cu(repetitions=(5,5,5))
   sys.find.neighbors(method='cutoff', cutoff=3, shell_thickness=1, cells=True)
   assert len(sys.atoms.neighbors.distance[0]) == 6
   assert np.abs(sys.atoms.neighbors.distance[0][0]-3.61) <= 1E-3