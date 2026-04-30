import pyscal3.core as pc
import os
import numpy as np
from ase.build import bulk

def test_voronoi_props():
    nx = 5
    sys = pc.System.create.lattice.bcc(repetitions = [nx, nx, nx], lattice_constant=3.127)
    sys.find.neighbors(method="voronoi")

    assert sys.atoms.voronoi.vertex.numbers[0][0] == 6
    assert (sys.atoms.voronoi.vertex.positions[0][0][0]+1.5635 < 1E-4)
    assert (sys.atoms.voronoi.volume[0]-15.288104691499992 < 1E-5)
    assert (sys.atoms.voronoi.face.perimeters[0][0]-6.6333687143110005 < 1E-5)
    assert sys.atoms.voronoi.face.vertices[0][0] == 6

def test_voronoi_vector():
    sys = pc.System.create.lattice.fcc(repetitions=(4,4,4))
    sys.find.neighbors(method='voronoi')
    sys.calculate.voronoi_vector()
    assert sys.atoms.voronoi.vector[0][1] == 12

#def test_voronoi_vertices():
#    nx = np.random.randint(1, 10)
#    nverts = (nx**3*12)
#    sys = Structure().lattice.bcc(repetitions = [nx, nx, nx])
#    sys.find_neighbors(method='voronoi', cutoff=0.1)
#    assert len(sys.atoms.voronoi.vertex.unique_positions) == nverts