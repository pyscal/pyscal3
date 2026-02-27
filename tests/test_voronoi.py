"""Tests for Voronoi tessellation."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_voronoi_props_bcc():
    atoms = make_crystal("bcc", lattice_constant=3.127, repetitions=(5, 5, 5))
    pyscal3.find_neighbors(atoms, method="voronoi")

    # Check vertex numbers (stored in arrays since uniform shape)
    vn = atoms.arrays["pyscal_vertex_numbers"]
    assert vn[0][0] == 6

    # Check vertex positions
    vp = atoms.arrays["pyscal_vertex_vectors"]
    assert abs(vp[0][0] + 1.5635) < 1e-4

    # Check voronoi volume
    vol = atoms.arrays["pyscal_voronoi_volume"]
    assert abs(vol[0] - 15.288104691499992) < 1e-5

    # Check face perimeters  
    fp = atoms.arrays["pyscal_face_perimeters"]
    assert abs(fp[0][0] - 6.6333687143110005) < 1e-5

    # Check face vertices
    fv = atoms.arrays["pyscal_face_vertices"]
    assert fv[0][0] == 6


def test_voronoi_vector_fcc():
    atoms = make_crystal("fcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="voronoi")
    vv = pyscal3.voronoi_vector(atoms)
    assert vv[0][1] == 12  # 12 square faces for FCC
