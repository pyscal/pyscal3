"""voronoi_vector must apply the area cutoff to the face it belongs to
(regression: weights were indexed by the number of accepted faces)."""
import numpy as np

import pyscal3
from pyscal3.structures import make_crystal


def _rows(atoms, key):
    if key in atoms.arrays:
        return [np.asarray(r) for r in atoms.arrays[key]]
    return [np.asarray(r) for r in atoms.info[key]]


def _reference_vorovector(fv, vn, vv, w, edge_cutoff=0.05, area_cutoff=0.01):
    vv = np.asarray(vv, float).reshape(-1, 3)
    out = np.zeros(4, dtype=int)
    st = 1
    for fi, vno in enumerate(fv):
        idx = np.asarray(vn[st:st + vno], dtype=int)
        st += vno + 1
        pts = vv[idx]
        edges = np.linalg.norm(pts - np.roll(pts, 1, axis=0), axis=1)
        if w[fi] > area_cutoff:
            ec = int(np.sum(edges / edges.sum() > edge_cutoff))
            if 3 <= ec <= 6:
                out[ec - 3] += 1
    return out


def test_voronoi_vector_matches_reference_on_noisy_fcc():
    np.random.seed(7)
    atoms = make_crystal("fcc", lattice_constant=4.0, repetitions=(4, 4, 4), noise=0.25)
    pyscal3.find_neighbors(atoms, method="voronoi")
    vv = pyscal3.voronoi_vector(atoms)
    fvs = _rows(atoms, "pyscal_face_vertices")
    vns = _rows(atoms, "pyscal_vertex_numbers")
    vvs = _rows(atoms, "pyscal_vertex_vectors")
    ws = _rows(atoms, "pyscal_neighborweight")
    ref = np.array([_reference_vorovector(fv, vn, v, w) for fv, vn, v, w in zip(fvs, vns, vvs, ws)])
    # noisy structure must actually contain faces below the area cutoff
    assert any(np.any(w <= 0.01) for w in ws)
    assert np.array_equal(vv, ref)


def test_voronoi_vector_area_cutoff_removes_small_faces():
    np.random.seed(7)
    atoms = make_crystal("bcc", lattice_constant=3.0, repetitions=(4, 4, 4), noise=0.15)
    pyscal3.find_neighbors(atoms, method="voronoi")
    strict = pyscal3.voronoi_vector(atoms, area_cutoff=0.05).sum(axis=1)
    loose = pyscal3.voronoi_vector(atoms, area_cutoff=0.0).sum(axis=1)
    assert np.all(strict <= loose)
    assert strict.sum() < loose.sum()
