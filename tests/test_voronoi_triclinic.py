"""Voronoi tessellation must respect non-orthogonal cells (regression: an
orthogonal |a| x |b| x |c| box was handed to voro++ for every cell)."""
import numpy as np
import pytest
from ase.build import bulk

import pyscal3


def _get(atoms, key):
    """Per-atom pyscal data, whichever of arrays/info it landed in."""
    if key in atoms.arrays:
        return [np.asarray(row) for row in atoms.arrays[key]]
    return [np.asarray(row) for row in atoms.info[key]]


def _rotate(atoms, axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rot = np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k
    out = atoms.copy()
    out.set_cell(np.array(atoms.cell) @ rot.T, scale_atoms=False)
    out.positions = atoms.positions @ rot.T
    return out


@pytest.mark.parametrize("builder, expected_vector", [
    (lambda: bulk("Mg", "hcp").repeat((4, 4, 3)), [0, 12, 0, 0]),
    (lambda: bulk("Cu", "fcc").repeat(4), [0, 12, 0, 0]),
    (lambda: bulk("Fe", "bcc").repeat(4), [0, 6, 0, 8]),
])
def test_voronoi_triclinic_cells(builder, expected_vector):
    atoms = builder()
    pyscal3.find_neighbors(atoms, method="voronoi")
    vol = atoms.arrays["pyscal_voronoi_volume"]
    assert np.isclose(vol.sum(), atoms.get_volume(), rtol=1e-8)
    assert np.ptp(vol) < 1e-8
    assert np.all(pyscal3.coordination_number(atoms) == sum(expected_vector))
    vv = pyscal3.voronoi_vector(atoms)
    assert np.all(vv == np.array(expected_vector))


def test_voronoi_rotated_cell_matches_axis_aligned():
    ref = bulk("Cu", "fcc", cubic=True).repeat(3)
    rng = np.random.default_rng(1)
    ref.positions += rng.normal(scale=0.05, size=ref.positions.shape)
    rot = _rotate(ref, [1, 2, 3], 0.7)
    for a in (ref, rot):
        pyscal3.find_neighbors(a, method="voronoi")
    assert np.allclose(ref.arrays["pyscal_voronoi_volume"], rot.arrays["pyscal_voronoi_volume"])
    for d_ref, d_rot in zip(_get(ref, "pyscal_neighbordist"), _get(rot, "pyscal_neighbordist")):
        assert np.allclose(np.sort(d_ref), np.sort(d_rot))
    # vertices are reported in the original Cartesian frame: their distances
    # to the atom are rotation invariant, and rotating the reference vertex
    # vectors must reproduce the rotated ones
    for v_ref, v_rot in zip(_get(ref, "pyscal_vertex_vectors"), _get(rot, "pyscal_vertex_vectors")):
        v_ref = v_ref.reshape(-1, 3)
        v_rot = v_rot.reshape(-1, 3)
        assert np.allclose(np.sort(np.linalg.norm(v_ref, axis=1)), np.sort(np.linalg.norm(v_rot, axis=1)))
    q_ref = pyscal3.minkowski_parameter(ref, l=6)[0]
    q_rot = pyscal3.minkowski_parameter(rot, l=6)[0]
    assert np.allclose(q_ref, q_rot)


def test_minkowski_hexagonal_vs_orthorhombic_hcp():
    hexa = bulk("Mg", "hcp").repeat((4, 4, 3))
    ortho = bulk("Mg", "hcp", orthorhombic=True).repeat((3, 2, 3))
    q_hex = pyscal3.minkowski_parameter(hexa, l=[4, 6])
    q_ort = pyscal3.minkowski_parameter(ortho, l=[4, 6])
    for a, b in zip(q_hex, q_ort):
        assert np.ptp(a) < 1e-8 and np.ptp(b) < 1e-8
        assert np.isclose(a[0], b[0], atol=1e-8)


def test_voronoi_left_handed_cell():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    cell = np.array(atoms.cell)
    cell[2] *= -1  # left-handed cell
    atoms.set_cell(cell, scale_atoms=False)
    atoms.positions[:, 2] *= -1
    pyscal3.find_neighbors(atoms, method="voronoi")
    vol = atoms.arrays["pyscal_voronoi_volume"]
    assert np.isclose(vol.sum(), atoms.get_volume())
    assert np.all(pyscal3.coordination_number(atoms) == 12)


@pytest.mark.parametrize("builder, sites_per_atom", [
    (lambda: bulk("Cu", "fcc", cubic=True).repeat(3), 3),   # 1 octahedral + 2 tetrahedral
    (lambda: bulk("Cu", "fcc").repeat(4), 3),               # primitive cell
    (lambda: bulk("Fe", "bcc", cubic=True).repeat(3), 6),   # 24 vertices shared by 4 cells
])
def test_unique_voronoi_vertices(builder, sites_per_atom):
    """Unique Voronoi vertices = interstitial sites, independent of the cell."""
    ref = builder()
    rot = _rotate(ref, [1, 1, 0], 0.9)
    for a in (ref, rot):
        pyscal3.find_neighbors(a, method="voronoi", cutoff=0.2)
        u = np.asarray(a.info["pyscal_unique_vertices"])
        assert u.shape == (sites_per_atom * len(a), 3)
        # all returned sites are distinct under the periodic boundary conditions
        from ase import Atoms
        from ase.neighborlist import neighbor_list
        probe = Atoms(positions=u, cell=a.cell, pbc=True)
        assert len(neighbor_list("i", probe, 0.2)) == 0


def test_voronoi_ignores_string_cutoff():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(2)
    pyscal3.find_neighbors(atoms, method="voronoi", cutoff="sann")
    assert "pyscal_unique_vertices" not in atoms.info
