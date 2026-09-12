"""Atoms without a cell or with non-periodic directions (molecules,
clusters, slabs) must be handled without periodic images."""
import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, molecule, fcc111
from ase.neighborlist import neighbor_list

import pyscal3
from pyscal3._bridge import effective_periodic_cell


def _ase_cn(atoms, cutoff):
    return np.bincount(neighbor_list("i", atoms, cutoff), minlength=len(atoms))


def _sets(atoms):
    nb = atoms.arrays.get("pyscal_neighbors")
    if nb is None:
        nb = atoms.info["pyscal_neighbors"]
    return [frozenset(int(j) for j in row) for row in nb]


def _ase_sets(atoms, cutoff):
    i, j = neighbor_list("ij", atoms, cutoff)
    out = [set() for _ in range(len(atoms))]
    for a, b in zip(i, j):
        out[a].add(int(b))
    return [frozenset(s) for s in out]


def test_molecule_without_cell():
    water = molecule("H2O")
    pyscal3.find_neighbors(water, method="cutoff", cutoff=1.2)
    assert _sets(water) == [frozenset({1, 2}), frozenset({0}), frozenset({0})]
    d = pyscal3.get_distance(water, water.positions[0], water.positions[1])
    assert np.isclose(d, np.linalg.norm(water.positions[0] - water.positions[1]))


def test_cluster_pbc_false_has_surface():
    cluster = bulk("Cu", "fcc", cubic=True).repeat(3)
    cluster.set_pbc(False)
    cluster.center(vacuum=0.0)   # cell exactly bounding the atoms
    pyscal3.find_neighbors(cluster, method="cutoff", cutoff=3.0)
    assert _sets(cluster) == _ase_sets(cluster, 3.0)
    cn = pyscal3.coordination_number(cluster)
    assert cn.max() == 12 and cn.min() < 12


def test_cluster_adaptive_and_number():
    cluster = bulk("Cu", "fcc", cubic=True).repeat(3)
    cluster.set_pbc(False)
    pyscal3.find_neighbors(cluster, method="cutoff", cutoff=0)
    assert pyscal3.coordination_number(cluster).max() == 12
    pyscal3.find_neighbors(cluster, method="number", nmax=12)
    assert np.all(pyscal3.coordination_number(cluster) == 12)


@pytest.mark.parametrize("vacuum", [0.5, 10.0])
def test_slab_mixed_pbc(vacuum):
    """fcc(111) slab periodic in-plane only; a thin vacuum must not let
    the surfaces see their periodic images."""
    slab = fcc111("Cu", size=(4, 4, 4), vacuum=vacuum)
    assert list(slab.pbc) == [True, True, False]
    pyscal3.find_neighbors(slab, method="cutoff", cutoff=3.0)
    assert _sets(slab) == _ase_sets(slab, 3.0)
    values, counts = np.unique(pyscal3.coordination_number(slab), return_counts=True)
    assert values.tolist() == [9, 12] and counts.tolist() == [32, 32]


def test_wire_periodic_in_one_direction():
    wire = bulk("Cu", "fcc", cubic=True).repeat((6, 2, 2))
    wire.set_pbc([True, False, False])
    pyscal3.find_neighbors(wire, method="cutoff", cutoff=3.0)
    assert _sets(wire) == _ase_sets(wire, 3.0)


def test_effective_cell_geometry():
    slab = fcc111("Cu", size=(2, 2, 3), vacuum=1.0)
    cell, periodic = effective_periodic_cell(slab, pad=3.0)
    assert periodic.tolist() == [True, True, False]
    assert np.allclose(cell[:2], np.array(slab.cell)[:2])
    normal = np.cross(cell[0], cell[1])
    normal /= np.linalg.norm(normal)
    assert np.allclose(np.abs(cell[2] @ normal), np.linalg.norm(cell[2]))   # orthogonal
    extent = np.ptp(slab.positions @ normal)
    assert np.isclose(np.linalg.norm(cell[2]), extent + 6.0)


def test_density_descriptors_require_periodic_cell():
    cluster = bulk("Cu", "fcc", cubic=True).repeat(2)
    cluster.set_pbc(False)
    pyscal3.find_neighbors(cluster, method="cutoff", cutoff=3.0)
    with pytest.raises(ValueError, match="periodic"):
        pyscal3.radial_distribution_function(cluster, rmax=4.0)
    with pytest.raises(ValueError, match="periodic"):
        pyscal3.entropy(cluster, rm=3.0)
