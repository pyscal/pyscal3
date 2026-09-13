"""Where pyscal stores its results: ASE writers must keep working and no
stale neighbor data may survive a new neighbor search."""
import io

import numpy as np
import pytest
from ase.build import bulk
from ase.io import write

import pyscal3


def _perturbed():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    del atoms[0]                                   # ragged neighbor lists
    atoms.positions += np.random.default_rng(0).normal(scale=0.05, size=atoms.positions.shape)
    return atoms


@pytest.mark.parametrize("builder", [lambda: bulk("Cu", "fcc", cubic=True).repeat(3), _perturbed])
@pytest.mark.parametrize("fmt", ["extxyz", "vasp", "lammps-data"])
def test_ase_write_after_pyscal(builder, fmt):
    atoms = builder()
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    pyscal3.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
    pyscal3.coordination_number(atoms)
    write(io.StringIO(), atoms, format=fmt)


def test_neighbor_vectors_live_in_info():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    assert "pyscal_diff" not in atoms.arrays
    assert np.asarray(atoms.info["pyscal_diff"]).shape == (len(atoms), 12, 3)
    for key in ("pyscal_neighbors", "pyscal_neighbordist", "pyscal_theta", "pyscal_phi", "pyscal_cutoff"):
        assert key in atoms.arrays
    # descriptors that need the vectors still work
    res = pyscal3.ace(atoms, nmax=2, lmax=2, nu_max=2)
    assert res["full"].shape[0] == len(atoms)


def test_switching_neighbor_method_clears_stale_keys():
    atoms = bulk("Fe", "bcc", cubic=True).repeat(3)
    pyscal3.find_neighbors(atoms, method="voronoi")
    assert "pyscal_voronoi_volume" in atoms.arrays
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    leftovers = [k for k in list(atoms.arrays) + list(atoms.info)
                 if "voro" in k or "face" in k or "vertex" in k]
    assert leftovers == []
    with pytest.raises(ValueError, match="Voronoi"):
        pyscal3.voronoi_vector(atoms)
    assert atoms.info["pyscal_neighbor_method"] == "cutoff"


def test_descriptor_results_survive_new_search():
    atoms = bulk("Cu", "fcc", cubic=True).repeat(3)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    q6 = pyscal3.steinhardt_parameter(atoms, l=6)[0]
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    assert np.allclose(atoms.arrays["pyscal_q6"], q6)
