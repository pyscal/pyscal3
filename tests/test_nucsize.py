"""Tests for solid finding and clustering."""
from pathlib import Path
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal

DATA = Path(__file__).resolve().parent / "files"


def test_find_solids_cluster():
    from ase.io import read
    atoms = read(str(DATA / "cluster.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.63)
    val = pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)
    assert val == 176


def test_find_clusters_separate():
    from ase.io import read
    atoms = read(str(DATA / "cluster.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.63)
    pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=False)

    # manually cluster using stored solid mask
    condition = atoms.arrays.get("pyscal_solid", np.zeros(len(atoms))) > 0
    val = pyscal3.find_clusters(atoms, condition=condition, largest=True)
    assert val == 176


def test_find_clusters_with_cutoff():
    from ase.io import read
    atoms = read(str(DATA / "cluster.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.63)
    pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=False)

    condition = atoms.arrays.get("pyscal_solid", np.zeros(len(atoms))) > 0
    val = pyscal3.find_clusters(atoms, condition=condition, largest=True, cutoff=3.63)
    assert val == 176


def test_find_solids_fraction():
    atoms = make_crystal("bcc", lattice_constant=3.20, repetitions=(2, 2, 2))
    assert len(atoms) == 16
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.63)
    val = pyscal3.find_solids(atoms, bonds=0.8, threshold=0.5, avgthreshold=0.6, cluster=True)
    assert val == 16
