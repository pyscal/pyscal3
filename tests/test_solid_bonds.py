"""The solid-bond count must follow the threshold and the `right` flag
(regression: a dangling else made every neighbour count as a bond)."""
from pathlib import Path
import numpy as np
from ase.io import read

import pyscal3
from pyscal3.structures import make_crystal

DATA = Path(__file__).resolve().parent / "files"


def _rows(atoms, key):
    if key in atoms.arrays:
        return [np.asarray(r) for r in atoms.arrays[key]]
    return [np.asarray(r) for r in atoms.info[key]]


def test_bond_count_follows_threshold_liquid():
    atoms = read(str(DATA / "conf.lqd.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=False)
    bonds = atoms.arrays["pyscal_bonds"]
    sij = _rows(atoms, "pyscal_sij")
    expected = np.array([np.sum(s > 0.5) for s in sij])
    assert np.array_equal(bonds, expected)
    # in a liquid most bonds are not solid
    nn = np.array([len(s) for s in sij])
    assert bonds.mean() < 0.5 * nn.mean()


def test_bond_count_right_false():
    atoms = read(str(DATA / "conf.lqd.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6,
                        cluster=False, right=False)
    bonds = atoms.arrays["pyscal_bonds"]
    sij = _rows(atoms, "pyscal_sij")
    expected = np.array([np.sum(s < 0.5) for s in sij])
    assert np.array_equal(bonds, expected)


def test_bond_count_perfect_crystal():
    atoms = make_crystal("fcc", lattice_constant=4.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=False)
    assert np.all(atoms.arrays["pyscal_bonds"] == 12)
    assert np.all(atoms.arrays["pyscal_solid"] == 1)


def test_bonds_criterion_actually_filters():
    """Requiring more solid bonds than a liquid atom has must remove it."""
    atoms = read(str(DATA / "conf.lqd.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.0, cluster=False)
    solid6 = atoms.arrays["pyscal_solid"].copy()
    pyscal3.find_solids(atoms, bonds=1, threshold=0.5, avgthreshold=0.0, cluster=False)
    solid1 = atoms.arrays["pyscal_solid"].copy()
    assert solid6.sum() < solid1.sum()
