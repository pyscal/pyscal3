"""Tests for disorder parameter."""
from pathlib import Path
import numpy as np
import pyscal3

DATA = Path(__file__).resolve().parent / "files"


def test_ordered_disorder():
    """FCC should have low disorder."""
    from ase.io import read
    atoms = read(str(DATA / "conf.fcc.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)

    # Calculate q6 first (disorder is based on it)
    pyscal3.steinhardt_parameter(atoms, l=6)
    dis = pyscal3.disorder(atoms, q=6)
    assert np.mean(dis) < 0.50

    dis_avg = pyscal3.disorder(atoms, q=6, averaged=True)
    assert np.mean(dis_avg) < 0.50


def test_disordered_disorder():
    """Liquid should have high disorder."""
    from ase.io import read
    atoms = read(str(DATA / "conf.lqd.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    pyscal3.steinhardt_parameter(atoms, l=6)

    dis = pyscal3.disorder(atoms, q=6)
    assert np.mean(dis) > 1.00

    dis_avg = pyscal3.disorder(atoms, q=6, averaged=True)
    assert np.mean(dis_avg) > 1.00
