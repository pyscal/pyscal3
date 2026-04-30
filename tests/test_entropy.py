"""Tests for entropy parameter."""

from pathlib import Path
import numpy as np
import pyscal3

DATA = Path(__file__).resolve().parent / "files"


def test_entropy_fcc():
    from ase.io import read

    atoms = read(str(DATA / "conf.fcc.Al.dump"), format="lammps-dump-text")
    box = atoms.cell
    lat = np.linalg.norm(box[0]) / 5

    # Use a fixed cutoff to avoid adaptive cutoff variability
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=lat * 0.854)
    ent = pyscal3.entropy(atoms, rm=1.4 * lat, average=True, local=True)
    # Entropy should be strongly negative for perfect FCC
    assert np.mean(ent) < -3.0
    # All individual values should be close together (ordered structure)
    assert np.std(ent) < 0.1
