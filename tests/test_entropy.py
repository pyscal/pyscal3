"""Tests for entropy parameter."""

from pathlib import Path
import numpy as np
import pytest
import pyscal3
from pyscal3.structures import make_crystal

DATA = Path(__file__).resolve().parent / "files"


def _rows(atoms, key):
    if key in atoms.arrays:
        return [np.asarray(r) for r in atoms.arrays[key]]
    return [np.asarray(r) for r in atoms.info[key]]


def _reference_entropy(atoms, rm, sigma=0.2, rstart=0.001, h=0.001, local=False):
    """Piaggi-Parrinello entropy (without the 2 pi k_B prefactor) per atom."""
    dists = _rows(atoms, "pyscal_neighbordist")
    cutoffs = atoms.arrays["pyscal_cutoff"]
    nsteps = int((rm - rstart) / h)
    r = rstart + h * np.arange(nsteps + 1)
    rho_global = len(atoms) / atoms.get_volume()
    out = np.zeros(len(atoms))
    for i, d in enumerate(dists):
        rho = len(d) / (4.0 / 3.0 * np.pi * cutoffs[i] ** 3) if local else rho_global
        g = np.exp(-(r[:, None] - d[None, :]) ** 2 / (2 * sigma**2)).sum(axis=1)
        g = g / (4 * np.pi * rho * r**2 * np.sqrt(2 * np.pi * sigma**2))
        integrand = np.where(g > 1e-30, (g * np.log(np.where(g > 1e-30, g, 1.0)) - g + 1) * r**2, r**2)
        out[i] = -rho * h * (0.5 * integrand[0] + integrand[1:-1].sum() + 0.5 * integrand[-1])
    return out


@pytest.mark.parametrize("local", [False, True])
def test_entropy_matches_reference(local):
    np.random.seed(1)
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3), noise=0.15)
    rm = 1.4 * 4.05
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=rm)
    ent = pyscal3.entropy(atoms, rm=rm, local=local)
    ref = _reference_entropy(atoms, rm, local=local)
    assert np.allclose(ent, ref, atol=1e-8)
    assert np.std(ent) > 1e-4   # noisy structure: per-atom values differ


def test_entropy_local_uses_each_atoms_density():
    """Regression: local mode used atom 0's density for every atom."""
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    del atoms[0]                      # create a vacancy -> unequal coordination
    rm = 1.4 * 4.05
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=rm)
    ent = pyscal3.entropy(atoms, rm=rm, local=True)
    assert np.allclose(ent, _reference_entropy(atoms, rm, local=True), atol=1e-8)


def test_entropy_solid_vs_liquid():
    from ase.io import read

    solid = read(str(DATA / "conf.fcc.Al.dump"), format="lammps-dump-text")
    lat = np.linalg.norm(solid.cell[0]) / 5
    rm = 1.4 * lat
    pyscal3.find_neighbors(solid, method="cutoff", cutoff=rm)
    ent_solid = pyscal3.entropy(solid, rm=rm, average=True)
    assert np.std(ent_solid) < 0.05

    liquid = read(str(DATA / "conf.lqd.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(liquid, method="cutoff", cutoff=rm)
    ent_liquid = pyscal3.entropy(liquid, rm=rm, average=True)
    assert np.mean(ent_solid) < np.mean(ent_liquid)


def test_entropy_warns_when_rm_exceeds_cutoff():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    with pytest.warns(UserWarning, match="neighbor cutoff"):
        pyscal3.entropy(atoms, rm=5.0)
