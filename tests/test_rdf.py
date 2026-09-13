"""Tests for radial distribution function."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_rdf_bcc():
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(10, 10, 10))
    rdf, r = pyscal3.radial_distribution_function(atoms, rmax=2)
    args = np.argsort(rdf)[::-1]
    assert r[args[0]] - 0.86 < 1e-5


def test_rdf_fcc():
    atoms = make_crystal("fcc", lattice_constant=1.0, repetitions=(10, 10, 10))
    rdf, r = pyscal3.radial_distribution_function(atoms, rmax=2)
    args = np.argsort(rdf)[::-1]
    assert r[args[0]] - 0.70 < 1e-5


def _coordination_from_rdf(rdf, r, rho, rcut):
    dr = r[1] - r[0]
    mask = r + dr <= rcut + 1e-9
    shell_vols = 4.0 / 3.0 * np.pi * ((r[mask] + dr) ** 3 - r[mask] ** 3)
    return np.sum(rdf[mask] * rho * shell_vols)


def test_rdf_normalisation_fcc():
    """Integrating rho * g(r) * 4 pi r^2 must count neighbours."""
    from ase.build import bulk

    atoms = bulk("Cu", "fcc", cubic=True).repeat(6)
    rdf, r = pyscal3.radial_distribution_function(atoms, rmin=0, rmax=6.0, bins=120)
    rho = len(atoms) / atoms.get_volume()
    # first shell only (2.55 A): 12 neighbours; up to 6 A: 12+6+24+12+24 = 78
    assert np.isclose(_coordination_from_rdf(rdf, r, rho, 3.0), 12.0, atol=1e-6)
    assert np.isclose(_coordination_from_rdf(rdf, r, rho, 6.0), 78.0, atol=1e-6)


def test_rdf_ideal_gas_tends_to_one():
    from ase import Atoms

    rng = np.random.default_rng(0)
    L = 20.0
    atoms = Atoms("Ar%d" % 4000, positions=rng.uniform(0, L, size=(4000, 3)),
                  cell=[L, L, L], pbc=True)
    rdf, r = pyscal3.radial_distribution_function(atoms, rmin=0, rmax=6.0, bins=30)
    assert abs(rdf[5:].mean() - 1.0) < 0.05
