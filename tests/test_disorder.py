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


def _reference_disorder(atoms, l=6):
    """Kawasaki-Onuki disorder from the stored q_lm, normalised overlaps."""
    q = atoms.arrays["pyscal_q%d_real" % l] + 1j * atoms.arrays["pyscal_q%d_imag" % l]
    qn = q / np.linalg.norm(q, axis=1)[:, None]
    nb = atoms.info.get("pyscal_neighbors", atoms.arrays.get("pyscal_neighbors"))
    out = np.zeros(len(atoms))
    for i, row in enumerate(nb):
        row = np.asarray(row)
        s_ij = np.real(np.sum(qn[i][None, :] * np.conj(qn[row]), axis=1))
        out[i] = np.mean(2.0 - 2.0 * s_ij)
    return out


def test_disorder_matches_reference_liquid():
    from ase.io import read

    atoms = read(str(DATA / "conf.lqd.dump"), format="lammps-dump-text")
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    dis = pyscal3.disorder(atoms, q=6)
    assert np.allclose(dis, _reference_disorder(atoms, 6), atol=1e-10)


def test_disorder_perfect_crystal_is_zero():
    from pyscal3.structures import make_crystal

    atoms = make_crystal("fcc", lattice_constant=4.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    dis = pyscal3.disorder(atoms, q=6)
    assert np.allclose(dis, 0.0, atol=1e-10)


def test_disorder_uses_the_current_neighbors():
    # q_lm stored by an earlier steinhardt_parameter call with other
    # neighbors must not be reused
    import os
    from ase.io import read
    path = os.path.join(os.path.dirname(__file__), "..", "examples", "conf.lqd.Al.dump")
    stale = read(path, format="lammps-dump-text")
    pyscal3.find_neighbors(stale, method="cutoff", cutoff=3.5)
    pyscal3.steinhardt_parameter(stale, l=6)
    pyscal3.find_neighbors(stale, method="voronoi")
    fresh = read(path, format="lammps-dump-text")
    pyscal3.find_neighbors(fresh, method="voronoi")
    assert np.allclose(pyscal3.disorder(stale, q=6), pyscal3.disorder(fresh, q=6))
    assert np.allclose(stale.arrays["pyscal_q6"], fresh.arrays["pyscal_q6"])
