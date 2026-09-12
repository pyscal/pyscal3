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
