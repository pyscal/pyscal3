"""Tests for Steinhardt bond-order parameters."""
import numpy as np
import pyscal3
from pyscal3.structures import make_crystal


def test_q4_q6_bcc():
    """q4 and q6 for BCC should match known values."""
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
    assert round(np.mean(q4), 2) == 0.51
    assert round(np.mean(q6), 2) == 0.63


def test_q4_q6_bcc_averaged():
    """Averaged q4 and q6 for BCC."""
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
    assert round(np.mean(q4), 2) == 0.51
    assert round(np.mean(q6), 2) == 0.63


def test_q3_bcc_averaged():
    """Averaged q3 for BCC should be zero."""
    atoms = make_crystal("bcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    [q3] = pyscal3.steinhardt_parameter(atoms, l=[3], averaged=True)
    assert round(np.mean(q3), 2) == 0.00


def test_q12_fcc():
    """q12 for FCC (normal and averaged)."""
    atoms = make_crystal("fcc", lattice_constant=1.0, repetitions=(4, 4, 4))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    [q12] = pyscal3.steinhardt_parameter(atoms, l=[12])
    assert round(np.mean(q12), 2) == 0.60

    [q12_avg] = pyscal3.steinhardt_parameter(atoms, l=[12], averaged=True)
    assert round(np.mean(q12_avg), 2) == 0.60


def test_q4_q6_via_ase_bulk():
    """q4/q6 for W BCC via ase.build.bulk."""
    from ase.build import bulk
    atoms = bulk("W", a=1.0, cubic=True).repeat([4, 4, 4])
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0.9)

    q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
    assert round(np.mean(q4), 2) == 0.51
    assert round(np.mean(q6), 2) == 0.63


def test_numpy_integer_l_accepted():
    from pyscal3.structures import make_crystal

    atoms = make_crystal("fcc", lattice_constant=4.0, repetitions=(3, 3, 3))
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    q_np = pyscal3.steinhardt_parameter(atoms, l=np.int64(6))[0]
    q_py = pyscal3.steinhardt_parameter(atoms, l=6)[0]
    assert np.allclose(q_np, q_py)
    q_arr = pyscal3.steinhardt_parameter(atoms, l=np.array([4, 6]))
    assert len(q_arr) == 2
    w = pyscal3.wigner_w_parameter(atoms, l=np.int32(6))[0]
    assert np.allclose(w, -0.01316, atol=1e-4)


def test_qlm_matches_scipy_spherical_harmonics():
    """q_lm of each atom is the mean of Y_lm over its neighbour vectors.

    The reference is scipy.special.sph_harm_y (with the Condon-Shortley
    phase), evaluated at the stored neighbour angles, for every m including
    the stored real and imaginary parts.
    """
    from ase import Atoms
    from scipy.special import sph_harm_y

    rng = np.random.default_rng(11)
    length = (120 * 11.8) ** (1 / 3)
    atoms = Atoms("Cu120", positions=rng.uniform(0, length, (120, 3)),
                  cell=[length] * 3, pbc=True)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.5)
    theta = atoms.info["pyscal_theta"]
    phi = atoms.info["pyscal_phi"]
    for l in range(1, 11):
        (q,) = pyscal3.steinhardt_parameter(atoms, l)
        real = np.asarray(atoms.info.get("pyscal_q%d_real" % l, atoms.arrays.get("pyscal_q%d_real" % l)))
        imag = np.asarray(atoms.info.get("pyscal_q%d_imag" % l, atoms.arrays.get("pyscal_q%d_imag" % l)))
        for i in range(len(atoms)):
            m = np.arange(-l, l + 1)
            ylm = sph_harm_y(l, m[:, None], np.asarray(theta[i])[None, :],
                             np.asarray(phi[i])[None, :]).mean(axis=1)
            np.testing.assert_allclose(real[i], ylm.real, rtol=0, atol=1e-12)
            np.testing.assert_allclose(imag[i], ylm.imag, rtol=0, atol=1e-12)
            expected = np.sqrt(4 * np.pi / (2 * l + 1) * np.sum(np.abs(ylm) ** 2))
            np.testing.assert_allclose(q[i], expected, rtol=1e-12)
