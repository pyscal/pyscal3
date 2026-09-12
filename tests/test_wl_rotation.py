"""W_l must be rotationally invariant (regression: q_lm for m<0 lacked the
(-1)^m phase, which cancels in q_l but not in the 3j contraction)."""
import numpy as np
import pytest

import pyscal3
from pyscal3.structures import make_crystal


def _rotate(atoms, axis, angle):
    axis = np.asarray(axis, float) / np.linalg.norm(axis)
    k = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rot = np.eye(3) + np.sin(angle) * k + (1 - np.cos(angle)) * k @ k
    out = atoms.copy()
    out.set_cell(np.array(atoms.cell) @ rot.T, scale_atoms=False)
    out.positions = atoms.positions @ rot.T
    return out


@pytest.mark.parametrize("structure, lc, what4, what6", [
    ("fcc", 4.05, -0.15932, -0.01316),
    ("bcc", 3.16, 0.15932, 0.01316),
])
@pytest.mark.parametrize("axis, angle", [([1, 2, 3], 0.65), ([0, 1, 1], 1.9)])
def test_what_rotation_invariant(structure, lc, what4, what6, axis, angle):
    atoms = _rotate(make_crystal(structure, lattice_constant=lc, repetitions=(3, 3, 3)), axis, angle)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)
    w4, w6 = pyscal3.wigner_w_parameter(atoms, l=[4, 6])
    assert np.allclose(w4, what4, atol=1e-4)
    assert np.allclose(w6, what6, atol=1e-4)


def test_what6_perturbed_rotation_invariant():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3), noise=0.1)
    rot = _rotate(atoms, [3, -1, 2], 1.1)
    for a in (atoms, rot):
        pyscal3.find_neighbors(a, method="cutoff", cutoff=3.5)
    w_ref = pyscal3.wigner_w_parameter(atoms, l=6, averaged=True)[0]
    w_rot = pyscal3.wigner_w_parameter(rot, l=6, averaged=True)[0]
    assert np.allclose(w_ref, w_rot, atol=1e-10)
