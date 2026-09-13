"""Grain boundary construction (csl module)."""
import numpy as np

import pyscal3
from pyscal3 import csl
from pyscal3.structures import make_grain_boundary


def test_sigma5_100_tilt_angle():
    thetas = csl.get_theta_m_n_list([1, 0, 0], 5)
    assert len(thetas) >= 1
    assert np.isclose(np.degrees(thetas[0][0]), 36.8699, atol=1e-3)


def test_make_grain_boundary_fcc_sigma5():
    gb = make_grain_boundary(axis=[1, 0, 0], structure="fcc", lattice_constant=4.05,
                             sigma=5, gb_plane=(0, 1, 3), element="Al")
    assert len(gb) > 0
    assert gb.get_volume() > 0
    assert set(gb.get_chemical_symbols()) == {"Al"}
    # bulk-like density: number of atoms consistent with fcc Al
    assert np.isclose(len(gb) / gb.get_volume(), 4 / 4.05**3, rtol=0.05)
    pyscal3.find_neighbors(gb, method="cutoff", cutoff=0)
    cn = pyscal3.coordination_number(gb)
    assert cn.max() <= 14 and np.median(cn) == 12
