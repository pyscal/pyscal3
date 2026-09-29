"""Cross-check of the matscipy-neighbours search against the previous C++ search.

Temporary: removed together with pyscal3.neighbors._legacy_find_neighbors.
The expected differences are
- the order of neighbours inside a row for the fixed cutoff and the shell,
- phi of neighbour vectors on the z axis, where it is undefined,
- which of several equidistant candidates are kept at a cut (ties),
- empty rows of small cells: stored as (n, 0) arrays, where the old search
  put them into atoms.info after padding the cell with ghost atoms.
"""
import warnings

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, fcc111

import pyscal3
from pyscal3.neighbors import _legacy_find_neighbors

KEYS = ("neighbors", "neighbordist", "neighborweight", "diff", "r", "theta", "phi",
        "cutoff", "temp_neighbors", "temp_neighbordist")


def _structures():
    rng = np.random.default_rng(1)
    s = {}
    s["fcc_perfect_108"] = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(3)
    s["fcc_perfect_864"] = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6)
    a = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6)
    a.rattle(0.1, seed=2)
    s["fcc_rattled_864"] = a
    a = bulk("Fe", "bcc", a=2.87).repeat((7, 6, 5))
    a.rattle(0.05, seed=3)
    s["bcc_triclinic_210"] = a
    a = bulk("Ti", "hcp", a=2.95, c=4.68).repeat((6, 6, 4))
    a.rattle(0.05, seed=4)
    s["hcp_288"] = a
    a = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4)
    a.rattle(0.05, seed=5)
    a.set_cell(-a.cell, scale_atoms=False)
    a.wrap()
    s["left_handed_256"] = a
    a = fcc111("Al", size=(6, 6, 6), vacuum=8.0, a=4.05)
    a.rattle(0.05, seed=6)
    s["slab_216"] = a
    pos = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4).positions
    s["cluster_256"] = Atoms("Cu256", positions=pos + rng.normal(0, 0.05, pos.shape))
    length = (500 * 11.8) ** (1 / 3)
    s["random_500"] = Atoms("Cu500", positions=rng.uniform(0, length, (500, 3)),
                            cell=[length] * 3, pbc=True)
    s["fcc_prim_1"] = bulk("Cu", "fcc", a=3.61)
    return s


STRUCTURES = _structures()
SETTINGS = {
    "cutoff_3": dict(method="cutoff", cutoff=3.0),
    "cutoff_5": dict(method="cutoff", cutoff=5.0),
    "shell": dict(method="cutoff", cutoff=2.6, shell_thickness=1.0),
    "adaptive": dict(method="cutoff", cutoff="adaptive"),
    "sann": dict(method="cutoff", cutoff="sann"),
    "number_12": dict(method="number", nmax=12),
    "number_8": dict(method="number", nmax=8),
    "number_noassign": dict(method="number", nmax=12, assign_neighbor=False),
}


def _get(atoms, key):
    k = "pyscal_" + key
    if k in atoms.arrays:
        return "arrays", atoms.arrays[k]
    if k in atoms.info:
        return "info", atoms.info[k]
    return None, None


def _rows(atoms):
    _, nb = _get(atoms, "neighbors")
    cols = [_get(atoms, k)[1] for k in ("diff", "neighbordist", "theta", "phi", "r")]
    return [[(int(j), *map(float, v), float(d), float(t), float(p), float(r))
             for j, v, d, t, p, r in zip(nb[a], *(c[a] for c in cols))]
            for a in range(len(atoms))]


def _key(t):
    return (t[0], round(t[1], 6), round(t[2], 6), round(t[3], 6))


@pytest.mark.parametrize("setting", sorted(SETTINGS))
@pytest.mark.parametrize("name", sorted(STRUCTURES))
def test_same_result_as_legacy_search(name, setting):
    kw = SETTINGS[setting]
    old, new = STRUCTURES[name].copy(), STRUCTURES[name].copy()
    errors = []
    for atoms, func in ((old, _legacy_find_neighbors), (new, pyscal3.find_neighbors)):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                func(atoms, **kw)
            errors.append(None)
        except RuntimeError as e:
            errors.append(type(e))
    assert errors[0] == errors[1]
    if errors[0] is not None:
        return

    n = len(old)
    ties = setting == "number_8" and "perfect" in name or name == "fcc_prim_1"
    for key in KEYS:
        where_old, v_old = _get(old, key)
        where_new, v_new = _get(new, key)
        empty = all(len(row) == 0 for row in v_old) if key != "cutoff" else False
        if empty:
            assert where_new == "arrays" and np.shape(v_new) == (n, 0)
            continue
        assert where_old == where_new, key
        if where_old == "arrays" and key != "cutoff":
            assert np.shape(v_old) == np.shape(v_new), key
            assert np.asarray(v_old).dtype == np.asarray(v_new).dtype, key
    np.testing.assert_allclose(np.asarray(_get(new, "cutoff")[1], float),
                               np.asarray(_get(old, "cutoff")[1], float), rtol=1e-13, atol=0)
    for r_old, r_new in zip(_rows(old), _rows(new)):
        if ties:
            np.testing.assert_allclose(sorted(t[4] for t in r_new), sorted(t[4] for t in r_old),
                                       rtol=0, atol=1e-10)
            continue
        r_old, r_new = sorted(r_old, key=_key), sorted(r_new, key=_key)
        assert [_key(t) for t in r_old] == [_key(t) for t in r_new]
        if r_old:
            a, b = np.array(r_old), np.array(r_new)
            np.testing.assert_allclose(b[:, 1:6], a[:, 1:6], rtol=0, atol=1e-10)
            np.testing.assert_allclose(b[:, 7], a[:, 7], rtol=0, atol=1e-10)
            dphi = np.abs(b[:, 6] - a[:, 6])
            dphi = np.minimum(dphi, 2 * np.pi - dphi)
            on_axis = np.hypot(a[:, 1], a[:, 2]) < 1e-9
            assert (dphi[~on_axis] < 1e-10).all()
    _, t_old = _get(old, "temp_neighbors")
    _, t_new = _get(new, "temp_neighbors")
    _, d_old = _get(old, "temp_neighbordist")
    _, d_new = _get(new, "temp_neighbordist")
    for a in range(n):
        assert sorted(zip(np.round(d_old[a], 9), map(int, t_old[a]))) == \
            sorted(zip(np.round(d_new[a], 9), map(int, t_new[a])))
