"""The neighbour search against ASE and against the definitions of each method.

Pair lists come from ase.neighborlist.neighbor_list. The adaptive, SANN and
number methods are rebuilt here from those candidates with the rules the
C++ code implements, so these tests pin the behaviour independently of the
search code.
"""
import warnings

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, fcc111
from ase.neighborlist import neighbor_list

import pyscal3
from pyscal3._bridge import effective_periodic_cell, guess_cutoff


def _structures():
    rng = np.random.default_rng(7)
    s = {}
    a = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4)
    a.rattle(0.1, seed=1)
    s["fcc_rattled"] = a
    a = bulk("Fe", "bcc", a=2.87).repeat((5, 4, 4))
    a.rattle(0.05, seed=2)
    s["bcc_triclinic"] = a
    a = bulk("Ti", "hcp", a=2.95, c=4.68).repeat((4, 4, 3))
    a.rattle(0.05, seed=3)
    s["hcp"] = a
    a = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(3)
    a.rattle(0.05, seed=4)
    a.set_cell(-a.cell, scale_atoms=False)
    a.wrap()
    s["left_handed"] = a
    a = fcc111("Al", size=(4, 4, 5), vacuum=8.0, a=4.05)
    a.rattle(0.05, seed=5)
    s["slab"] = a
    pos = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(3).positions
    # far from the origin, on the negative side
    s["cluster"] = Atoms("Cu108", positions=pos + rng.normal(0, 0.05, pos.shape) - 50.0)
    length = (200 * 11.8) ** (1 / 3)
    s["random"] = Atoms("Cu200", positions=rng.uniform(0, length, (200, 3)),
                        cell=[length] * 3, pbc=True)
    a = bulk("Cu", "fcc", a=3.61)
    a.rattle(0.02, seed=6)
    s["one_atom"] = a
    return s


STRUCTURES = _structures()


def _rows(atoms, key):
    k = "pyscal_" + key
    return atoms.arrays[k] if k in atoms.arrays else atoms.info[k]


def _found(atoms):
    """Per atom: sorted list of (j, dx, dy, dz, d, theta, phi, r) from pyscal."""
    cols = [_rows(atoms, k) for k in ("neighbors", "diff", "neighbordist", "theta", "phi", "r")]
    out = []
    for nb, df, dd, th, ph, rr in zip(*cols):
        out.append(sorted((int(j), *map(float, v), float(d), float(t), float(p), float(r))
                          for j, v, d, t, p, r in zip(nb, df, dd, th, ph, rr)))
    return out


def _ase_pairs(atoms, rmax):
    """i, j, d, v with v = r_i - r_j, from ASE."""
    i, j, D = neighbor_list("ijD", atoms, rmax * (1 + 1e-9) + 1e-12)
    v = -D
    return i, j, np.sqrt((v ** 2).sum(axis=1)), v


def _expected(atoms, i, j, d, v):
    out = [[] for _ in range(len(atoms))]
    for a, b, dd, vv in zip(i, j, d, v):
        r = float(np.sqrt((vv ** 2).sum()))
        out[a].append((int(b), *map(float, vv), float(dd), float(np.arccos(vv[2] / r)),
                       float(np.arctan2(vv[1], vv[0])), r))
    return [sorted(row) for row in out]


def _assert_same(found, expected):
    assert len(found) == len(expected)
    for f, e in zip(found, expected):
        key = lambda t: (t[0], round(t[1], 6), round(t[2], 6), round(t[3], 6))
        f, e = sorted(f, key=key), sorted(e, key=key)
        assert [key(t) for t in f] == [key(t) for t in e]
        if f:
            f, e = np.array(f), np.array(e)
            np.testing.assert_allclose(f[:, 1:6], e[:, 1:6], rtol=0, atol=1e-10)
            np.testing.assert_allclose(f[:, 7], e[:, 7], rtol=0, atol=1e-10)
            dphi = np.abs(f[:, 6] - e[:, 6])
            dphi = np.minimum(dphi, 2 * np.pi - dphi)
            on_axis = np.hypot(e[:, 1], e[:, 2]) < 1e-9     # phi undefined on the z axis
            assert (dphi[~on_axis] < 1e-10).all()


@pytest.mark.parametrize("name", sorted(STRUCTURES))
@pytest.mark.parametrize("rc", [3.0, 5.0])
def test_cutoff_matches_ase(name, rc):
    atoms = STRUCTURES[name].copy()
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=rc)
    i, j, d, v = _ase_pairs(atoms, rc)
    keep = d < rc
    _assert_same(_found(atoms), _expected(atoms, i[keep], j[keep], d[keep], v[keep]))
    has = np.bincount(i[keep], minlength=len(atoms)) > 0
    np.testing.assert_array_equal(atoms.arrays["pyscal_cutoff"], np.where(has, rc, 0.0))


@pytest.mark.parametrize("name", sorted(STRUCTURES))
def test_shell_matches_ase(name):
    atoms = STRUCTURES[name].copy()
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.6, shell_thickness=1.0)
    i, j, d, v = _ase_pairs(atoms, 3.6)
    keep = (d >= 2.6) & (d <= 3.6)
    _assert_same(_found(atoms), _expected(atoms, i[keep], j[keep], d[keep], v[keep]))


def _candidates(atoms, prefactor):
    """Candidates within prefactor * (V / N)^(1/3) of the working cell, per atom
    sorted by distance and then by neighbour index."""
    pad = guess_cutoff(atoms, prefactor)
    cell, _ = effective_periodic_cell(atoms, pad if pad > 0 else 10.0)
    guess = prefactor * np.cbrt(abs(np.linalg.det(cell)) / len(atoms))
    i, j, d, v = _ase_pairs(atoms, guess)
    keep = d <= guess
    i, j, d, v = i[keep], j[keep], d[keep], v[keep]
    key = np.rint(d / 1e-10)      # distances equal to 1e-10 count as ties
    rows = []
    for a in range(len(atoms)):
        sel = np.nonzero(i == a)[0]
        order = sel[np.lexsort((j[sel], key[sel]))]
        rows.append(order)
    return guess, j, d, v, rows


def _select(atoms, j, d, v, rows, counts):
    ii = np.concatenate([np.full(c, a) for a, c in enumerate(counts)]).astype(int)
    idx = np.concatenate([r[:c] for r, c in zip(rows, counts)]).astype(int)
    return _expected(atoms, ii, j[idx], d[idx], v[idx])


def _tie_at_boundary(d, rows, counts):
    """True if some atom has equal distances on both sides of its cut."""
    return any(0 < c < len(r) and abs(d[r[c]] - d[r[c - 1]]) < 1e-9
               for r, c in zip(rows, counts))


def _assert_selection(atoms, j, d, v, rows, counts):
    """Compare with the selected candidates. Where equal distances straddle
    the cut, either choice is valid, so only the distances are compared."""
    expected = _select(atoms, j, d, v, rows, counts)
    if not _tie_at_boundary(d, rows, counts):
        _assert_same(_found(atoms), expected)
        return
    for f, e in zip(_found(atoms), expected):
        np.testing.assert_allclose(sorted(t[4] for t in f), sorted(t[4] for t in e),
                                   rtol=0, atol=1e-10)


def _check_temp(atoms, j, d, rows):
    tj, td = _rows(atoms, "temp_neighbors"), _rows(atoms, "temp_neighbordist")
    for a, r in enumerate(rows):
        np.testing.assert_array_equal(np.asarray(tj[a], dtype=int), j[r])
        np.testing.assert_allclose(np.asarray(td[a], dtype=float), d[r], rtol=0, atol=1e-10)


@pytest.mark.parametrize("name", sorted(STRUCTURES))
def test_adaptive_follows_its_definition(name):
    atoms = STRUCTURES[name].copy()
    guess, j, d, v, rows = _candidates(atoms, 2)
    if min(len(r) for r in rows) < 6:
        with pytest.raises(RuntimeError):
            pyscal3.find_neighbors(atoms, method="cutoff", cutoff="adaptive")
        return
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff="adaptive")
    dcut = np.array([1.2 * (1.0 / 6) * sum(d[r[:6]].tolist()) for r in rows])
    kept = [r[d[r] < dc] for r, dc in zip(rows, dcut)]
    _assert_selection(atoms, j, d, v, kept, [len(r) for r in kept])
    np.testing.assert_allclose(atoms.arrays["pyscal_cutoff"], dcut, rtol=1e-14)
    _check_temp(atoms, j, d, rows)


def _sann_counts(d, rows):
    counts, cuts = [], []
    for r in rows:
        dist = d[r].tolist()
        m, summ = 3, dist[0] + dist[1] + dist[2]
        dcut = summ / (m - 2)
        while m < len(dist) and dcut >= dist[m]:
            summ += dist[m]
            m += 1
            dcut = summ / (m - 2)
        if m == len(dist):
            return None, None
        counts.append(m)
        cuts.append(dcut)
    return counts, np.array(cuts)


@pytest.mark.parametrize("name", sorted(STRUCTURES))
def test_sann_follows_its_definition(name):
    atoms = STRUCTURES[name].copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff="sann")
    # the threshold that converged: the first multiple of 2 that works
    for k in range(1, 10):
        guess, j, d, v, rows = _candidates(atoms, 2 * k)
        if min(len(r) for r in rows) < 3:
            continue
        counts, cuts = _sann_counts(d, rows)
        if counts is not None:
            break
    _assert_selection(atoms, j, d, v, rows, counts)
    np.testing.assert_allclose(atoms.arrays["pyscal_cutoff"], cuts, rtol=1e-14)
    _check_temp(atoms, j, d, rows)


@pytest.mark.parametrize("name", sorted(STRUCTURES))
@pytest.mark.parametrize("nmax", [8, 12])
def test_number_follows_its_definition(name, nmax):
    atoms = STRUCTURES[name].copy()
    guess, j, d, v, rows = _candidates(atoms, 2)
    if min(len(r) for r in rows) < nmax:
        with pytest.raises(RuntimeError):
            pyscal3.find_neighbors(atoms, method="number", nmax=nmax)
        return
    pyscal3.find_neighbors(atoms, method="number", nmax=nmax)
    _assert_selection(atoms, j, d, v, rows, [nmax] * len(atoms))
    np.testing.assert_allclose(atoms.arrays["pyscal_cutoff"], guess, rtol=1e-14)
    _check_temp(atoms, j, d, rows)


def test_number_ties_are_broken_by_index():
    """8 of 12 equidistant neighbours: the ones with the smallest index."""
    atoms = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4)
    pyscal3.find_neighbors(atoms, method="number", nmax=8)
    nb = atoms.arrays["pyscal_neighbors"]
    i, j, d = neighbor_list("ijd", atoms, 2.6)
    for a in range(len(atoms)):
        np.testing.assert_array_equal(nb[a], np.sort(j[i == a])[:8])


def test_storage_uniform_rows():
    atoms = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    for key in ("neighbors", "neighbordist", "neighborweight", "r", "theta", "phi"):
        assert atoms.arrays["pyscal_" + key].shape == (256, 12)
    assert atoms.arrays["pyscal_neighbors"].dtype == np.int64
    assert atoms.arrays["pyscal_neighbordist"].dtype == np.float64
    assert "pyscal_diff" not in atoms.arrays
    assert atoms.info["pyscal_diff"].shape == (256, 12, 3)
    assert atoms.arrays["pyscal_cutoff"].shape == (256,)
    assert atoms.arrays["pyscal_temp_neighbors"].shape == (256, 0)


def test_storage_ragged_rows():
    atoms = STRUCTURES["fcc_rattled"].copy()
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.3)
    counts = [len(row) for row in atoms.info["pyscal_neighbors"]]
    assert len(set(counts)) > 1
    for key in ("neighbors", "neighbordist", "neighborweight", "r", "theta", "phi", "diff"):
        rows = atoms.info["pyscal_" + key]
        assert isinstance(rows, list) and isinstance(rows[0], list)
        assert [len(row) for row in rows] == counts
        assert "pyscal_" + key not in atoms.arrays
    assert isinstance(atoms.info["pyscal_neighbors"][0][0], int)
    assert len(atoms.info["pyscal_diff"][0][0]) == 3


@pytest.mark.parametrize("repeat", [3, 6])
def test_storage_empty_rows_do_not_depend_on_cell_size(repeat):
    """No neighbours at all: (n, 0) arrays, for small and large cells alike."""
    atoms = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(repeat)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.6, shell_thickness=0.5)
    for key in ("neighbors", "neighbordist", "diff", "temp_neighbors"):
        assert atoms.arrays["pyscal_" + key].shape == (len(atoms), 0)


def test_storage_can_be_written_by_ase(tmp_path):
    atoms = STRUCTURES["fcc_rattled"].copy()
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0)
    atoms.write(tmp_path / "out.extxyz")


def test_cutoff_is_strict_and_shell_is_inclusive():
    atoms = Atoms("Cu2", positions=[[0, 0, 0], [2.5, 0, 0]], cell=[20, 20, 20], pbc=False)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.5)
    assert atoms.arrays["pyscal_neighbors"].shape == (2, 0)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=2.0, shell_thickness=0.5)
    assert atoms.arrays["pyscal_neighbors"].tolist() == [[1], [0]]


def test_cells_argument_is_ignored():
    ref = STRUCTURES["fcc_rattled"].copy()
    pyscal3.find_neighbors(ref, method="cutoff", cutoff=3.0)
    for cells in (True, False):
        atoms = STRUCTURES["fcc_rattled"].copy()
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.0, cells=cells)
        np.testing.assert_array_equal(atoms.arrays["pyscal_neighbors"],
                                      ref.arrays["pyscal_neighbors"])


@pytest.mark.parametrize("enabled", [True, False])
def test_garbage_collector_state_is_restored(enabled):
    """Ragged rows are built with the collector paused, then it is restored."""
    import gc

    atoms = STRUCTURES["fcc_rattled"].copy()
    was = gc.isenabled()
    try:
        gc.enable() if enabled else gc.disable()
        pyscal3.find_neighbors(atoms, method="cutoff", cutoff=4.3)
        assert isinstance(atoms.info["pyscal_neighbors"], list)
        assert gc.isenabled() == enabled
        pyscal3.common_neighbor_analysis(atoms)
        assert gc.isenabled() == enabled
    finally:
        gc.enable() if was else gc.disable()
