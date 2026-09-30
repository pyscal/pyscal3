"""The thread-count setting of pyscal3."""
import os
import subprocess
import sys

import pyscal3


def test_set_and_get_num_threads():
    before = pyscal3.get_num_threads()
    try:
        pyscal3.set_num_threads(3)
        assert pyscal3.get_num_threads() == 3
        pyscal3.set_num_threads(1)
        assert pyscal3.get_num_threads() == 1
        pyscal3.set_num_threads(None)
        assert pyscal3.get_num_threads() >= 1
    finally:
        pyscal3.set_num_threads(before)


def _threads_at_import(env):
    code = "import pyscal3; print(pyscal3.get_num_threads())"
    full = {k: v for k, v in os.environ.items() if k not in ("PYSCAL_NUM_THREADS", "OMP_NUM_THREADS")}
    full.update(env)
    out = subprocess.run([sys.executable, "-c", code], env=full, capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    return int(out.stdout.strip())


def test_default_from_environment():
    assert _threads_at_import({"PYSCAL_NUM_THREADS": "2"}) == 2
    assert _threads_at_import({"OMP_NUM_THREADS": "3"}) == 3
    assert _threads_at_import({"PYSCAL_NUM_THREADS": "2", "OMP_NUM_THREADS": "3"}) == 2
    assert _threads_at_import({}) >= 1


import contextlib
import warnings

import numpy as np
import pytest
from ase.build import bulk


@contextlib.contextmanager
def threads(n):
    before = pyscal3.get_num_threads()
    pyscal3.set_num_threads(n)
    try:
        yield
    finally:
        pyscal3.set_num_threads(before)


def _crystal():
    # large enough that the per-atom and per-pair loops are split into chunks
    atoms = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(10)
    atoms.rattle(0.1, seed=4)
    return atoms


def _stored(atoms):
    out = {}
    for store in (atoms.arrays, atoms.info):
        for key, value in store.items():
            if key.startswith("pyscal_") and not isinstance(value, (bool, str)):
                out[key] = value
    return out


def _assert_same(a, b):
    assert sorted(a) == sorted(b)
    for key in a:
        x, y = a[key], b[key]
        if isinstance(x, list):
            assert x == y, key
        else:
            np.testing.assert_array_equal(np.asarray(x), np.asarray(y), err_msg=key)


@pytest.mark.parametrize("kwargs", [
    dict(method="cutoff", cutoff=3.0),
    dict(method="cutoff", cutoff=4.3),
    dict(method="cutoff", cutoff=2.6, shell_thickness=1.0),
    dict(method="cutoff", cutoff="adaptive"),
    dict(method="cutoff", cutoff="sann"),
    dict(method="number", nmax=12),
])
def test_neighbors_do_not_depend_on_threads(kwargs):
    results = []
    for n in (1, 4):
        atoms = _crystal()
        with threads(n), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pyscal3.find_neighbors(atoms, **kwargs)
        results.append(_stored(atoms))
    _assert_same(*results)


def _descriptors(atoms, reference):
    a = atoms.copy()
    out = {}
    out["q"] = np.concatenate(pyscal3.steinhardt_parameter(a, [4, 6]))
    out["avg_q"] = np.concatenate(pyscal3.steinhardt_parameter(a, [4, 6], averaged=True))
    out["w"] = np.concatenate(pyscal3.wigner_w_parameter(a, [4, 6]))
    out["avg_w"] = np.concatenate(pyscal3.wigner_w_parameter(a, [6], averaged=True))
    out["disorder"] = pyscal3.disorder(a, averaged=True)
    b = a.copy()
    pyscal3.find_solids(b, cluster=True)
    out["sij"] = np.concatenate([np.ravel(r) for r in b.arrays.get("pyscal_sij", b.info.get("pyscal_sij"))])
    out["solid"] = b.arrays["pyscal_solid"]
    out["cluster"] = b.arrays["pyscal_cluster"]
    cp, cos = pyscal3.chi_params(a, angles=True)
    out["chi"] = cp
    out["cosines"] = np.concatenate([np.asarray(c) for c in cos])
    out["angular"] = pyscal3.angular_criteria(a)
    out["entropy"] = pyscal3.entropy(a, rm=3.0, average=True)
    numbers = a.get_atomic_numbers()
    numbers[::3] = 28
    a.set_atomic_numbers(numbers)
    out["sro"] = np.asarray(pyscal3.short_range_order(a, average=False))
    a.arrays["pyscal_prop"] = np.arange(len(a), dtype=float)
    out["avg"] = pyscal3.average_over_neighbors(a, "prop")
    out["strain"] = pyscal3.atomic_strain(a, reference)
    out["d2min"] = pyscal3.d2min(a, reference)
    out["slip"] = pyscal3.slip_vector(a, reference)
    c = atoms.copy()
    out["centro"] = pyscal3.centrosymmetry(c, nmax=12)
    return out


@pytest.mark.parametrize("cutoff", [3.0, 4.3])
def test_descriptors_do_not_depend_on_threads(cutoff):
    results = []
    for n in (1, 4):
        with threads(n):
            atoms, ref = _crystal(), _crystal()
            atoms.positions += np.random.default_rng(1).normal(0, 0.02, atoms.positions.shape)
            pyscal3.find_neighbors(atoms, method="cutoff", cutoff=cutoff)
            pyscal3.find_neighbors(ref, method="cutoff", cutoff=cutoff)
            results.append(_descriptors(atoms, ref))
    for key in results[0]:
        np.testing.assert_array_equal(results[0][key], results[1][key], err_msg=key)
