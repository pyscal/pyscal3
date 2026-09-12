"""average_over_neighbors key handling and the entropy `averaged` alias."""
import numpy as np
import pytest

import pyscal3
from pyscal3.structures import make_crystal


def _atoms():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3), noise=0.1)
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=3.5)
    pyscal3.steinhardt_parameter(atoms, l=6)
    return atoms


def test_average_over_neighbors_accepts_prefixed_key():
    atoms = _atoms()
    a = pyscal3.average_over_neighbors(atoms, "q6")
    b = pyscal3.average_over_neighbors(atoms, "pyscal_q6")
    assert np.allclose(a, b)
    # equals the built-in neighbour average (which includes the atom itself)
    avg = pyscal3.steinhardt_parameter(atoms, l=6, averaged=True)[0]
    assert a.shape == avg.shape


def test_average_over_neighbors_custom_array_and_missing_key():
    atoms = _atoms()
    atoms.arrays["myprop"] = np.arange(len(atoms), dtype=float)
    res = pyscal3.average_over_neighbors(atoms, "myprop", include_self=False)
    nb = atoms.arrays.get("pyscal_neighbors")
    if nb is None:
        nb = atoms.info["pyscal_neighbors"]
    expected = np.array([np.mean(atoms.arrays["myprop"][list(row)]) for row in nb])
    assert np.allclose(res, expected)
    with pytest.raises(KeyError):
        pyscal3.average_over_neighbors(atoms, "does_not_exist")


def test_entropy_averaged_alias():
    atoms = make_crystal("fcc", lattice_constant=4.05, repetitions=(3, 3, 3))
    rm = 1.4 * 4.05
    pyscal3.find_neighbors(atoms, method="cutoff", cutoff=rm)
    a = pyscal3.entropy(atoms, rm=rm, average=True)
    b = pyscal3.entropy(atoms, rm=rm, averaged=True)
    assert np.allclose(a, b)
    assert "pyscal_average_entropy" in atoms.arrays
