import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc
import pyscal3.ase as psa
from ase.build import bulk


def test_calculate_centrosymmetry():
    ti_hcp = bulk("Ti", orthorhombic=True)
    q = psa.calculate_centrosymmetry(ti_hcp, nmax=12)
    assert np.round(np.mean(np.array(q)), decimals=2) == 8.7


def test_calculate_cna():
    al_fcc = bulk("Al")
    fe_bcc = bulk("Fe")
    ti_hcp = bulk("Ti")

    cna = psa.calculate_cna(al_fcc)
    assert cna["fcc"] == 1

    cna = psa.calculate_cna(fe_bcc)
    assert cna["bcc"] == 1

    cna = psa.calculate_cna(ti_hcp)
    assert cna["hcp"] == 2


def test_calculate_steinhardt_parameter():
    w_bcc = bulk("W", a=1.0, cubic=True).repeat([4, 4, 4])
    q = psa.calculate_steinhardt_parameter(
        ase_atoms=w_bcc,
        q=[4, 6],
        nmax=12,
        method='cutoff',
        cutoff=0.9,
    )
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51, "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63, "Calculated q4 value is wrong!"

    w_bcc = bulk("W", a=1.0, cubic=True).repeat([4, 4, 4])
    q = psa.calculate_steinhardt_parameter(
        ase_atoms=w_bcc,
        q=[4, 6],
        nmax=12,
        method='cutoff',
        cutoff=0.9,
        averaged=True,
    )
    assert np.round(np.mean(np.array(q[0])), decimals=2) == 0.51, "Calculated q4 value is wrong!"
    assert np.round(np.mean(np.array(q[1])), decimals=2) == 0.63, "Calculated q4 value is wrong!"


def test_calculate_diamond_structure():
    # the ASE diamond structure is wrong so we use pyscal to generate the structure
    c_diamond = pc.System.create.lattice.diamond(
        element="C",
        repetitions=(2,2,2),
        lattice_constant=4.00,
    ).convert_to.ase()
    diamond_dict = psa.calculate_diamond_structure(c_diamond)
    assert diamond_dict['cubic diamond'] == len(c_diamond)


def test_calculate_radial_distribution_function():
    w_bcc = bulk("W", a=1.00).repeat([10, 10, 10])
    al_fcc = bulk("Al", a=1.00).repeat([10, 10, 10])

    rdf, r = psa.calculate_radial_distribution_function(w_bcc, rmax=2)
    args = np.argsort(rdf)[::-1]
    assert(r[args[0]] - 0.86 < 1E-5)

    rdf, r = psa.calculate_radial_distribution_function(al_fcc, rmax=2)
    args = np.argsort(rdf)[::-1]
    assert(r[args[0]] - 0.70 < 1E-5)


def test_get_symmetry():
    w_bcc = bulk("W", cubic=True)
    al_fcc = bulk("Al", cubic=True)
    assert psa.get_symmetry(w_bcc)["international_symbol"] == "Im-3m"
    assert psa.get_symmetry(al_fcc)["international_symbol"] == "Fm-3m"
