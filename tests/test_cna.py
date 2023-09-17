import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc
from ase.build import bulk

#test conventional cna for fcc
def test_cna_cutoff():
    sys = pc.System.create.lattice.fcc(repetitions=(7,7,7), lattice_constant=4.00)
    res = sys.calculate.common_neighbor_analysis(lattice_constant=4.00)
    assert res["fcc"] == 7*7*7*4

def test_cna_a1():
    sys = pc.System.create.lattice.fcc(repetitions=(7,7,7), lattice_constant=4.00)
    res = sys.calculate.common_neighbor_analysis(lattice_constant=None)
    assert res["fcc"] == 7*7*7*4

#now test adaptive
def test_cna_adaptive():
    sys = pc.System.create.lattice.fcc(repetitions=(7,7,7), lattice_constant=4.00)
    sys.calculate.common_neighbor_analysis(lattice_constant=None)
    assert sys.atoms.structure[0] == 1

    sys = pc.System.create.lattice.hcp(repetitions=(7,7,7), lattice_constant=4.00)
    sys.calculate.common_neighbor_analysis(lattice_constant=None)
    assert sys.atoms.structure[0] == 2

    sys = pc.System.create.lattice.bcc(repetitions=(7,7,7), lattice_constant=4.00)
    sys.calculate.common_neighbor_analysis(lattice_constant=None)
    assert sys.atoms.structure[0] == 3


def test_ase_bulks():
    
    al_fcc = bulk("Al")
    fe_bcc = bulk("Fe")
    ti_hcp = bulk("Ti")

    sys = pc.System()
    sys.read.ase(al_fcc)
    cna = sys.calculate.common_neighbor_analysis()
    assert cna["fcc"] == 1

    sys = pc.System()
    sys.read.ase(fe_bcc)
    cna = sys.calculate.common_neighbor_analysis()
    assert cna["bcc"] == 1

    sys = pc.System()
    sys.read.ase(ti_hcp)
    cna = sys.calculate.common_neighbor_analysis()
    assert cna["hcp"] == 2

def test_cna_diamond():
    sys = pc.System.create.lattice.diamond(repetitions=(7,7,7), lattice_constant=4.00)
    sys.calculate.diamond_structure()
    assert sys.atoms.structure[0] == 5