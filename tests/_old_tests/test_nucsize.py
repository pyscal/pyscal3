import pytest
import os
import numpy as np
import pyscal3.core as pc

def test_complex_system():
    sys = pc.System('tests/files/cluster.dump')
    sys.find.neighbors(method='cutoff', cutoff=3.63)
    assert 176 == sys.find.solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)

def test_cluster():
    sys = pc.System('tests/files/cluster.dump')
    sys.find.neighbors(method='cutoff', cutoff=3.63)
    sys.find.solids(cluster=False)
    val = sys.find.clusters(sys.atoms.steinhardt.order.sij.solid, largest = True)
    assert 176 == val

def test_cluster_cutoff():
    sys = pc.System('tests/files/cluster.dump')
    sys.find.neighbors(method='cutoff', cutoff=3.63)
    sys.find.solids(cluster=False)
    val = sys.find.clusters(sys.atoms.steinhardt.order.sij.solid, largest = True, cutoff=3.63)
    assert 176 == val

def test_system_nucsize_fraction():
    #create some atoms
    sys = pc.System.create.lattice.bcc(repetitions = [2,2,2], lattice_constant=3.20)    

    #test that atoms are set properly
    assert sys.natoms == 16

    #now calculate nucsize
    sys.find.neighbors(method='cutoff', cutoff=3.63)
    assert 16 == sys.find.solids(bonds=0.8, threshold=0.5, avgthreshold=0.6, cluster=True)