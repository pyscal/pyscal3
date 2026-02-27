import pytest
import os,sys,inspect
import numpy as np
import pyscal3.core as pc

def test_delete():
    sys = pc.System.create.element.Fe()
    assert sys.natoms == 2

    sys.delete(ids=[1])
    assert sys.natoms == 1

    sys = pc.System.create.element.Fe()
    sys.delete(indices=[0])
    assert sys.natoms == 1

    def condition(atom):
    	if atom.ids[0] == 1:
    		return True

    sys = pc.System.create.element.Fe()
    sys.delete(condition=condition)
    assert sys.natoms == 1

    sys = pc.System.create.element.Fe()
    sys.apply_selection(indices=[0])
    sys.delete(selection=True)
    assert sys.natoms == 1

