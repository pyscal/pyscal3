from pyscal3.operations.symmetry import get_symmetry
from pyscal3.core import System

def test_bcc():
	sys = System.create.element.Fe()
	sym = get_symmetry(sys)
	assert sym["international_symbol"] == 'Im-3m'	

def test_fcc():
	sys = System.create.element.Cu()
	sym = get_symmetry(sys)
	assert sym["international_symbol"] == 'Fm-3m'	

def test_hcp():
	sys = System.create.element.Mg()
	sym = get_symmetry(sys)
	assert sym["international_symbol"] == 'P6_3/mmc'	

def test_diamond():
	sys = System.create.element.Si()
	sym = get_symmetry(sys)
	assert sym["international_symbol"] == 'Fd-3m'	