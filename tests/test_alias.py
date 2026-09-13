"""`import pyscal` must be an exact alias of pyscal3 without re-importing
submodules (regression: submodules were executed twice, replacing
pyscal3.descriptors and loading the C++ extension a second time)."""
import sys

import pyscal3


def test_alias_module_identity():
    import pyscal

    assert pyscal is pyscal3
    assert pyscal.__version__ == pyscal3.__version__


def test_submodules_are_not_duplicated():
    import pyscal3.descriptors
    import pyscal3.csystem
    before = pyscal3.steinhardt_parameter

    import pyscal.descriptors
    import pyscal.csystem
    from pyscal.structures import make_crystal
    from pyscal import find_neighbors

    assert pyscal.descriptors is pyscal3.descriptors
    assert pyscal.csystem is pyscal3.csystem
    assert sys.modules["pyscal.descriptors"] is sys.modules["pyscal3.descriptors"]
    assert pyscal3.descriptors.steinhardt_parameter is before
    assert make_crystal is pyscal3.structures.make_crystal
    assert find_neighbors is pyscal3.find_neighbors
