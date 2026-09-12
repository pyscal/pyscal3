"""pyscal — alias package for `pyscal3`.

The library is distributed on PyPI as ``pyscal3`` (the name ``pyscal`` on PyPI
refers to a separate, unrelated package). This module makes ``import pyscal``
and ``from pyscal.<sub> import …`` work as exact synonyms for ``pyscal3``.
"""
import importlib as _importlib
import sys as _sys

import pyscal3 as _pyscal3

# Register the submodules under the alias name *before* replacing this module,
# so that ``import pyscal.descriptors`` resolves to the already imported
# ``pyscal3.descriptors`` instead of executing the module a second time (which
# would create a duplicate copy and load the C++ extension twice).
_SUBMODULES = (
    "_bridge",
    "csl",
    "csystem",
    "descriptors",
    "neighbors",
    "structures",
    "structures.creator",
    "structures.grain_boundary",
    "trajectory",
)
for _name in _SUBMODULES:
    _sys.modules[__name__ + "." + _name] = _importlib.import_module("pyscal3." + _name)

# Replace this module with pyscal3 in sys.modules so that:
#   import pyscal               -> returns the pyscal3 module
#   from pyscal import find_neighbors
#   from pyscal.structures import make_crystal
# all behave identically to the pyscal3 equivalents.
_sys.modules[__name__] = _pyscal3
