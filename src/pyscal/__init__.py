"""pyscal — alias package for `pyscal3`.

The library is distributed on PyPI as ``pyscal3`` (the name ``pyscal`` on PyPI
refers to a separate, unrelated package). This module makes ``import pyscal``
and ``from pyscal.<sub> import …`` work as exact synonyms for ``pyscal3``.
"""
import sys as _sys
import pyscal3 as _pyscal3

# Replace this module with pyscal3 in sys.modules so that:
#   import pyscal               -> returns the pyscal3 module
#   from pyscal import find_neighbors
#   from pyscal.structures import make_crystal
# all behave identically to the pyscal3 equivalents.
_sys.modules[__name__] = _pyscal3
