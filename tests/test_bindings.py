"""Every function exposed by the C++ extension is used by the Python layer."""
import pathlib
import re

import pyscal3.csystem as pc

SRC = pathlib.Path(pyscal3.__file__).parent if False else pathlib.Path(pc.__file__).parent


def test_all_bound_functions_are_used():
    source = "\n".join(p.read_text() for p in SRC.rglob("*.py"))
    bound = [name for name in dir(pc) if not name.startswith("_")]
    assert bound, "extension exposes no functions"
    unused = [name for name in bound if not re.search(r"\bpc\.%s\(" % re.escape(name), source)]
    assert unused == []
