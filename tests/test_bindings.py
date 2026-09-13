"""Every function exposed by the C++ extension is used by the Python layer."""
import pathlib
import re

import pyscal3.csystem as pc

SRC = pathlib.Path(pc.__file__).parent


def test_all_bound_functions_are_used():
    # explicit encoding: the sources are UTF-8, the platform default is not (cp1252 on Windows)
    source = "\n".join(p.read_text(encoding="utf-8") for p in SRC.rglob("*.py"))
    bound = [name for name in dir(pc) if not name.startswith("_")]
    assert bound, "extension exposes no functions"
    unused = [name for name in bound if not re.search(r"\bpc\.%s\(" % re.escape(name), source)]
    assert unused == []
