"""The version is defined once, in pyscal3/__init__.py."""
import importlib.metadata
import re

import pyscal3


def test_version_matches_installed_metadata():
    assert re.fullmatch(r"\d+\.\d+\.\d+([a-z]+\d+)?", pyscal3.__version__)
    assert importlib.metadata.version("pyscal3") == pyscal3.__version__
