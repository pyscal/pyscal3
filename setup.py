"""
Build configuration for pyscal3 C++ extension.

All project metadata lives in pyproject.toml (PEP 621).
This file only defines the C++ extension module.
"""
import sys

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

# MSVC uses /O2; GCC/Clang use -O3
if sys.platform == "win32":
    extra_compile_args = ["/O2"]
else:
    extra_compile_args = ["-O3"]

setup(
    ext_modules=[
        Pybind11Extension(
            "pyscal3.csystem",
            [
                "src/pyscal3/neighbor.cpp",
                "src/pyscal3/sh.cpp",
                "src/pyscal3/solids.cpp",
                "src/pyscal3/voronoi.cpp",
                "src/pyscal3/cna.cpp",
                "src/pyscal3/centrosymmetry.cpp",
                "src/pyscal3/entropy.cpp",
                "src/pyscal3/puremath.cpp",
                "src/pyscal3/system_binding.cpp",
                "lib/voro++/voro++.cc",
            ],
            language="c++",
            include_dirs=["lib/voro++"],
            extra_compile_args=extra_compile_args,
        ),
    ],
    cmdclass={"build_ext": build_ext},
    headers=["src/pyscal3/system.h"],
)
