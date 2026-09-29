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
                # neighbour search, see lib/matscipy-neighbours/VENDORED.md
                "lib/matscipy-neighbours/error.cc",
                "lib/matscipy-neighbours/tools.cc",
                "lib/matscipy-neighbours/memory_space.cc",
                "lib/matscipy-neighbours/cell_list.cc",
                "lib/matscipy-neighbours/neighbour_list.cc",
                "lib/matscipy-neighbours/first_neighbours.cc",
                "lib/matscipy-neighbours/triplet_list.cc",
            ],
            language="c++",
            cxx_std=17,
            include_dirs=["lib/voro++", "lib/matscipy-neighbours"],
            extra_compile_args=extra_compile_args,
        ),
    ],
    cmdclass={"build_ext": build_ext},
)
