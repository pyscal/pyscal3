from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()

setup(
    name='pyscal3',
    version='3.1.3',
    author='Sarath Menon',
    author_email='sarath.menon@pyscal.org',
    description='Python library written in C++ for calculation of local atomic structural environment',
    long_description=readme,
    long_description_content_type = "text/markdown",
    # tell setuptools to look for any packages under 'src'
    packages=find_packages('src'),
    # tell setuptools that all packages will be under the 'src' directory
    # and nowhere else
    package_dir={'':'src'},
    headers=["src/pyscal3/system.h"],
    ext_modules=[
        Pybind11Extension(
            "pyscal3.csystem",
            ["src/pyscal3/neighbor.cpp", "src/pyscal3/sh.cpp", 
            "src/pyscal3/solids.cpp", "src/pyscal3/voronoi.cpp",
            "src/pyscal3/cna.cpp", "src/pyscal3/centrosymmetry.cpp",
            "src/pyscal3/system_binding.cpp", "lib/voro++/voro++.cc"],
            language='c++',
            include_dirs=['lib/voro++'],
            extra_compile_args=['-O3'],
        ),
    ],
    # add custom build_ext command
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    download_url = 'https://anaconda.org/conda-forge/pyscal',
    url = 'https://pyscal.org',
    install_requires=['pybind11', 'numpy', 'ase', 'pyyaml'],
    classifiers=[
        'Programming Language :: Python :: 3'
    ],
    include_package_data=True,
)
