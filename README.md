# pyscal — python Structural Environment Calculator

Complete documentation with examples available [here](https://pyscal.org/).

## What changed in version 4

pyscal grew over the years from a small wrapper around a few C++ routines into a full structural-analysis toolkit. With each addition the original `System` / `Atoms` object model became more rigid: every new descriptor had to be threaded through the same class hierarchy, neighbor lists were tied to a particular `System` instance, and interoperating with the rest of the atomistic-simulation Python ecosystem (ASE, pymatgen, OVITO, LAMMPS dump readers) required converting back and forth between formats.

Version 4 is a ground-up rewrite around two ideas:

- **ASE `Atoms` is the data structure.** pyscal3 no longer ships its own atoms class. Any `ase.Atoms` object — read from a LAMMPS dump, a POSCAR, a CIF, an XYZ, built with `ase.build.bulk`, or constructed by hand — can be passed directly to a pyscal3 function. Results are written back as `atoms.arrays["pyscal_*"]` for per-atom quantities and `atoms.info["pyscal_*"]` for global ones, so they travel with the atoms object and remain accessible to the rest of the ASE ecosystem.

- **Functional API.** Every descriptor is a top-level function in `pyscal3` that takes an `Atoms` object as its first argument. Neighbor lists are computed once with `pyscal3.find_neighbors(atoms, ...)` and reused by every subsequent descriptor. There is no system state to keep in sync, no method-chaining order to remember, and no class to subclass when adding a new descriptor.

The C++ core is unchanged in spirit — voro++ for tessellation, hand-written kernels for Steinhardt-type invariants, pybind11 bindings — but is now exposed through a much thinner Python layer.

A typical session looks like:

```python
import pyscal3
from ase.build import bulk

atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)   # adaptive
q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
labels = pyscal3.identify_ackland_jones(atoms)

print(atoms.arrays["pyscal_q6"].mean())   # ≈ 0.57 for fcc
```

## What is included

This release consolidates a number of new descriptors alongside the existing ones:

- Neighbors — fixed cutoff, adaptive cutoff, SANN, Voronoi
- Steinhardt $q_l$ and averaged $\bar{q}_l$
- Wigner $W_l$ — third-order rotational invariants for distinguishing cubic structures
- Minkowski structure metrics — Voronoi-area-weighted, parameter-free $q_l$
- Solid/liquid classification and clustering
- Disorder parameter
- Angular and $\chi$ parameters; Ackland-Jones structural classifier
- Coordination measures — coordination number, effective coordination, generalised coordination, local density
- Angular distribution function and bond-length distribution function
- Voronoi tessellation — structural vector and Voronoi volume
- Centrosymmetry parameter
- Common neighbor analysis (CNA, adaptive CNA, diamond variants)
- Entropy parameter (Piaggi–Parrinello)
- Warren–Cowley short-range order
- Atomic deformation — strain tensor, von Mises invariant, $D^2_{\min}$, slip vector
- Wigner–Seitz defect analysis — vacancies, interstitials, antisites against a reference
- Atomic Cluster Expansion (ACE) descriptors up to body order four

A complete list with mathematical definitions and example notebooks is available in the [documentation](https://pyscal.org/).

## Installation

pyscal3 can be installed from conda-forge:

```
conda install -c conda-forge pyscal3
```

or from PyPI:

```
pip install pyscal3
```

To build from source:

```
git clone https://github.com/pyscal/pyscal3.git
cd pyscal3
pip install .
```

A C++ compiler with C++17 support is required when building from source.

## Citing the work

If you use pyscal3 in your work, please cite the [following article](https://joss.theoj.org/papers/10.21105/joss.01824):

Sarath Menon, Grisell Díaz Leines and Jutta Rogal (2019). pyscal: A python module for structural analysis of atomic environments. *Journal of Open Source Software* 4(43), 1824. https://doi.org/10.21105/joss.01824

## Works using pyscal

For a list of publications that used pyscal, see [here](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=315020929885190486&as_sdt=5).
