
# pyscal3 - python Structural Environment Calculator


**pyscal3** is a Python library for computing local atomic structural environments from atomistic simulation data. It provides fast C++-backed calculations of [Steinhardt's bond-orientational order parameters](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784), common neighbor analysis, Voronoi tessellation, and many more descriptors — all through a clean functional API built on [ASE](https://wiki.fysik.dtu.dk/ase/) (Atomic Simulation Environment).

## Key features

- **ASE-first API** — all functions take and return standard `ase.Atoms` objects. Results are stored in `atoms.arrays` and `atoms.info` with the `pyscal_` prefix.
- **Fast C++ core** — neighbor finding, Steinhardt parameters, CNA, and other descriptors run in compiled C++ via pybind11.
- **Comprehensive descriptor suite:**
  - [Steinhardt bond-orientational order parameters](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.28.784) and their [averaged](https://aip.scitation.org/doi/full/10.1063/1.2977970) and [disorder](https://doi.org/10.1063/1.3656762) variants
  - [Voronoi tessellation](http://math.lbl.gov/voro++) with face-area-weighted Steinhardt parameters
  - [Common neighbor analysis](https://iopscience.iop.org/article/10.1088/0965-0393/20/4/045021) (CNA, adaptive CNA, diamond CNA)
  - [Solid/liquid classification](https://link.springer.com/chapter/10.1007/b99429) with clustering
  - [Centrosymmetry parameter](https://doi.org/10.1103/PhysRevB.58.11085)
  - [Angular parameters](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.47.15717) and [Ackland-Jones](https://doi.org/10.1103/PhysRevB.73.054104) chi parameters
  - [Cowley short-range order](https://doi.org/10.1103/PhysRev.120.1648)
  - [Entropy parameter](https://doi.org/10.1063/1.4998408) for structure distinction
  - Radial distribution function
- **Structure creation** — built-in routines for common crystals, elements, general lattices, and grain boundaries.
- **Trajectory support** — memory-efficient analysis of large LAMMPS dump trajectories.

## Quick start

```python
from ase.build import bulk
import pyscal3

# Create a structure
atoms = bulk("Cu", "fcc", cubic=True).repeat(3)

# Find neighbors
pyscal3.find_neighbors(atoms, method="cutoff", cutoff=0)

# Calculate descriptors
q = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
print(atoms.arrays["pyscal_q6"])

# Common neighbor analysis
result = pyscal3.common_neighbor_analysis(atoms)
print(result)  # {'fcc': 108, 'hcp': 0, 'bcc': 0, 'ico': 0, 'others': 0}
```

Results are stored directly on the ASE Atoms object:
- **Per-atom data** → `atoms.arrays["pyscal_<name>"]` (NumPy arrays)
- **System-level or ragged data** → `atoms.info["pyscal_<name>"]`

This makes it easy to combine pyscal3 with ASE's I/O, visualization, and analysis tools.

## Citing

If you use pyscal in your work, please cite:

> Sarath Menon, Grisell Díaz Leines and Jutta Rogal (2019). pyscal: A python module for structural analysis of atomic environments. *Journal of Open Source Software*, 4(43), 1824, https://doi.org/10.21105/joss.01824