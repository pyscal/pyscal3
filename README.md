# pyscal

pyscal computes descriptors of the local atomic structure from atomistic simulation data.
It works on [ASE](https://wiki.fysik.dtu.dk/ase/) `Atoms` objects, so any structure that ASE can read or build can be analysed directly, and the results are stored on the same object.
The calculations run in C++ on all CPU cores.

Documentation with examples: [pyscal.org](https://pyscal.org).

```python
import pyscal
from ase.io import read

atoms = read("dump.lammpstrj", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)        # adaptive cutoff
q4, q6 = pyscal.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
counts = pyscal.common_neighbor_analysis(atoms)               # {'fcc': ..., 'hcp': ..., ...}
atoms.arrays["pyscal_structure"]                              # label of each atom
```

## What pyscal computes

- **Neighbors:** fixed cutoff, adaptive cutoff, SANN, a fixed number of neighbors, Voronoi.
- **Crystal structure of each atom:** common neighbor analysis (adaptive and conventional), diamond structure identification, Ackland–Jones classification, Voronoi structure vectors.
- **Orientational order:** Steinhardt parameters $q_l$ and their averaged form, Wigner $W_l$, Minkowski structure metrics, angular and $\chi$ parameters.
- **Solid and liquid:** solid–liquid classification and clustering, disorder parameters, the entropy fingerprint.
- **Defects and deformation:** centrosymmetry, atomic strain, von Mises strain, $D^2_\mathrm{min}$, slip vector, Wigner–Seitz analysis of vacancies, interstitials and antisites.
- **Chemistry and coordination:** Warren–Cowley short range order, coordination numbers, local density, radial, angular and bond length distributions.
- **Machine learning:** Atomic Cluster Expansion (ACE) descriptors.
- **Structures and trajectories:** builders for crystals and grain boundaries, and a lazy reader for LAMMPS trajectories.

## Installation

pyscal is distributed as `pyscal3`, because the name `pyscal` on PyPI belongs to an unrelated package.
After installation, `import pyscal` and `import pyscal3` both work.

```
pip install pyscal3
```

or

```
conda install -c conda-forge pyscal3
```

Building from source (`pip install .` in a clone) needs a C++17 compiler.

pyscal 4 has a new interface built around ASE. Scripts written for pyscal 3 need small changes, described in [Migrating from pyscal 3](https://pyscal.org/docs/guide/migration.html). To keep the old interface, install `pyscal3<4`.

## Third-party code

pyscal includes two C++ libraries in `lib/`, compiled into its extension module:

- [voro++](https://math.lbl.gov/voro++/) by Chris H. Rycroft, for Voronoi tessellation.
- [matscipy-neighbours](https://github.com/libAtoms/matscipy-neighbours) by the libAtoms developers, for the neighbor search. MIT licence, see `lib/matscipy-neighbours/LICENSE.md` and `VENDORED.md` in the same directory.

## Citing

If you use pyscal in your work, please cite the [following article](https://joss.theoj.org/papers/10.21105/joss.01824):

Sarath Menon, Grisell Díaz Leines and Jutta Rogal (2019). pyscal: A python module for structural analysis of atomic environments. *Journal of Open Source Software* 4(43), 1824. https://doi.org/10.21105/joss.01824

For a list of publications that used pyscal, see [Google Scholar](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=315020929885190486&as_sdt=5).

## Contributing

Bug reports, questions and contributions are welcome. See [Contributing](https://pyscal.org/docs/contributing.html).
