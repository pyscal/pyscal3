# pyscal

pyscal computes descriptors of the local atomic structure from atomistic simulation data.
It works on [ASE](https://wiki.fysik.dtu.dk/ase/) `Atoms` objects, so any structure that ASE can read or build can be analysed directly, and the results are stored on the same object.
The calculations run in C++ on all CPU cores.

```python
import pyscal
from ase.io import read

atoms = read("dump.lammpstrj", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)
q4, q6 = pyscal.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
pyscal.common_neighbor_analysis(atoms)
```

pyscal is installed with `pip install pyscal3` or `conda install -c conda-forge pyscal3`, and imported as `pyscal`.

::::{grid} 1 2 2 2
:gutter: 3

:::{grid-item-card} Get started
:link: tour
:link-type: doc
Install pyscal and go through a complete analysis, from reading a file to plotting the results.
:::

:::{grid-item-card} Finding neighbors
:link: guide/neighbors
:link-type: doc
The neighbor methods, how they compare, and which one to use.
:::

:::{grid-item-card} Descriptors
:link: descriptors/index
:link-type: doc
What each descriptor measures, when to use it, and examples on realistic structures.
:::

:::{grid-item-card} API reference
:link: api
:link-type: doc
All functions and their parameters.
:::
::::

## What pyscal computes

- **Crystal structure of each atom:** common neighbor analysis, diamond structure identification, Ackland–Jones classification, Voronoi structure vectors.
- **Orientational order:** Steinhardt parameters $q_l$ and their averaged form, Wigner $W_l$, Minkowski structure metrics, angular and $\chi$ parameters.
- **Solid and liquid:** solid–liquid classification and clustering, disorder parameters, the entropy fingerprint.
- **Defects and deformation:** centrosymmetry, atomic strain, von Mises strain, $D^2_\mathrm{min}$, slip vector, Wigner–Seitz analysis of vacancies, interstitials and antisites.
- **Chemistry and coordination:** Warren–Cowley short range order, coordination numbers, local density, radial, angular and bond length distributions.
- **Machine learning:** Atomic Cluster Expansion (ACE) descriptors.

## Citing pyscal

If you use pyscal in your work, please cite:

> S. Menon, G. Díaz Leines and J. Rogal, pyscal: A python module for structural analysis of atomic environments, *Journal of Open Source Software* **4**, 1824 (2019). [doi:10.21105/joss.01824](https://doi.org/10.21105/joss.01824)

Please also cite the original publication of each method you use. They are listed on the page of each descriptor.
