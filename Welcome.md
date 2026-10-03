# pyscal

pyscal is a Python library for the calculation of local atomic structural
environments, including Steinhardt bond-orientational order parameters, from
atomistic simulation data. Any [ASE](https://wiki.fysik.dtu.dk/ase/) `Atoms`
object can be analysed directly; results are written back to the `Atoms`
object. The core routines are written in C++ and exposed through pybind11.

```python
import pyscal
from ase.build import bulk

atoms = bulk("Cu", "fcc", cubic=True).repeat(4)
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)   # adaptive cutoff
q4, q6 = pyscal.steinhardt_parameter(atoms, l=[4, 6])
```

Complete documentation is available at [pyscal.org](https://pyscal.org).

## Examples

1. [Getting started](docs/tour.md): loading structures with ASE, the pyscal workflow, and where results are stored.
2. [Creating structures](docs/guide/structures.md): built-in crystal types, elements by name, custom lattices and grain boundaries.
3. [Finding neighbors](docs/guide/neighbors.md): fixed, adaptive and SANN cutoffs, Voronoi and number-based neighbor methods.
4. [Steinhardt parameters](docs/descriptors/steinhardt.md): bond-orientational order parameters $q_l$ and their neighbor-averaged variants.
5. [Common neighbor analysis](docs/descriptors/cna.md): adaptive and conventional CNA, diamond structure identification.
6. [Voronoi tessellation](docs/descriptors/voronoi.md): Voronoi structure vector and Voronoi volumes.
7. [Disorder parameter](docs/descriptors/disorder.md): structural disorder from Steinhardt parameter correlations.
8. [Angular and $\chi$ parameters](docs/descriptors/angular.md): angular criteria for tetrahedral ordering and $\chi$ parameters.
9. [Centrosymmetry parameter](docs/descriptors/centrosymmetry.md): detecting defects and broken symmetry in crystals.
10. [Entropy parameter](docs/descriptors/entropy.md): pair-entropy fingerprint for distinguishing solid and liquid.
11. [Short-range order](examples/11_short_range_order.ipynb): Warren-Cowley parameters for alloys.
12. [Solid/liquid clustering](docs/descriptors/solid_liquid.md): identifying solid atoms in a melt and clustering by arbitrary conditions.
13. [Trajectory module](docs/guide/files.md): lazy access to multi-frame LAMMPS dump files.
14. [Wigner $W_l$ parameters](docs/descriptors/wigner_w.md): third-order bond-orientational invariants.
15. [Minkowski structure metrics](docs/descriptors/minkowski.md): Voronoi-area-weighted Steinhardt parameters.
16. [Ackland-Jones classification](docs/descriptors/ackland_jones.md): fcc/bcc/hcp/icosahedral labels from angular histograms.
17. [Coordination variants](docs/descriptors/coordination.md): coordination number, effective and generalized coordination, local density.
18. [Angular and bond-length distributions](examples/18_angular_bond_distributions.ipynb): ADF and BLDF as local fingerprints.
19. [Deformation descriptors](docs/descriptors/deformation.md): atomic strain, von Mises invariant, $D^2_{\min}$ and slip vector.
20. [Wigner-Seitz defect analysis](docs/descriptors/wigner_seitz.md): vacancies, interstitials and antisites against a reference.
21. [ACE descriptors](docs/descriptors/ace.md): Atomic Cluster Expansion descriptors.
