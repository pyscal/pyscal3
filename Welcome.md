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
6. [Voronoi tessellation](examples/06_voronoi_tessellation.ipynb): Voronoi structure vector and Voronoi volumes.
7. [Disorder parameter](docs/descriptors/disorder.md): structural disorder from Steinhardt parameter correlations.
8. [Angular and $\chi$ parameters](docs/descriptors/angular.md): angular criteria for tetrahedral ordering and $\chi$ parameters.
9. [Centrosymmetry parameter](examples/09_centrosymmetry_parameter.ipynb): detecting defects and broken symmetry in crystals.
10. [Entropy parameter](docs/descriptors/entropy.md): pair-entropy fingerprint for distinguishing solid and liquid.
11. [Short-range order](examples/11_short_range_order.ipynb): Warren-Cowley parameters for alloys.
12. [Solid/liquid clustering](docs/descriptors/solid_liquid.md): identifying solid atoms in a melt and clustering by arbitrary conditions.
13. [Trajectory module](docs/guide/files.md): lazy access to multi-frame LAMMPS dump files.
14. [Wigner $W_l$ parameters](examples/14_wigner_w_parameters.ipynb): third-order bond-orientational invariants.
15. [Minkowski structure metrics](examples/15_minkowski_structure_metrics.ipynb): Voronoi-area-weighted Steinhardt parameters.
16. [Ackland-Jones classification](examples/16_ackland_jones_classification.ipynb): fcc/bcc/hcp/icosahedral labels from angular histograms.
17. [Coordination variants](examples/17_coordination_variants.ipynb): coordination number, effective and generalized coordination, local density.
18. [Angular and bond-length distributions](examples/18_angular_bond_distributions.ipynb): ADF and BLDF as local fingerprints.
19. [Deformation descriptors](examples/21_deformation_descriptors.ipynb): atomic strain, von Mises invariant, $D^2_{\min}$ and slip vector.
20. [Wigner-Seitz defect analysis](examples/22_wigner_seitz_defects.ipynb): vacancies, interstitials and antisites against a reference.
21. [ACE descriptors](examples/28_ace_descriptors.ipynb): Atomic Cluster Expansion descriptors.
