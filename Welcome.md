# pyscal

pyscal computes descriptors of the local atomic structure from atomistic simulation data.
It works on [ASE](https://wiki.fysik.dtu.dk/ase/) `Atoms` objects, and stores its results on the same object.
The complete documentation is at [pyscal.org](https://pyscal.org).

The pages below are the documentation as runnable notebooks.
Open one with a right click and *Open With → Jupytext Notebook*, then run its cells.

## Get started

- [A first analysis](docs/tour.md): read a structure, find neighbors, compute descriptors, plot and save the results.

## User guide

- [Finding neighbors](docs/guide/neighbors.md): the neighbor methods, how they compare, and which one to use.
- [Building structures](docs/guide/structures.md): crystals, elements, custom lattices and grain boundaries.
- [Reading files and trajectories](docs/guide/files.md): file formats, periodic boundaries, LAMMPS trajectories, writing results.
- [Large systems](docs/guide/large_systems.md): threads, memory and store_rows.

## Descriptors

- [Common neighbor analysis](docs/descriptors/cna.md): fcc, hcp, bcc, icosahedral and diamond labels.
- [Ackland–Jones classification](docs/descriptors/ackland_jones.md): labels from bond angles.
- [Steinhardt parameters](docs/descriptors/steinhardt.md): $q_l$ and the averaged $\bar{q}_l$.
- [Wigner $W_l$ parameters](docs/descriptors/wigner_w.md): third order invariants.
- [Minkowski structure metrics](docs/descriptors/minkowski.md): $q_l$ weighted by Voronoi face areas.
- [Voronoi tessellation](docs/descriptors/voronoi.md): Voronoi vectors and volumes.
- [Angular criteria and $\chi$ parameters](docs/descriptors/angular.md): tetrahedral order and bond angle histograms.
- [Solid–liquid classification](docs/descriptors/solid_liquid.md): solid atoms and clusters in a liquid.
- [Disorder parameter](docs/descriptors/disorder.md): local disorder from bond correlations.
- [Entropy parameter](docs/descriptors/entropy.md): a pair entropy fingerprint.
- [Centrosymmetry](docs/descriptors/centrosymmetry.md): stacking faults, vacancies and surfaces.
- [Atomic deformation](docs/descriptors/deformation.md): strain, $D^2_\mathrm{min}$ and slip vector.
- [Wigner–Seitz analysis](docs/descriptors/wigner_seitz.md): vacancies, interstitials and antisites.
- [Short range order](docs/descriptors/sro.md): Warren–Cowley parameters.
- [Coordination](docs/descriptors/coordination.md): coordination numbers and local density.
- [Distribution functions](docs/descriptors/distributions.md): radial, angular and bond length distributions.
- [ACE descriptors](docs/descriptors/ace.md): Atomic Cluster Expansion.
