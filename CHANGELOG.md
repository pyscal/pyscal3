# Changelog

## 4.0.0 (unreleased)

### A new API

pyscal 4 is a rewrite around two ideas: ASE `Atoms` is the data structure,
and every descriptor is a top-level function. The `System` / `Atoms`
classes of pyscal 3 are gone. Neighbor lists are computed once with
`find_neighbors(atoms, ...)` and reused by all descriptors; per-atom results
are written to `atoms.arrays["pyscal_*"]` (uniform shapes) or
`atoms.info["pyscal_*"]` (ragged data and arrays with more than two
dimensions). `import pyscal` and `import pyscal3` are synonyms.

New descriptors: Wigner `W_l`, Minkowski structure metrics, Ackland-Jones
classification, effective and generalized coordination numbers, local
density, angular and bond-length distribution functions, atomic strain,
von Mises strain, `D^2_min`, slip vector, Wigner-Seitz defect analysis,
and ACE descriptors up to body order four.

### Behaviour changes relative to pyscal 3.x

- `radial_distribution_function` now returns a properly normalised `g(r)`
  (it tends to 1 for an ideal gas and its integral counts neighbors); the
  previous values were scaled by an arbitrary constant.
- `short_range_order` computes the Warren-Cowley parameter
  `alpha_AB = 1 - p_AB / c_B` for atoms of the reference species and
  averages over those atoms; species may be given as symbols or atomic
  numbers and default to the two most abundant species. The old mean over
  all atoms was identically zero.
- SANN neighbor lists now contain the correct `m` nearest atoms (an
  off-by-one skipped the 4th nearest neighbor).
- The disorder parameter compares each atom with its actual neighbors;
  values on disordered structures differ from 3.x.
- `find_solids` honours the `bonds` criterion (a dangling `else` had made it
  inert).
- `W_l` values are now rotation invariant (a missing `(-1)^m` phase made them
  depend on the crystal orientation).
- `entropy(local=True)` uses each atom's own local density, the trapezoidal
  integration includes all grid points, and a warning is issued when `rm`
  exceeds the neighbor cutoff.
- `make_crystal` without `element` assigns the lattice type as atomic
  number (type 1 -> Z=1, ...) instead of Z=0, so sublattices remain
  distinguishable.
- Neighbor vectors (`pyscal_diff`) are stored in `atoms.info`, so
  `ase.io.write` keeps working after a pyscal calculation.

### Neighbor search

- `find_neighbors` uses the matscipy-neighbours library (libAtoms, MIT), whose
  C++ core is included in `lib/matscipy-neighbours`. For 131 072 atoms, a
  cutoff with 12 neighbors per atom takes 0.14 s instead of 3.2 s, the
  adaptive, SANN and number methods are 5 to 6 times faster, and a 5 A
  cutoff, where atoms have different numbers of neighbors, takes 2.3 s
  instead of 7.1 s. The stored keys, formats and values are unchanged, except
  for the two points below.
- Candidates at the same distance (to within 1e-10) are ordered by atom index,
  so `method="number"` picks the same neighbors every time when a shell is
  split, for example `nmax=8` in fcc. Before, the choice was left to the
  sorting routine.
- If no atom has a neighbor, the neighbor keys are stored as (n, 0) arrays in
  `atoms.arrays` for any cell size. Small cells used to put them in
  `atoms.info` as lists.
- The `cells` argument of `find_neighbors` is ignored.
- The bonds are also stored as flat arrays in `atoms.info`
  (`pyscal_bond_offsets`, `pyscal_bond_neighbors`, `pyscal_bond_distance`,
  `pyscal_bond_weight`, `pyscal_bond_vector`, `pyscal_bond_theta`,
  `pyscal_bond_phi`, and `pyscal_candidate_*` for the candidate methods), which
  the descriptors read. `find_neighbors(..., store_rows=False)` stores only
  these and skips the per-atom row keys: for 131 072 atoms with a 5 A cutoff
  the search takes 0.31 s instead of 1.33 s, and the atoms can be written to
  extxyz.

### Descriptor speed

The C++ descriptor routines take the neighbor data as flat arrays instead of
nested lists, and the spherical harmonics of a bond are evaluated for all m at
once. For 131 072 atoms with 12 neighbors each (5 A cutoff: ragged neighbor
lists), on one core:

- `steinhardt_parameter(atoms, [4, 6])`: 0.09 s instead of 1.27 s (0.53 s
  instead of 2.26 s); averaged q_l, `wigner_w_parameter` and `disorder`
  improve similarly.
- `find_solids` with clustering: 0.08 s instead of 0.57 s.
- `chi_params`: 0.03 s instead of 0.88 s; `angular_criteria`: 0.015 s instead
  of 0.61 s; `short_range_order` and `average_over_neighbors`: a few ms instead
  of 0.08 s.
- `atomic_strain`, `von_mises_strain`, `d2min` and `slip_vector` (32 000 atoms):
  6 ms instead of about 0.9 s (0.19 s instead of about 3 s).
- `ace` (4000 atoms, 5 A cutoff): 0.37 s instead of 47 s.
- `common_neighbor_analysis`: 1.84 s instead of 79 s for 1 000 188 atoms (0.45
  s instead of 10.3 s for 256 000), and 10 to 12 times faster than before for
  structures with surfaces or liquid. `diamond_structure` improves similarly.
  CNA no longer builds a padded supercell; the labels are unchanged.

Most results are bitwise the same as before. q_l, the q_lm parts, W_l and
everything derived from them (disorder, `find_solids`) can differ in the last
digit (at most about 1e-15), and the strain family by up to about 1e-12
relative for badly conditioned fits.

### Threads

The neighbor search, CNA and `diamond_structure`, and the per-atom loops of
the C++ descriptors run on all CPUs available to the process, with the GIL
released. `pyscal.set_num_threads(n)` and `pyscal.get_num_threads()` set and
return the number of threads; the default can also be set with the
environment variable `PYSCAL_NUM_THREADS`, or else `OMP_NUM_THREADS`. Results
are bitwise the same for any number of threads. The Voronoi tessellation, the
cluster search of `find_clusters` and `ace` stay serial. The two pair loops of
the included matscipy-neighbours code run on pyscal's threads instead of
OpenMP (see `lib/matscipy-neighbours/VENDORED.md`).

For 1 000 188 fcc atoms (3 A cutoff) on the 14 cores of an Apple M4 Pro,
compared with one thread: `find_neighbors` is 7 times faster,
`steinhardt_parameter`, `wigner_w_parameter`, `chi_params` and `entropy`
7 to 10 times, `atomic_strain` 9 times, `common_neighbor_analysis` and
`diamond_structure` 7 to 8 times, and `find_solids` with clustering 5 times.
On 14 threads, `find_neighbors` takes as long as OVITO's
`CutoffNeighborFinder` for 1 000 188 atoms (0.12 s) and up to 1.3 times
longer between 60 000 and 260 000 atoms; `common_neighbor_analysis` is
faster than OVITO's adaptive CNA at every size tested (0.28 s against
0.33 s for 1 000 188 atoms), and `steinhardt_parameter` including the
neighbor search is 3 times faster than freud (0.18 s against 0.56 s).

### Fixes

- `average_over_neighbors` averages each column of a property with several
  values per atom, such as `pyscal_ace`, and returns an array of the same
  shape. It used to return the mean over all values of each atom.
- `identify_ackland_jones` follows the method of Ackland and Jones (2006), as
  in the original implementation of LAMMPS `compute ackland/atom`: it
  chooses its own neighbors from the six nearest atoms, assigns the
  structure with the smallest deviation from the ideal angle counts, and
  labels atoms that match no structure as other. Before, it used the stored
  neighbor list and a simplified decision tree that labelled any unmatched
  atom with an angle near 139 degrees as hcp, so a liquid came out as 99.6 %
  hcp (now 83 % other). It no longer needs `find_neighbors`, stores integer
  labels in `pyscal_structure` like `common_neighbor_analysis` (names are
  still returned), and stores its eight angle counts in
  `pyscal_ackland_chi` instead of computing `pyscal_chiparams`.
- `effective_coordination_number` iterates the weighted mean bond length to
  self-consistency, as Hoppe defines it and as its docstring said. It used
  to stop after one step, so values were too low when bond lengths differ:
  11.28 instead of 11.63 for perfect bcc with 14 neighbors.
- `von_mises_strain` uses the definition of Shimizu, Ogata and Li, with a
  factor 1/6 on the differences of the normal strains. With the factor 1/2
  used before, the value depended on the orientation of the axes: a pure
  shear gave a value up to 1.7 times larger after a rotation.
- `disorder` computes the $q_{lm}$ from the current neighbors. It used to
  reuse values stored by an earlier call with another neighbor list.
- `make_crystal(noise=...)` adds the random displacements once to every atom.
  Before, they were added again to each replicated copy, so atoms in later
  copies moved up to twice as far. A new `seed` argument makes them
  reproducible.
- `common_neighbor_analysis` and `diamond_structure` no longer label every
  atom as "others" when a single atom, for example an isolated atom next to a
  surface, has too few neighbor candidates.
- `find_clusters` and `find_solids` no longer crash with a segmentation fault
  for clusters of a few hundred thousand atoms (the cluster search was
  recursive and overflowed the stack).
- `ace` returns bitwise the same descriptors on repeated calls; complex
  products in the B basis could round differently from call to call.
- `diamond_structure` no longer labels an atom as hexagonal diamond when
  its second shell has the icosahedral CNA signature; such an atom is now
  treated like any other non-diamond site.
- Cell lists are built on fractional coordinates: triclinic, hexagonal and
  rotated cells (primitive fcc, hcp, ASE `fcc111` slabs, LAMMPS triclinic
  boxes) gave wrong neighbors or hung above 250 atoms.
- Ghost-atom padding is driven by the search radius, so cutoffs above half
  the box width and skewed cells no longer lose neighbors silently.
- Voronoi tessellation of triclinic cells (volumes, faces, Minkowski
  metrics, Voronoi vectors) and Voronoi vertex clean-up in such cells.
- Non-periodic directions and cell-less `Atoms` (molecules, clusters,
  slabs) are handled without spurious periodic images.
- Unwrapped coordinates far outside the primary cell.
- Left-handed cells, clustering after Voronoi neighbor finding, ACE
  `nu=3` coupling with Wigner 3j symbols, Voronoi-vector area cutoff,
  LAMMPS triclinic dumps with scaled coordinates, numpy integer orders,
  file handles in `Trajectory`, and the `pyscal` alias no longer imports
  submodules twice.

### Packaging

- Requires Python >= 3.10, numpy, scipy >= 1.15, ase and pyyaml; wheels for
  CPython 3.10-3.14 on Linux x86_64, macOS (x86_64 and arm64) and Windows.
- License metadata corrected to BSD 3-Clause (matching the LICENSE file).
