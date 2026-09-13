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

### Fixes

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
