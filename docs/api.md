# API reference

All functions are available from `import pyscal` (or, equivalently, `import pyscal3`).
They take an ASE `Atoms` object as their first argument, return their results, and store them on the `Atoms` object under keys that start with `pyscal_`.
The pages of the [user guide](guide/neighbors) and the [descriptors](descriptors/index) explain the methods and show examples.

## Overview

| Function | Description |
|---|---|
| {py:func}`find_neighbors <pyscal3.find_neighbors>` | Find neighbors of all atoms |
| {py:func}`get_distance <pyscal3.get_distance>` | Get the distance between two positions respecting periodic boundaries |
| {py:func}`common_neighbor_analysis <pyscal3.common_neighbor_analysis>` | Calculate Common Neighbor Analysis (CNA) or Adaptive CNA |
| {py:func}`diamond_structure <pyscal3.diamond_structure>` | Identify diamond structure using extended CNA |
| {py:func}`identify_ackland_jones <pyscal3.identify_ackland_jones>` | Classify atomic environments with the Ackland–Jones method |
| {py:func}`voronoi_vector <pyscal3.voronoi_vector>` | Calculate the Voronoi structure identification vector (n3, n4, n5, n6) |
| {py:func}`steinhardt_parameter <pyscal3.steinhardt_parameter>` | Calculate Steinhardt bond order parameters q_l |
| {py:func}`wigner_w_parameter <pyscal3.wigner_w_parameter>` | Calculate the third-order Steinhardt invariant W_l |
| {py:func}`minkowski_parameter <pyscal3.minkowski_parameter>` | Calculate Minkowski structure metrics $q_l$ weighted by Voronoi face areas |
| {py:func}`angular_criteria <pyscal3.angular_criteria>` | Calculate angular criteria for diamond structure identification |
| {py:func}`chi_params <pyscal3.chi_params>` | Calculate chi-parameter vector for structure identification |
| {py:func}`find_solids <pyscal3.find_solids>` | Distinguish solid and liquid atoms |
| {py:func}`find_clusters <pyscal3.find_clusters>` | Cluster atoms based on a boolean condition |
| {py:func}`disorder <pyscal3.disorder>` | Calculate the disorder parameter |
| {py:func}`entropy <pyscal3.entropy>` | Calculate the entropy parameter for each atom |
| {py:func}`centrosymmetry <pyscal3.centrosymmetry>` | Calculate the centrosymmetry parameter |
| {py:func}`atomic_strain <pyscal3.atomic_strain>` | Calculate the atomic (Green-Lagrange) strain tensor for each atom |
| {py:func}`von_mises_strain <pyscal3.von_mises_strain>` | Compute the von Mises shear strain invariant from the atomic strain |
| {py:func}`d2min <pyscal3.d2min>` | Compute the D^2_min non-affine displacement (Falk & Langer 1998) |
| {py:func}`slip_vector <pyscal3.slip_vector>` | Compute the slip vector for each atom |
| {py:func}`wigner_seitz_analysis <pyscal3.wigner_seitz_analysis>` | Wigner-Seitz cell analysis for vacancy/interstitial detection |
| {py:func}`identify_defect_atoms <pyscal3.identify_defect_atoms>` | Identify which atoms are at vacancies, interstitials, or antisites |
| {py:func}`short_range_order <pyscal3.short_range_order>` | Calculate the Warren-Cowley short-range order parameter |
| {py:func}`coordination_number <pyscal3.coordination_number>` | Return the simple coordination number (integer neighbor count) |
| {py:func}`effective_coordination_number <pyscal3.effective_coordination_number>` | Calculate the effective coordination number (ECoN) |
| {py:func}`generalized_coordination_number <pyscal3.generalized_coordination_number>` | Calculate the generalized coordination number (GCN) |
| {py:func}`local_density <pyscal3.local_density>` | Estimate the local atomic number density |
| {py:func}`radial_distribution_function <pyscal3.radial_distribution_function>` | Calculate radial distribution function g(r) |
| {py:func}`angular_distribution_function <pyscal3.angular_distribution_function>` | Calculate the angular distribution function (ADF) |
| {py:func}`bond_length_distribution <pyscal3.bond_length_distribution>` | Calculate the bond-length distribution function (BLDF) |
| {py:func}`ace <pyscal3.ace>` | Compute Atomic Cluster Expansion (ACE) descriptors |
| {py:func}`average_over_neighbors <pyscal3.average_over_neighbors>` | Average a per-atom property over each atom's neighbors |
| {py:func}`set_num_threads <pyscal3.set_num_threads>` | Set the number of threads used by pyscal3 |
| {py:func}`get_num_threads <pyscal3.get_num_threads>` | Return the number of threads used by pyscal3 |
| [`structures`](#structure-creation) | Builders for crystals, elements, custom lattices and grain boundaries |
| {py:class}`Trajectory <pyscal3.Trajectory>` | Lazy reader for LAMMPS dump trajectories |

## Neighbors

See also [Finding neighbors](guide/neighbors).

```{eval-rst}
.. autofunction:: pyscal3.find_neighbors
```

```{eval-rst}
.. autofunction:: pyscal3.get_distance
```

## Crystal structure

See also [Common neighbor analysis](descriptors/cna).

```{eval-rst}
.. autofunction:: pyscal3.common_neighbor_analysis
```

```{eval-rst}
.. autofunction:: pyscal3.diamond_structure
```

```{eval-rst}
.. autofunction:: pyscal3.identify_ackland_jones
```

```{eval-rst}
.. autofunction:: pyscal3.voronoi_vector
```

## Orientational order

See also [Steinhardt parameters](descriptors/steinhardt).

```{eval-rst}
.. autofunction:: pyscal3.steinhardt_parameter
```

```{eval-rst}
.. autofunction:: pyscal3.wigner_w_parameter
```

```{eval-rst}
.. autofunction:: pyscal3.minkowski_parameter
```

```{eval-rst}
.. autofunction:: pyscal3.angular_criteria
```

```{eval-rst}
.. autofunction:: pyscal3.chi_params
```

## Solid and liquid

See also [Solid–liquid classification](descriptors/solid_liquid).

```{eval-rst}
.. autofunction:: pyscal3.find_solids
```

```{eval-rst}
.. autofunction:: pyscal3.find_clusters
```

```{eval-rst}
.. autofunction:: pyscal3.disorder
```

```{eval-rst}
.. autofunction:: pyscal3.entropy
```

## Defects and deformation

See also [Atomic deformation](descriptors/deformation).

```{eval-rst}
.. autofunction:: pyscal3.centrosymmetry
```

```{eval-rst}
.. autofunction:: pyscal3.atomic_strain
```

```{eval-rst}
.. autofunction:: pyscal3.von_mises_strain
```

```{eval-rst}
.. autofunction:: pyscal3.d2min
```

```{eval-rst}
.. autofunction:: pyscal3.slip_vector
```

```{eval-rst}
.. autofunction:: pyscal3.wigner_seitz_analysis
```

```{eval-rst}
.. autofunction:: pyscal3.identify_defect_atoms
```

## Chemistry, coordination and distributions

See also [Coordination](descriptors/coordination).

```{eval-rst}
.. autofunction:: pyscal3.short_range_order
```

```{eval-rst}
.. autofunction:: pyscal3.coordination_number
```

```{eval-rst}
.. autofunction:: pyscal3.effective_coordination_number
```

```{eval-rst}
.. autofunction:: pyscal3.generalized_coordination_number
```

```{eval-rst}
.. autofunction:: pyscal3.local_density
```

```{eval-rst}
.. autofunction:: pyscal3.radial_distribution_function
```

```{eval-rst}
.. autofunction:: pyscal3.angular_distribution_function
```

```{eval-rst}
.. autofunction:: pyscal3.bond_length_distribution
```

## Machine learning

See also [ACE descriptors](descriptors/ace.md).

```{eval-rst}
.. autofunction:: pyscal3.ace
```

## Utilities

See also [Large systems](guide/large_systems).

```{eval-rst}
.. autofunction:: pyscal3.average_over_neighbors
```

```{eval-rst}
.. autofunction:: pyscal3.set_num_threads
```

```{eval-rst}
.. autofunction:: pyscal3.get_num_threads
```

## Structure creation

See also [Building structures](guide/structures).

```{eval-rst}
.. autofunction:: pyscal3.structures.make_crystal
```

```{eval-rst}
.. autofunction:: pyscal3.structures.make_element
```

```{eval-rst}
.. autofunction:: pyscal3.structures.make_general_lattice
```

```{eval-rst}
.. autofunction:: pyscal3.structures.make_grain_boundary
```

```{eval-rst}
.. autofunction:: pyscal3.structures.available_structures
```

```{eval-rst}
.. autofunction:: pyscal3.structures.available_elements
```

## Trajectory

See also [Reading files and trajectories](guide/files).

```{eval-rst}
.. autoclass:: pyscal3.Trajectory
   :members: get_block, load, unload
```

```{eval-rst}
.. autoclass:: pyscal3.trajectory.Timeslice
   :members: to_atoms, to_file, to_dict
```
