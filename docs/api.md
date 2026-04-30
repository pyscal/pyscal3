# API Reference

All public functions are available directly from `import pyscal3`.
Results are stored on the ASE `Atoms` object with the `pyscal_` prefix.

## Neighbor Finding

```{eval-rst}
.. autofunction:: pyscal3.find_neighbors
```

```{eval-rst}
.. autofunction:: pyscal3.get_distance
```

## Structural Descriptors

### Steinhardt Parameters

```{eval-rst}
.. autofunction:: pyscal3.steinhardt_parameter
```

### Wigner $W_l$ Parameters

```{eval-rst}
.. autofunction:: pyscal3.wigner_w_parameter
```

### Minkowski Structure Metrics

```{eval-rst}
.. autofunction:: pyscal3.minkowski_parameter
```

### Disorder Parameter

```{eval-rst}
.. autofunction:: pyscal3.disorder
```

### Common Neighbor Analysis

```{eval-rst}
.. autofunction:: pyscal3.common_neighbor_analysis
```

### Diamond Structure

```{eval-rst}
.. autofunction:: pyscal3.diamond_structure
```

### Centrosymmetry Parameter

```{eval-rst}
.. autofunction:: pyscal3.centrosymmetry
```

### Voronoi Vector

```{eval-rst}
.. autofunction:: pyscal3.voronoi_vector
```

### Entropy Parameter

```{eval-rst}
.. autofunction:: pyscal3.entropy
```

### Short-Range Order

```{eval-rst}
.. autofunction:: pyscal3.short_range_order
```

### Radial Distribution Function

```{eval-rst}
.. autofunction:: pyscal3.radial_distribution_function
```

### Angular Criteria

```{eval-rst}
.. autofunction:: pyscal3.angular_criteria
```

### Chi Parameters

```{eval-rst}
.. autofunction:: pyscal3.chi_params
```

### Ackland-Jones Classification

```{eval-rst}
.. autofunction:: pyscal3.identify_ackland_jones
```

### Coordination Numbers

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

### Angular and Bond Length Distributions

```{eval-rst}
.. autofunction:: pyscal3.angular_distribution_function
```

```{eval-rst}
.. autofunction:: pyscal3.bond_length_distribution
```

### Atomic Deformation

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

### Wigner-Seitz Defect Analysis

```{eval-rst}
.. autofunction:: pyscal3.wigner_seitz_analysis
```

```{eval-rst}
.. autofunction:: pyscal3.identify_defect_atoms
```

### ACE Descriptors

```{eval-rst}
.. autofunction:: pyscal3.ace
```

## Solid/Liquid Classification

```{eval-rst}
.. autofunction:: pyscal3.find_solids
```

```{eval-rst}
.. autofunction:: pyscal3.find_clusters
```

## Utilities

```{eval-rst}
.. autofunction:: pyscal3.average_over_neighbors
```

## Structure Creation

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

```{eval-rst}
.. autoclass:: pyscal3.Trajectory
   :members:
   :undoc-members:
```
