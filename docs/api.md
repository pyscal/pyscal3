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
