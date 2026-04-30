# Minkowski structure metrics

The standard Steinhardt parameters weight all neighbors equally, which makes them sensitive to the precise choice of cutoff. Mickel et al. [1] proposed a parameter-free variant in which each neighbor contribution is weighted by the relative area of its Voronoi facet,

$$
q_l(i) = \sqrt{\frac{4\pi}{2l+1} \sum_{m=-l}^{l} \bigg|\sum_{j=1}^{N(i)} \frac{A_{ij}}{A(i)} Y_{lm}(\theta_{ij}, \phi_{ij}) \bigg|^2}
$$

where $A_{ij}$ is the facet area shared by atoms $i$ and $j$, and $A(i) = \sum_j A_{ij}$. The resulting *Minkowski structure metrics* are continuous functions of the atomic positions and remove the discontinuities that affect cutoff-based Steinhardt parameters near coordination changes.

Higher powers of the area weight may be used as discussed in [1]; the exponent is set with the `voroexp` keyword. Neighbors must be obtained through the Voronoi method.

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='voronoi')
q4, q6 = pyscal.minkowski_parameter(atoms, l=[4, 6])
```

Averaged versions in the spirit of Lechner and Dellago can be obtained by setting `averaged=True`. Per-atom values are stored as `atoms.arrays['pyscal_q4']`, `atoms.arrays['pyscal_q6']`, etc.

## References

1. Mickel, W., Kapfer, S. C., Schröder-Turk, G. E. & Mecke, K. Shortcomings of the bond orientational order parameters for the analysis of disordered particulate matter. The Journal of Chemical Physics 138, 044501 (2013).
