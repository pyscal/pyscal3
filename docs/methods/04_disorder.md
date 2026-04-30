# Disorder parameter

Kawasaki and Onuki [1] proposed a disorder variable based on Steinhardt's order paramaters which can be used to distinguish between ordered and disordered structures

The disorder variable for an atom is defined as,

$$
D_j = \frac{1}{n_b^j} \sum_{k \in neighbors } [S_{jj} + S_{kk} - 2S_{jk}]
$$

where S is given by,

$$
S_{jk} = \sum_{-l \leq m \leq l} q_{lm}^j (q_{lm}^k)^*
$$

l = 6 was used in the original publication as it is a good indicator of crystallinity. However, l = 4 can also be used for treating bcc structures. An averaged disorder parameter for each atom can also be calculated in pyscal,

$$
\bar{D}_j = \frac{1}{n_b^j} \sum_{k \in neighbors } D_j
$$

In pyscal3, the disorder parameter can be calculated by the following code-block,

``` python
import pyscal3
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
pyscal3.steinhardt_parameter(atoms, l=6)
pyscal3.disorder(atoms, q=6, averaged=True)
```

The value of `q` can be set to any integer for which Steinhardt's parameters have been computed. The per-atom disorder is stored as `atoms.arrays['pyscal_disorder']` (and `pyscal_avg_disorder` for the averaged variant).

## References

1. Kawasaki, T. & Onuki, A. Construction of a disorder variable from Steinhardt order parameters in binary mixtures at high densities in three dimensions. Journal of Chemical Physics 135, (2011).
