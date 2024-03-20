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

In pyscal, disorder parameter can be calculated by the following code-block,

``` python
from pyscal3 import System
sys = System('conf.dump')
sys.find.neighbors(method='cutoff', cutoff=0)
q = fcc.calculate.steinhardt_parameter(6)
sys.calculate.disorder(averaged=True, q=6)
```

The value of q can be replaced with whichever. The calculated values can be accessed by, `sys.atoms.steinhardt.disorder`

## References

1. Kawasaki, T. & Onuki, A. Construction of a disorder variable from Steinhardt order parameters in binary mixtures at high densities in three dimensions. Journal of Chemical Physics 135, (2011).
