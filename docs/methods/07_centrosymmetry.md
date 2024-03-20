# Centrosymmetry parameter

Centrosymmetry parameter (CSP) was introduced by Kelchner et al. [1] to identify defects in crystals. The parameter measures the loss of local symmetry. For an atom with $N$ nearest neighbors, the parameter is given by,

$$
\mathrm{CSP} = \sum_{i=1}^{N/2} \big | \textbf{r}_i + \textbf{r}_{i+N/2} \big |^2
$$

$\textbf{r}_i$ and $\textbf{r}_{i+N/2}$ are vectors from the central atom to two opposite pairs of neighbors. There are two main methods to identify the opposite pairs of neighbors as described in [this publication](https://arxiv.org/abs/2003.08879). The first of the approaches is called Greedy Edge Selection (GES) [2]
and is implemented in [LAMMPS](https://lammps.sandia.gov/) and [Ovito](https://www.ovito.org/). GES algorithm calculates a weight $w_{ij} = |\textbf{r}_i + \textbf{r}_j|$ for all combinations of neighbors around an atom and calculates CSP over the smallest $N/2$ weights.

A centrosymmetry parameter calculation using GES algorithm can be carried out as follows-

``` python
from pyscal3 import System
sys = System('conf.dump')
csm = sys.calculate.centrosymmetry(nmax = 12)
```

`nmax` parameter specifies the number of nearest neighbors to be considered for the calculation of CSP. 

## References

1. Kelchner, C. L., Plimpton, S. J. & Hamilton, J. C. Dislocation nucleation and defect structure during surface indentation. Phys. Rev. B 58, 11085–11088 (1998).
2. Stukowski, A. Structure identification methods for atomistic simulations of crystalline materials. Modelling and Simulation in Materials Science and Engineering 20, (2012).

