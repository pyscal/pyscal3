# Classification of atoms as solid or liquid


pyscal can also be used to distinguish solid and liquid atoms. The classification is based on Steinhardt's parameters,
specifically $q_6$. The method defines two neighboring atoms $i$ and $j$ as having solid bonds if a parameter $s_{ij}$ [1],

$$
s_{ij} = \sum_{m=-6}^6 q_{6m}(i) q_{6m}^*(j) \geq \mathrm{threshold}
$$

Additionally, a second order parameter is used to improve the distinction in solid-liquid boundaries [2]. This is defined by the criteria,

$$
\langle s_{ij} \rangle > \mathrm{avgthreshold}
$$

If a particle has $n$ number of bonds with $s_{ij} \geq \mathrm{threshold}$ and the above condition is also satisfied, it is considered as a solid. The solid atoms can be clustered to find the largest solid cluster of atoms. 

Finding solid atoms in liquid start with reading in a file and calculation of neighbors.

``` python
from pyscal3 import System
sys = System('conf.dump')
sys.find.neighbors(method='cutoff', cutoff=4)
```

Once again, there are various methods for finding neighbors. Once the neighbors are calculated, solid atoms can be found directly by,

``` python
sys.find.solids(bonds=6, threshold=0.5, avgthreshold=0.6, cluster=True)
```

`bonds` set the number of minimum bonds a particle should have (as defined above), `threshold` and `avgthreshold` are the same quantities that appear in the equations above. Setting the keyword `cluster` to True returns the size of the largest solid cluster. It is also possible to check if each atom is solid or not.

``` python
sys.atoms.solid
```

## References

1. Auer, S. & Frenkel, D. Numerical Simulation of Crystal Nucleation in Colloids. in Advanced Computer Simulation: Approaches for Soft Matter Sciences I (eds. Dr. Holm, C. & Prof. Dr. Kremer, K.) 149–208 (Springer Berlin Heidelberg, Berlin, Heidelberg, 2005). doi:10.1007/b99429.
2. Bokeloh, J., Wilde, G., Rozas, R. E., Benjamin, R. & Horbach, J. Nucleation barriers for the liquid-to-crystal transition in simple metals: Experiment vs. simulation. European Physical Journal: Special Topics 223, 511–526 (2014).

