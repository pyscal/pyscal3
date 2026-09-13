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

Finding solid atoms in liquid starts with reading in a file and calculating neighbors.

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=4)
```

Once again, there are various methods for finding neighbors. Once the neighbors are calculated, solid atoms can be found directly by,

``` python
largest = pyscal.find_solids(atoms, bonds=6, threshold=0.5,
                              avgthreshold=0.6, cluster=True)
```

`bonds` sets the number of minimum bonds a particle should have (as defined above), `threshold` and `avgthreshold` are the same quantities that appear in the equations above. Setting the keyword `cluster` to `True` returns the size of the largest solid cluster. The per-atom solid/liquid label is stored as `atoms.arrays['pyscal_solid']`.

The intermediate quantities are stored as well, so the distribution of bond correlations can be inspected directly: $s_{ij}$ for every atom and each of its neighbors as `pyscal_sij` (in `atoms.arrays` when all atoms have the same number of neighbors, otherwise in `atoms.info`), the per-atom average $\langle s_{ij} \rangle$ as `atoms.arrays['pyscal_avg_sij']`, and the number of solid bonds as `atoms.arrays['pyscal_bonds']`.

``` python
import numpy as np

pyscal.find_solids(atoms, bonds=6, threshold=0.5, avgthreshold=0.6, cluster=False)
sij = atoms.arrays.get('pyscal_sij', atoms.info.get('pyscal_sij'))
all_sij = np.concatenate([np.asarray(row) for row in sij])   # every i-j pair
hist, edges = np.histogram(all_sij, bins=50, range=(-0.5, 1.0))
```

Clustering can use a different cutoff than the neighbor search through the `cutoff` keyword of `find_solids` (or `find_clusters`); by default the neighbor cutoff is used.

## References

1. Auer, S. & Frenkel, D. Numerical Simulation of Crystal Nucleation in Colloids. in Advanced Computer Simulation: Approaches for Soft Matter Sciences I (eds. Dr. Holm, C. & Prof. Dr. Kremer, K.) 149–208 (Springer Berlin Heidelberg, Berlin, Heidelberg, 2005). doi:10.1007/b99429.
2. Bokeloh, J., Wilde, G., Rozas, R. E., Benjamin, R. & Horbach, J. Nucleation barriers for the liquid-to-crystal transition in simple metals: Experiment vs. simulation. European Physical Journal: Special Topics 223, 511–526 (2014).

