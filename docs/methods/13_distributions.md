# Angular and bond length distributions

The angular distribution function (ADF) and the bond length distribution function (BLDF) are the two simplest two- and three-body distributions one can build from a neighbor list. They are diagnostic fingerprints of the local environment and complement the radial distribution function $g(r)$.

## Angular distribution function

For each central atom, every pair of neighbors $(j, k)$ defines a bond angle

$$
\theta_{ijk} = \arccos \bigg( \frac{\pmb{r}_{ij} \cdot \pmb{r}_{ik}}{|\pmb{r}_{ij}|\, |\pmb{r}_{ik}|} \bigg)
$$

The ADF is the histogram of $\theta_{ijk}$ over all atoms and all neighbor pairs.

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=3.5)
hist, angles = pyscal.angular_distribution_function(atoms, bins=180)
```

The histogram and bin centers are also stored as `atoms.info['pyscal_adf']` and `atoms.info['pyscal_adf_angles']`.

## Bond length distribution function

The BLDF is the distribution of pair distances $r_{ij}$ restricted to the bonds in the current neighbor list. It differs from $g(r)$ in that it carries no shell-volume normalisation, so it more directly reflects the geometry of the chosen coordination shell.

``` python
hist, r = pyscal.bond_length_distribution(atoms, bins=100)
```

The histogram and bin centers are stored as `atoms.info['pyscal_bldf']` and `atoms.info['pyscal_bldf_r']`.
