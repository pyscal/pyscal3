# Angular parameters

## Angular criteria for identification of diamond structure

Angular parameter introduced by Uttormark et al is used to measure the tetrahedrality of local atomic structure. An atom belonging to diamond structure has four nearest neighbors which gives rise to six three body angles around the atom. The angular parameter $A$ is then defined as,

$$ 
A = \sum_{i=1}^6 (\cos(\theta_i)+\frac{1}{3})^2
$$

An atom belonging to a diamond structure shows an angular parameter close to 0. It can be calculated in pyscal3 with -

``` python
import pyscal3
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
pyscal3.angular_criteria(atoms)
```

The per-atom value is stored as `atoms.arrays['pyscal_angular']`.

## $\chi$ parameters for structural identification

$\chi$ parameters introduced by Ackland and Jones [1] measures all local angles created by an atom with its neighbors and creates a histogram of these angles to produce vector which can be used to identify structures. After finding the neighbors of an atom, $\cos \theta_{ijk}$ for atoms j and k which are neighbors of i is calculated for all combinations of i, j and k. The set of all calculated cosine values are then added to a histogram with the following bins - \[-1.0, -0.945, -0.915, -0.755, -0.705, -0.195, 0.195, 0.245, 0.795, 1.0\]. Compared to $\chi$ parameters from $\chi_0$ to $\chi_7$ in the associated publication, the vector calculated in pyscal contains values from $\chi_0$ to $\chi_8$ which is due to an additional $\chi$ parameter which measures the number of neighbors between cosines -0.705 to -0.195. The $\chi$ vector is characteristic of the local atomic environment and can be used to identify crystal structures, details of which can be found in the publication[1].

$\chi$ parameters can be calculated in pyscal3 using,

``` python
import pyscal3
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
pyscal3.chi_params(atoms)
```

The per-atom $\chi$ vector is stored as `atoms.arrays['pyscal_chiparams']`. For a direct structural label (FCC/BCC/HCP/ICO/other) built from $\chi$ values, see [Ackland-Jones classification](11_ackland_jones.md).

## References

1. Ackland, G. J. & Jones, A. P. Applications of local crystal structure measures in experiment and simulation. Physical Review B - Condensed Matter and Materials Physics 73, 1–7 (2006).
