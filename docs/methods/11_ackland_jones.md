# Ackland-Jones structural classification

Building on the $\chi$ parameters described in the [angular parameters](05_angular.md) page, Ackland and Jones [1] proposed a decision tree on the histogram of bond angle cosines that assigns each atom to one of fcc, hcp, bcc, icosahedral, or "unknown". The classifier is parameter free once the neighbor list is fixed; the same set of $\chi$ bins,

$$
\chi_0, \chi_1, \ldots, \chi_7
$$

are used as features, and a small set of empirical inequalities then maps the bin counts to a structure label. The classifier is robust at finite temperatures because it relies on relative populations of angular bins rather than on individual bond lengths.

In pyscal3, the classification can be performed with,

``` python
import pyscal3
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
labels = pyscal3.identify_ackland_jones(atoms)
```

The integer label per atom is stored as `atoms.arrays['pyscal_ackland_label']` and the corresponding string name (`'fcc'`, `'bcc'`, `'hcp'`, `'ico'`, `'unknown'`) as `atoms.arrays['pyscal_structure']`.

## References

1. Ackland, G. J. & Jones, A. P. Applications of local crystal structure measures in experiment and simulation. Phys. Rev. B 73, 054104 (2006).
