# Coordination measures

The coordination number is the simplest and most widely used local descriptor of an atom's environment. pyscal exposes three variants that progressively relax the assumption of a sharp cutoff.

## Coordination number

Given a neighbor list, the coordination number $\mathrm{CN}(i)$ is just the number of neighbors of atom $i$. It is meaningful when neighbors come from a well-defined criterion such as a cutoff at the first $g(r)$ minimum.

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=3.5)
cn = pyscal.coordination_number(atoms)
```

Per-atom values are stored as `atoms.arrays['pyscal_cn']`.

## Effective coordination number

Hoppe [1] introduced an *effective* coordination number that smoothly weights each neighbor by its distance, removing the discontinuity of a hard cutoff. With $r_{ij}$ the bond length to neighbor $j$ and a self-consistently defined weighted average distance $\bar{r}(i)$,

$$
\mathrm{ECoN}(i) = \sum_{j} \exp \bigg(1 - \bigg(\frac{r_{ij}}{\bar{r}(i)}\bigg)^6 \bigg)
$$

with

$$
\bar{r}(i) = \frac{\sum_j r_{ij} \exp(1 - (r_{ij}/r_{\min}(i))^6)}{\sum_j \exp(1 - (r_{ij}/r_{\min}(i))^6)}
$$

iterated to convergence. ECoN reproduces the integer coordination number for ideal lattices and varies smoothly otherwise.

``` python
econ = pyscal.effective_coordination_number(atoms)
```

Per-atom values are stored as `atoms.arrays['pyscal_econ']`.

## Generalized coordination number

The generalized coordination number $\overline{\mathrm{CN}}(i)$ of Calle-Vallejo et al. [2] weights each neighbor by its own coordination,

$$
\overline{\mathrm{CN}}(i) = \sum_{j \in \mathrm{neigh}(i)} \frac{\mathrm{CN}(j)}{\mathrm{CN}_{\max}}
$$

where $\mathrm{CN}_{\max}$ is the bulk coordination of the underlying lattice (12 for fcc, 8 for bcc, etc.). It is widely used as a descriptor for catalytic activity at metal surfaces.

``` python
gcn = pyscal.generalized_coordination_number(atoms, cn_max=12)
```

Per-atom values are stored as `atoms.arrays['pyscal_gcn']`.

## References

1. Hoppe, R. The Coordination Number – an "Inorganic Chameleon". Angewandte Chemie International Edition 9, 25–34 (1970).
2. Calle-Vallejo, F. et al. Finding optimal surface sites on heterogeneous catalysts by counting nearest neighbors. Science 350, 185–189 (2015).
