# Atomic Cluster Expansion descriptors

The Atomic Cluster Expansion (ACE) of Drautz [1] provides a complete and systematically improvable basis for the local atomic environment. It is the foundation of several modern interatomic potentials and is increasingly used as a general-purpose feature set for machine learning of structure-property relations.

For each atom $i$ a set of two-body atomic basis functions is constructed,

$$
A_{nlm}(i) = \sum_{j \in \mathrm{neigh}(i)} R_{nl}(r_{ij})\, Y_{lm}(\hat{\pmb{r}}_{ij})
$$

where $R_{nl}$ are radial functions truncated by a cutoff and $Y_{lm}$ are spherical harmonics. Rotation-invariant features at body order $\nu$ are then built by contracting $\nu$ such atomic basis functions with appropriate generalised Clebsch-Gordan coefficients,

$$
B^{(\nu)}_{nl}(i) = \sum_{m_1 + \cdots + m_\nu = 0} C^{l_1 \cdots l_\nu}_{m_1 \cdots m_\nu} \prod_{k=1}^{\nu} A_{n_k l_k m_k}(i)
$$

The first three body orders give two-body ($\nu=1$), three-body ($\nu=2$), and four-body ($\nu=3$) descriptors. Increasing $n_{\max}$, $l_{\max}$, and $\nu_{\max}$ yields a systematically convergent representation.

In pyscal, the ACE descriptors up to $\nu=3$ can be calculated by,

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=5.0)
descriptors = pyscal.ace(atoms, nmax=4, lmax=4, nu_max=2, cutoff=5.0)
```

The neighbor list must be computed first; the radial basis is truncated at `cutoff` (by default the cutoff used for the neighbor list). The $\nu=3$ block couples three atomic basis functions with Wigner $3j$ symbols, so all blocks are rotationally invariant. The full descriptor matrix is stored as `atoms.arrays['pyscal_ace']` with shape $(N, n_{\mathrm{features}})$ and the individual blocks are returned in a dictionary (`'nu1'`, `'nu2'`, `'nu3'`, `'full'`). Hyperparameters used in the calculation are recorded in `atoms.info['pyscal_ace_params']`. Setting `normalize=False` disables the per-atom normalisation of the descriptor vector to unit length; choosing `nu_max=1` returns only the radial $\nu=1$ block.

## References

1. Drautz, R. Atomic cluster expansion for accurate and transferable interatomic potentials. Phys. Rev. B 99, 014104 (2019).
