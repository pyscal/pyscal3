# Common neighbor analysis

Common neighbor analysis (CNA) [1, 2] classifies the local environment of an atom from the bonding topology of its neighbors rather than from bond angles or orientational order parameters. For every pair of an atom $i$ and one of its neighbors $j$, three integers are recorded:

- the number of neighbors shared by $i$ and $j$ (the *common neighbors*),
- the number of bonds among those common neighbors, and
- the length of the longest chain formed by these bonds.

Each neighbor pair thus receives a signature such as 421 or 555. The multiset of signatures around an atom identifies the crystal structure:

| Structure | Signatures |
| --- | --- |
| fcc | 12 × 421 |
| hcp | 6 × 421 and 6 × 422 |
| bcc | 8 × 666 and 6 × 444 |
| icosahedral | 12 × 555 |

Atoms whose signatures match none of these patterns are labelled *other*.

## Conventional CNA

In conventional CNA a fixed cutoff decides which atoms are bonded. For fcc and hcp the cutoff lies between the first and second neighbor shell, $r_c = \tfrac{1}{2}\left(\tfrac{1}{\sqrt{2}} + 1\right) a \approx 0.854\,a$, and for bcc between the second and third shell, $r_c = \tfrac{1}{2}(1 + \sqrt{2})\,a \approx 1.207\,a$, where $a$ is the lattice constant. This works well for a single perfect crystal but breaks down when the lattice constant varies, for example under strain or in a two-phase system.

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
counts = pyscal.common_neighbor_analysis(atoms, lattice_constant=4.05)
```

## Adaptive CNA

Adaptive CNA [3] removes the lattice-constant dependence by choosing the cutoff per atom from its own environment. The 12 (fcc/hcp) or 14 (bcc) nearest neighbors are collected, and the cutoff is placed between the shells using the mean neighbor distance,

$$
r_c^{\mathrm{fcc}}(i) = \frac{1 + \sqrt{2}}{2}\, \frac{1}{12}\sum_{j=1}^{12} |\pmb{r}_{ij}|,
\qquad
r_c^{\mathrm{bcc}}(i) = \frac{1 + \sqrt{2}}{2}\, \frac{1}{14}\left(\frac{2}{\sqrt{3}}\sum_{j=1}^{8} |\pmb{r}_{ij}| + \sum_{j=9}^{14} |\pmb{r}_{ij}|\right).
$$

The fcc/hcp test is performed first with the 12-neighbor cutoff; atoms that remain unclassified are tested for bcc with the 14-neighbor cutoff. Adaptive CNA is used when no lattice constant is given and is the recommended mode.

``` python
counts = pyscal.common_neighbor_analysis(atoms)
```

`common_neighbor_analysis` finds its own neighbors, so `find_neighbors` does not have to be called first (an existing neighbor list is left untouched). It returns a dictionary with the number of atoms per structure and stores the per-atom label in `atoms.arrays['pyscal_structure']`: 0 other, 1 fcc, 2 hcp, 3 bcc, 4 icosahedral.

## Diamond structures

Cubic and hexagonal diamond lattices have only four nearest neighbors and are not resolved by the signatures above. pyscal follows the approach of Maras et al. [4]: the second-neighbor shell (the 12 neighbors of neighbors) of each atom is classified with the fcc/hcp CNA signatures, which distinguishes cubic diamond (fcc-like second shell) from hexagonal diamond (hcp-like second shell). Atoms whose second shell is imperfect but that neighbor a diamond atom are labelled as first or second neighbors of the respective diamond type.

``` python
counts = pyscal.diamond_structure(atoms)
```

The label codes stored in `atoms.arrays['pyscal_structure']` are 0 other, 1 cubic diamond, 2 and 3 first and second neighbors of cubic diamond, 4 hexagonal diamond, 5 and 6 first and second neighbors of hexagonal diamond.

## References

1. Honeycutt, J. D. & Andersen, H. C. Molecular dynamics study of melting and freezing of small Lennard-Jones clusters. J. Phys. Chem. 91, 4950–4963 (1987).
2. Faken, D. & Jónsson, H. Systematic analysis of local atomic structure combined with 3D computer graphics. Comput. Mater. Sci. 2, 279–286 (1994).
3. Stukowski, A. Structure identification methods for atomistic simulations of crystalline materials. Modelling Simul. Mater. Sci. Eng. 20, 045021 (2012).
4. Maras, E., Trushin, O., Stukowski, A., Ala-Nissila, T. & Jónsson, H. Global transition path search for dislocation formation in Ge on Si(001). Comput. Phys. Commun. 205, 13–21 (2016).
