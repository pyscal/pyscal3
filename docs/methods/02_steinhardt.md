# Steinhardt's parameters

Steinhardt's bond orientational order parameters [1] are a set of parameters based on [spherical harmonics](https://en.wikipedia.org/wiki/Spherical_harmonics) to explore the local atomic environment. These parameters have been used extensively for various uses such as distinction of crystal structures, identification of solid and liquid atoms and identification of defects.

These parameters, which are rotationally and translationally invariant are defined by,

$$ 
q_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l | q_{lm}(i) |^2 \Big )^{\frac{1}{2}} 
$$

where,

$$ 
q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\pmb{r}_{ij})
$$

in which $Y_{lm}$ are the spherical harmonics and $N(i)$ is the number of neighbours of particle $i$, $\pmb{r}_{ij}$ is the vector connecting particles $i$ and $j$, and $l$ and $m$ are both intergers with $m \in [-l,+l]$. Various parameters have found specific uses, such as $q_2$ and $q_6$ for identification of crystallinity, $q_6$ for identification of solidity, and $q_4$ and $q_6$ for distinction of crystal structures [2]. Commonly this method uses a cutoff radius to identify the neighbors of an atom. 
Once the cutoff is chosen and neighbors are calculated, the calculation of Steinhardt's parameters is straightforward.

``` python
import pyscal3
from ase.build import bulk

atoms = bulk('Cu', 'fcc', cubic=True).repeat(4)
pyscal3.find_neighbors(atoms, method='cutoff', cutoff=0)
q4, q6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6])
```

Results are also stored on the `Atoms` object as `atoms.arrays['pyscal_q4']` and `atoms.arrays['pyscal_q6']`.

## Averaged Steinhardt's parameters

At high temperatures, thermal vibrations affect the atomic positions. This in turn leads to overlapping distributions of $q_l$ parameters, which makes the identification of crystal structures difficult. To address this problem, the averaged version $\bar{q}_l$ of Steinhardt's parameters was introduced by Lechner and Dellago [3]. $\bar{q}_l$ is given by,

$$
\bar{q}_l (i) =  \Big(  \frac{4\pi}{2l+1}  \sum_{m=-l}^l \Big| \frac{1}{\tilde{N}(i)} \sum_{k=0}^{\tilde{N}(i)} q_{lm}(k) \Big|^2 \Big )^{\frac{1}{2}}
$$

where the sum from $k=0$ to $\tilde{N}(i)$ is over all the neighbors and the particle itself. The averaged parameters takes into account the first neighbor shell and also information from the neighboring atoms and thus reduces the overlap between the distributions. Commonly $\bar{q}_4$ and $\bar{q}_6$ are used in identification of crystal structures.
Averaged versions can be calculated by setting the keyword `averaged=True` as follows.

``` python
aq4, aq6 = pyscal3.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
```

## Voronoi weighted Steinhardt's parameters

In order to improve the resolution of crystal structures Mickel et al [2] proposed weighting the contribution of each neighbor to the Steinhardt parameters by the ratio of the area of the Voronoi facet shared between the neighbor and host atom. The weighted parameters are given by,

$$
q_{lm} (i) =  \frac{1}{N(i)} \sum_{j=1}^{N(i)} \frac{A_{ij}}{A} Y_{lm}(\pmb{r}_{ij})
$$

where $A_{ij}$ is the area of the Voronoi facet between atoms $i$ and $j$ and $A$ is the sum of the face areas of atom $i$. The dedicated function `minkowski_parameter` performs this calculation; see the [Minkowski parameters](09_minkowski.md) page for details.

``` python
pyscal3.find_neighbors(atoms, method='voronoi')
q4, q6 = pyscal3.minkowski_parameter(atoms, l=[4, 6])
```

## References

1. Steinhardt, P. J., Nelson, D. R. & Ronchetti, M. Bond-orientational order in liquids and glasses. Physical Review B 28, 784–805 (1983).
2. Mickel, W., Kapfer, S. C., Schröder-Turk, G. E. & Mecke, K. Shortcomings of the bond orientational order parameters for the analysis of disordered particulate matter. Journal of Chemical Physics 138, (2013).
3. Lechner, W. & Dellago, C. Accurate determination of crystal structures based on averaged local bond order parameters. Journal of Chemical Physics 129, (2008).
4. Haeberle, J., Sperl, M. & Born, P. Distinguishing noisy crystalline structures using bond orientational order parameters. Eur. Phys. J. E 42, 149 (2019).



