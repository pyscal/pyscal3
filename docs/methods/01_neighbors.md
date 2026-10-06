# Methods to calculate neighbors of a particle

pyscal includes different methods to explore the local environment of a particle that rely on the calculation of nearest neighbors. Various approaches to compute the neighbors of particles are discussed here.

## Fixed cutoff method

The most common method to calculate the nearest neighbors of an atom is using a cutoff radius. Commonly, a cutoff is selected as the first minimum of the radial distribution functions. Once a cutoff is selected, the neighbors of an atom are those that fall within this selected radius. The following code snippet will use the cutoff method to calculate neighbors. In this example, `conf.dump` is assumed to
be the input configuration of the system. A cutoff radius of 3 is assumed for calculation of neighbors.

``` python
import pyscal
from ase.io import read

atoms = read('conf.dump', format='lammps-dump-text')
pyscal.find_neighbors(atoms, method='cutoff', cutoff=3)
```

## Adaptive cutoff methods

A fixed cutoff radius can introduce limitations to explore the local environment of the particle in some cases:

-   At finite temperatures, when thermal fluctuations take place, the selection of a fixed cutoff may result in an inaccurate description of the local environment.
-   If there is more than one structure present in the system, for example, bcc and fcc, the selection of cutoff such that it includes the first shell of both structures can be difficult.

In order to achieve a more accurate description of the local environment, various adaptive approaches have been proposed. Two of the methods implemented in the module are discussed below.

### Solid angle based nearest neighbor algorithm (SANN)

SANN algorithm [1] determines the cutoff radius by counting the solid angles around an atom and equating it to $4\pi$. The algorithm solves the following equation iteratively.

$$
R_i^{(m)} = \frac{\sum_{j=1}^m r_{i,j}}{m-2} < r_{i, m+1}
$$

where $i$ is the host atom, $j$ are its neighbors with $r_{ij}$ is the distance between atoms $i$ and $j$. $R_i$ is the cutoff radius for each particle $i$ which is found by increasing the neighbor of neighbors $m$ iteratively. For a description of the algorithm and more details, please check the reference [1]. SANN algorithm can be used to find the neighbors by,

``` python
pyscal.find_neighbors(atoms, method='cutoff', cutoff='sann')
```

Since SANN algorithm involves sorting, a sufficiently large cutoff is used in the beginning to reduce the number entries to be sorted. This parameter is calculated by,

$$
r_{initial} = \mathrm{threshold} \times \bigg(\frac{\mathrm{Simulation~box~volume}}{\mathrm{Number~of~particles}}\bigg)^{\frac{1}{3}}
$$

a tunable `threshold` parameter can be set through function arguments.

### Adaptive cutoff method

An adaptive cutoff specific for each atom can also be found using an algorithm similar to adaptive common neighbor analysis [2]. This adaptive cutoff is calculated by first making a list of all neighbor
distances for each atom similar to SANN method. Once this list is available, then the cutoff is calculated from,

$$
r_{cut}(i) = \mathrm{padding}\times \bigg(\frac{1}{\mathrm{nlimit}} \sum_{j=1}^{\mathrm{nlimit}} r_{ij} \bigg)
$$

This method can be chosen by,

``` python
pyscal.find_neighbors(atoms, method='cutoff', cutoff=0)
```

The keyword `cutoff=0` selects the adaptive method. The `padding` and `nlimit` parameters in the above equation can be tuned using the respective keywords.

Either of the adaptive method can be used to find neighbors, which can then be used to calculate Steinhardt\'s parameters or their averaged version.

## Voronoi tessellation

[Voronoi tessellation](https://en.wikipedia.org/wiki/Voronoi_diagram) provides a completely parameter free geometric approach for calculation of neighbors. [Voro++](http://math.lbl.gov/voro++/) code is used for Voronoi tessellation. Neighbors can be calculated using this method by,

``` python
pyscal.find_neighbors(atoms, method='voronoi')
```

Finding neighbors using Voronoi tessellation also calculates a weight for each neighbor. The weight of a neighbor $j$ towards a host atom $i$ is given by,

$$
W_{ij} = \frac{A_{ij}}{\sum_{j=1}^N A_{ij}}
$$

where $A_{ij}$ is the area of Voronoi facet between atom $i$ and $j$, $N$ are all the neighbors identified through Voronoi tessellation. This weight can be used later for calculation of weighted Steinhardt's
parameters. Optionally, it is possible to choose the exponent for this weight. Option `voroexp` is used to set this option. For example if `voroexp=2`, the weight would be calculated as,

$$
W_{ij} = \frac{A_{ij}^2}{\sum_{j=1}^N A_{ij}^2}
$$

``` python
pyscal.find_neighbors(atoms, method='voronoi', voroexp=2)
```

## How neighbors are searched and stored

The fixed-cutoff, shell, adaptive, SANN and number searches use the [matscipy-neighbours](https://github.com/libAtoms/matscipy-neighbours) library, which is included in pyscal and runs on one thread. Voronoi neighbors come from Voro++.

- **Fixed cutoff.** An atom $j$ is a neighbor of $i$ if $r_{ij} < r_{cut}$.
- **Shell** (`shell_thickness > 0`). The condition is $r_{cut} \le r_{ij} \le r_{cut} + \mathrm{shell\_thickness}$.
- **Adaptive, SANN and number.** The candidates are the atoms with $r_{ij} \le r_{initial}$, with $r_{initial}$ as defined above. They are sorted by distance, and candidates at the same distance (to within $10^{-10}$) are sorted by atom index. The number method (`method='number'`) keeps the first `nmax` candidates, so in a perfect crystal the choice among equidistant neighbors is reproducible.
- **Periodic images.** When the cutoff exceeds half the width of the cell, several periodic images of the same atom, including the atom itself, can be neighbors. Each image is a separate entry, with its own distance and vector.
- **The `cells` argument** of `find_neighbors` is ignored. It is kept so that existing scripts keep working.

The results are stored in the `Atoms` object:
- When every atom has the same number of neighbors $k$, `atoms.arrays["pyscal_neighbors"]` is an $(n, k)$ integer array, and the distances, weights and angles are $(n, k)$ float arrays.
- When the numbers differ, the same keys are stored in `atoms.info` as lists of lists.
- The neighbor vectors $\mathbf{r}_i - \mathbf{r}_j$ are always stored in `atoms.info["pyscal_diff"]`.
- If no atom has a neighbor, all these keys are $(n, 0)$ arrays in `atoms.arrays`.

## References

1. van Meel, J. A., Filion, L., Valeriani, C. & Frenkel, D. A parameter-free, solid-angle based, nearest- neighbor algorithm. J Chem Phys 234107, (2012).
2. Stukowski, A. Structure identification methods for atomistic simulations of crystalline materials. Modelling and Simulation in Materials Science and Engineering 20, (2012).
