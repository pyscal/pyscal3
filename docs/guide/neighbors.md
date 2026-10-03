---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Finding neighbors

Almost every descriptor in pyscal is computed from the neighbors of each atom.
The neighbor list therefore affects every result, and choosing it is the most important decision in an analysis.
This page describes the methods of `find_neighbors`, compares them on a crystal and a liquid, and explains how the neighbors are stored.

```{code-cell} ipython3
:tags: [remove-cell]
import os, sys, warnings
sys.path.insert(0, os.path.abspath(".."))
from _plotstyle import EXAMPLES, COLOURS, METHODS, DARK, figure, label, note
os.chdir(EXAMPLES)
warnings.simplefilter("ignore")
%config InlineBackend.figure_format = "retina"
```

## The methods

| Method | Call | Neighbors of atom $i$ |
|---|---|---|
| Fixed cutoff | `find_neighbors(atoms, method="cutoff", cutoff=3.5)` | all atoms with $r_{ij} < r_c$ |
| Shell | `find_neighbors(atoms, method="cutoff", cutoff=3.5, shell_thickness=0.5)` | all atoms with $r_c \le r_{ij} \le r_c + \Delta r$ |
| Adaptive cutoff | `find_neighbors(atoms, method="cutoff", cutoff=0)` | all atoms within a cutoff $r_c(i)$ computed for each atom |
| SANN | `find_neighbors(atoms, method="cutoff", cutoff="sann")` | the nearest atoms whose solid angles add up to $4\pi$ |
| Number | `find_neighbors(atoms, method="number", nmax=12)` | the `nmax` nearest atoms |
| Voronoi | `find_neighbors(atoms, method="voronoi")` | all atoms that share a face of the Voronoi cell of $i$ |

Here $r_{ij}$ is the distance between atoms $i$ and $j$, using the nearest periodic image.

### Fixed cutoff

An atom $j$ is a neighbor of $i$ if $r_{ij} < r_c$.
The cutoff $r_c$ is usually placed at the first minimum of the radial distribution function $g(r)$, between the first and the second neighbor shell.
With `shell_thickness` $\Delta r > 0$, only the atoms in the shell $r_c \le r_{ij} \le r_c + \Delta r$ are neighbors.

### Adaptive cutoff

`cutoff=0` gives each atom its own cutoff, from the distances to its `nlimit` nearest atoms [1]:

$$
r_c(i) = \mathrm{padding} \times \frac{1}{\mathrm{nlimit}} \sum_{j=1}^{\mathrm{nlimit}} r_{ij},
$$

where the sum runs over the `nlimit` nearest atoms, sorted by distance.
The defaults are `nlimit=6` and `padding=1.2`.
The cutoff follows local changes of the density, for example near an interface or in a temperature gradient.

### SANN

The solid angle based nearest neighbor (SANN) method [2] has no parameter.
It adds the nearest atoms one by one, and stops at the smallest number $m$ for which

$$
R_i^{(m)} = \frac{1}{m - 2} \sum_{j=1}^{m} r_{ij} < r_{i, m+1},
$$

with the distances $r_{ij}$ sorted in increasing order.
The $m$ atoms closer than $R_i^{(m)}$ are the neighbors.

### Number

`method="number"` takes the `nmax` nearest atoms.
Atoms at the same distance, to within $10^{-10}$, are ordered by their index, so that in a perfect crystal the choice is the same on every run.

### Voronoi

`method="voronoi"` uses the Voronoi tessellation computed by [Voro++](https://math.lbl.gov/voro++/).
Two atoms are neighbors if their Voronoi cells share a face.
Each neighbor gets a weight from the area $A_{ij}$ of the shared face:

$$
w_{ij} = \frac{A_{ij}^{\,p}}{\sum_{k} A_{ik}^{\,p}},
$$

where the sum runs over all neighbors $k$ of $i$ and the exponent $p$ is set by `voroexp` (default 1).
Steinhardt parameters use these weights, which makes them less sensitive to the small faces that thermal motion creates (see [Minkowski structure metrics](../descriptors/minkowski)).

### Candidate radius

The adaptive, SANN and number methods first collect candidates within the radius

$$
r_\mathrm{initial} = \mathrm{threshold} \times \left( \frac{V}{N} \right)^{1/3},
$$

where $V$ is the volume of the cell, $N$ the number of atoms and `threshold` a parameter (default 2).
If an atom has too few candidates, `find_neighbors` warns and asks for a larger `threshold`.

## Choosing a cutoff from g(r)

The radial distribution function shows where the neighbor shells are.
We compute it for three MD snapshots from the `examples` folder: an fcc crystal, a bcc crystal and a liquid.

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.io import read

structures = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

rdf = {}
first_minimum = {}
for name, atoms in structures.items():
    # radial_distribution_function replaces the neighbor list, so use a copy
    g, r = pyscal.radial_distribution_function(atoms.copy(), rmax=6.0, bins=120)
    r = r + 0.5 * (r[1] - r[0])          # bin centres
    peak = np.argmax(g)
    minimum = peak + np.argmin(g[peak:peak + 40])
    rdf[name] = (r, g)
    first_minimum[name] = r[minimum]

first_minimum
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=3, ratio=0.33, wspace=0.12)
for k, name in enumerate(rdf):
    ax = mp[0, k]
    r, g = rdf[name]
    ax.plot(r, g, color=COLOURS[name], lw=1.8)
    ax.axvline(first_minimum[name], ls="--", color=DARK, lw=1)
    ax.axhline(1, ls=":", color=DARK, lw=0.8)
    label(ax, f"({'abc'[k]})  {name}")
    note(ax, f"first minimum\n{first_minimum[name]:.2f} Å", loc="upper right")
    ax.set_xlabel("$r$  (Å)")
    ax.set_xlim(2, 6)
    ax.set_ylim(0, 6.5)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("$g(r)$");
```

The dashed lines mark the first minimum of $g(r)$.
In the crystals, $g(r)$ drops almost to zero between the first and second shell, so a fixed cutoff there is well defined.
In bcc, the first minimum lies after the second shell, which is only 15 % farther away than the first, so the cutoff takes 8 + 6 = 14 neighbors.
In the liquid, the minimum is shallow, and the number of neighbors depends on where exactly the cutoff is placed.

## How the methods compare

We now find the neighbors of each structure with four methods and count the neighbors of each atom.
The fixed cutoff is the first minimum of $g(r)$ found above.

```{code-cell} ipython3
methods = {
    "fixed cutoff": lambda name: dict(method="cutoff", cutoff=first_minimum[name]),
    "adaptive": lambda name: dict(method="cutoff", cutoff=0),
    "SANN": lambda name: dict(method="cutoff", cutoff="sann"),
    "Voronoi": lambda name: dict(method="voronoi"),
}

counts = {}
for name, atoms in structures.items():
    for method, options in methods.items():
        pyscal.find_neighbors(atoms, **options(name))
        counts[name, method] = np.diff(atoms.info["pyscal_bond_offsets"])
```

`atoms.info["pyscal_bond_offsets"]` marks where the neighbors of each atom start in the flat neighbor arrays (see [How the neighbors are stored](#how-the-neighbors-are-stored)), so its differences are the numbers of neighbors.

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=3, ratio=0.36, wspace=0.12)
for k, name in enumerate(structures):
    ax = mp[0, k]
    for method, (colour, marker) in METHODS.items():
        n = counts[name, method]
        values, frequency = np.unique(n, return_counts=True)
        ax.plot(values, frequency / len(n), marker=marker, color=colour, mfc=colour,
                mec=DARK, mew=0.7, ms=5, lw=1.5, label=method)
    label(ax, f"({'abc'[k]})  {name}")
    ax.set_xlabel("Neighbors per atom")
    ax.set_xlim(4, 24)
    ax.set_ylim(0, 1.05)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("Fraction of atoms")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [hide-input]
import pandas as pd

pd.DataFrame({method: [counts[name, method].mean() for name in structures]
              for method in methods}, index=list(structures)).round(1)
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

def percent(name, method, n):
    return round(100 * np.mean(counts[name, method] == n))

glue("fixed_fcc_12", percent("fcc", "fixed cutoff", 12), display=False)
glue("fixed_bcc_14", percent("bcc", "fixed cutoff", 14), display=False)
glue("sann_fcc_12", percent("fcc", "SANN", 12), display=False)
glue("sann_bcc_14", percent("bcc", "SANN", 14), display=False)
```

The table gives the mean number of neighbors per atom.

- **Fixed cutoff.** With the cutoff at the first minimum of $g(r)$, {glue}`fixed_fcc_12` % of the fcc atoms have 12 neighbors and {glue}`fixed_bcc_14` % of the bcc atoms have 14, as in the perfect crystals. The others are atoms whose thermal displacement moved a neighbor across the cutoff.
- **Adaptive cutoff.** With the defaults, the adaptive cutoff takes the first shell in fcc but only part of the second shell in bcc, so the bcc atoms have between 6 and 14 neighbors. Larger values of `nlimit` or `padding` include the full second shell. In the liquid, it gives fewer neighbors than the other methods.
- **SANN.** Without any parameter, SANN gives 12 neighbors to {glue}`sann_fcc_12` % of the fcc atoms. In bcc it often misses one or two atoms of the second shell, and only {glue}`sann_bcc_14` % of the atoms have 14 neighbors.
- **Voronoi.** Thermal motion creates small extra faces of the Voronoi cells, so the Voronoi method finds more neighbors than the others, about 14 in fcc. The face-area weights make these extra neighbors count little in the Steinhardt parameters.

## Which method to use

- **A single crystal structure, or several with known lattice constants.** Use a fixed cutoff at the first minimum of $g(r)$. It is the fastest method and the easiest to reproduce.
- **Several phases, interfaces or large density changes.** Use SANN or the adaptive cutoff. Check the number of neighbors per atom, as above, before trusting the result.
- **No parameter at all, or face-area weights.** Use Voronoi neighbors. They are also required by `voronoi_vector`.
- **A fixed number of neighbors.** Use `method="number"`.

Some functions choose their own neighbors.
`common_neighbor_analysis` and `diamond_structure` search for their own neighbors and leave the stored list unchanged.
`radial_distribution_function`, `centrosymmetry` and `minkowski_parameter` replace the stored list with their own.
Call them on a copy (`atoms.copy()`), or before `find_neighbors`, if the list is needed afterwards.

## How the neighbors are stored

`find_neighbors` stores the neighbors on the `Atoms` object in two forms.

**Flat arrays**, in `atoms.info`.
`pyscal_bond_offsets` has $N + 1$ entries, and the neighbors of atom $i$ are the entries `offsets[i]` to `offsets[i + 1]` of the other arrays:

| Key | Content |
|---|---|
| `pyscal_bond_neighbors` | index $j$ of each neighbor |
| `pyscal_bond_distance` | distance $r_{ij}$ |
| `pyscal_bond_vector` | vector $\mathbf{r}_i - \mathbf{r}_j$, one row per neighbor |
| `pyscal_bond_weight` | weight $w_{ij}$ (Voronoi face weights, 1 for the other methods) |
| `pyscal_bond_theta`, `pyscal_bond_phi` | polar and azimuthal angle of $\mathbf{r}_i - \mathbf{r}_j$ |

The adaptive, SANN and number methods also store their candidates as `pyscal_candidate_offsets`, `pyscal_candidate_neighbors` and `pyscal_candidate_distance`.

```{code-cell} ipython3
atoms = structures["fcc"]
pyscal.find_neighbors(atoms, method="cutoff", cutoff=first_minimum["fcc"])

offsets = atoms.info["pyscal_bond_offsets"]
neighbors = atoms.info["pyscal_bond_neighbors"]
distances = atoms.info["pyscal_bond_distance"]

i = 0
neighbors[offsets[i]:offsets[i + 1]], distances[offsets[i]:offsets[i + 1]].round(2)
```

**Rows, one per atom.**
`atoms.arrays["pyscal_neighbors"]` (or `atoms.info["pyscal_neighbors"]`) holds the neighbors of each atom, and `pyscal_neighbordist`, `pyscal_neighborweight`, `pyscal_r`, `pyscal_theta`, `pyscal_phi` and `pyscal_diff` the values per neighbor.
When every atom has the same number $k$ of neighbors, these are $(N, k)$ arrays in `atoms.arrays`.
Otherwise they are lists of lists in `atoms.info`.

Building these lists is the slowest part of the search when the numbers of neighbors differ.
All descriptors read the flat arrays, so the rows can be skipped:

```{code-cell} ipython3
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0, store_rows=False)
```

With `store_rows=False`, `atoms` can also be written to an extended XYZ file with `atoms.write()`, which fails on lists of lists in `atoms.info`.

## Periodic boundaries

`find_neighbors` uses the cell and the periodic boundary conditions of `atoms` (`atoms.cell`, `atoms.pbc`).

- Directions with `pbc=False` have no periodic images.
- A structure without a cell, such as a molecule read from a plain XYZ file, is treated as isolated.
- When the cutoff is larger than half the width of the cell, an atom can be a neighbor several times, through different periodic images. Each image is a separate entry with its own distance and vector.

The search runs on all CPU cores. To use fewer, see [Number of threads](../install.md#number-of-threads).

## References

1. A. Stukowski, Structure identification methods for atomistic simulations of crystalline materials, *Modelling Simul. Mater. Sci. Eng.* **20**, 045021 (2012). [doi:10.1088/0965-0393/20/4/045021](https://doi.org/10.1088/0965-0393/20/4/045021)
2. J. A. van Meel, L. Filion, C. Valeriani and D. Frenkel, A parameter-free, solid-angle based, nearest-neighbor algorithm, *J. Chem. Phys.* **136**, 234107 (2012). [doi:10.1063/1.4729313](https://doi.org/10.1063/1.4729313)
