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

# Angular criteria and χ parameters

The two descriptors on this page are computed from the angles between the bonds of an atom to its neighbors.
The angular criterion $A$ of Uttormark et al. [1] measures how far the four nearest neighbors of an atom are from a regular tetrahedron.
It is used to find atoms in diamond structures, such as crystalline silicon or germanium.
The χ parameters of Ackland and Jones [2] are a histogram of the cosines of all bond angles of an atom.
They are the input of the [Ackland–Jones classification](ackland_jones).
For a structure label of diamond atoms, see [Common neighbor analysis](cna), and for the bond angle distribution of a whole system, see [Distribution functions](distributions).

```{code-cell} ipython3
:tags: [remove-cell]
import os, sys, warnings
sys.path.insert(0, os.path.abspath(".."))
from _plotstyle import EXAMPLES, COLOURS, DARK, figure, label, note
os.chdir(EXAMPLES)
warnings.simplefilter("ignore")
%config InlineBackend.figure_format = "retina"
```

## Definition

### Angular criterion

For an atom $i$, take its four nearest neighbors.
The bonds to these neighbors form six angles $\theta_{jk}$, one for each pair $j < k$ of the four neighbors.
The angular criterion is

$$
A(i) = \sum_{j < k} \left( \cos\theta_{jk} + \frac{1}{3} \right)^2,
$$

where $\cos\theta_{jk} = \mathbf{r}_{ij} \cdot \mathbf{r}_{ik} / (r_{ij} r_{ik})$, and $\mathbf{r}_{ij}$ is the vector from atom $i$ to its neighbor $j$.
The tetrahedral angle of 109.47° has $\cos\theta = -1/3$.
$A = 0$ when the four neighbors form a regular tetrahedron around the atom, and $A > 0$ otherwise.

### χ parameters

For an atom $i$ with $N(i)$ neighbors, compute $\cos\theta_{jk}$ for all $N(i)\,(N(i) - 1)/2$ pairs of neighbors.
The χ parameter $\chi_n(i)$ is the number of pairs whose cosine falls in bin $n$:

| Parameter | $\cos\theta$ | $\theta$ |
|---|---|---|
| $\chi_0$ | $[-1, -0.945)$ | 160.9° to 180° |
| $\chi_1$ | $[-0.945, -0.915)$ | 156.2° to 160.9° |
| $\chi_2$ | $[-0.915, -0.755)$ | 139.0° to 156.2° |
| $\chi_3$ | $[-0.755, -0.705)$ | 134.8° to 139.0° |
| $\chi_4$ | $[-0.705, -0.195)$ | 101.2° to 134.8° |
| $\chi_5$ | $[-0.195, 0.195)$ | 78.8° to 101.2° |
| $\chi_6$ | $[0.195, 0.245)$ | 75.8° to 78.8° |
| $\chi_7$ | $[0.245, 0.795)$ | 37.3° to 75.8° |
| $\chi_8$ | $[0.795, 1]$ | 0° to 37.3° |

Ackland and Jones [2] use eight bins, $\chi_0$ to $\chi_7$.
pyscal splits their fourth bin, $[-0.755, -0.195)$, at $-0.705$ into $\chi_3$ and $\chi_4$.
The other bins are the same, so $\chi_5$ to $\chi_8$ in pyscal are $\chi_4$ to $\chi_7$ in the paper.

## Usage

Find the neighbors first, then compute the descriptors:

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.build import bulk

silicon = bulk("Si", "diamond", a=5.43, cubic=True).repeat(4)
pyscal.find_neighbors(silicon, method="cutoff", cutoff=0)

A = pyscal.angular_criteria(silicon)
chi = pyscal.chi_params(silicon)
A[:4].round(3), chi[0]
```

`angular_criteria` returns $A$ for each atom, and `chi_params` returns an $(N, 9)$ array of integer counts.
Both store their results on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_angular"]` | $(N,)$ | $A$ |
| `atoms.arrays["pyscal_chiparams"]` | $(N, 9)$ | $\chi_0, \dots, \chi_8$ |
| `atoms.info["pyscal_cosines"]` | list of $N$ lists | $\cos\theta_{jk}$ for all pairs of neighbors, with `angles=True` |

With `angles=True`, `chi_params` returns the χ parameters and the list of cosines.

Both functions use the stored neighbor list and raise an error if `find_neighbors` has not been called.
`angular_criteria` uses only the four nearest of the stored neighbors, so any neighbor method that finds at least four neighbors gives the same result.
`chi_params` uses all stored neighbors, and the counts change with the neighbor method.

## Values for perfect crystals

```{code-cell} ipython3
import pandas as pd
from ase.cluster import Icosahedron
from pyscal.structures import make_crystal

icosahedron = Icosahedron("Cu", noshells=2)   # 13 atoms, the first is the centre
icosahedron.center(vacuum=5.0)

crystals = {
    "fcc": make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4)),
    "hcp": make_crystal("hcp", lattice_constant=2.51, repetitions=(4, 4, 4)),
    "bcc": make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4)),
    "ico": icosahedron,
    "diamond": make_crystal("diamond", lattice_constant=5.43, repetitions=(3, 3, 3)),
}

perfect = {}
for name, crystal in crystals.items():
    pyscal.find_neighbors(crystal, method="cutoff", cutoff=0)
    perfect[name] = pyscal.chi_params(crystal)[0]

table = pd.DataFrame(perfect, index=[f"χ{n}" for n in range(9)]).T
table["neighbors"] = table.sum(axis=1).map(lambda pairs: int((1 + np.sqrt(1 + 8 * pairs)) / 2))
table
```

The table gives the χ parameters of one atom in each perfect structure, with first shell neighbors: 12 in fcc, hcp and the icosahedron, 14 (8 + 6) in bcc and 4 in diamond.
For the icosahedron, the atom is the centre of a 13-atom cluster.
$\chi_0$ counts pairs of neighbors on nearly opposite sides of the atom.
There are 7 such pairs in bcc, 6 in fcc and the icosahedron, and 3 in hcp.
$\chi_2$ is nonzero only in hcp, from the six angles of 146.4° between neighbors in the layers above and below the atom.
$\chi_5$ counts angles near 90°, which are absent in the icosahedron.
These differences are the basis of the [Ackland–Jones classification](ackland_jones).

## Crystals at finite temperature

Thermal vibrations spread the angles, and pairs move into neighboring bins.
We compute the χ parameters for three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.
The neighbors are all atoms within the first minimum of the radial distribution function (see [Finding neighbors](../guide/neighbors)).

```{code-cell} ipython3
from ase.io import read

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}
cutoffs = {"fcc": 3.5, "bcc": 3.8, "liquid": 4.1}

chi_snapshots = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=cutoffs[name])
    chi_snapshots[name] = pyscal.chi_params(snapshot)
```

```{code-cell} ipython3
:tags: [hide-input]
bins = np.arange(9)
mp = figure(columns=2, ratio=0.42, wspace=0.12)
groups = [["fcc", "hcp", "bcc", "ico"], ["fcc", "bcc", "liquid"]]
for k, names in enumerate(groups):
    ax = mp[0, k]
    width = 0.8 / len(names)
    for m, name in enumerate(names):
        x = bins + (m - (len(names) - 1) / 2) * width
        if k == 0:
            ax.bar(x, perfect[name], width, color=COLOURS[name], ec=DARK, lw=0.6,
                   label=name)
        else:
            chi = chi_snapshots[name]
            ax.bar(x, chi.mean(axis=0), width, color=COLOURS[name], ec=DARK, lw=0.6,
                   yerr=chi.std(axis=0), ecolor=DARK,
                   error_kw=dict(elinewidth=0.8, capsize=1.5), label=name)
    ax.set_xticks(bins, [f"$\\chi_{n}$" for n in bins])
    ax.set_xlim(-0.6, 8.6)
    ax.set_ylim(0, 46)
    ax.set_xlabel("Bin")
    if k:
        ax.tick_params(labelleft=False)
label(mp[0, 0], "(a)  perfect structures")
label(mp[0, 1], "(b)  MD snapshots")
mp[0, 0].set_ylabel("Number of pairs")
handles, labels = [], []
for ax in (mp[0, 0], mp[0, 1]):
    for h, l in zip(*ax.get_legend_handles_labels()):
        if l not in labels:
            handles.append(h)
            labels.append(l)
mp.fig.legend(handles, labels, frameon=False, ncol=5, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

def percent(values):
    return round(100 * float(np.mean(values)))

glue("fcc_chi0_mean", round(float(chi_snapshots["fcc"][:, 0].mean()), 1), display=False)
glue("fcc_chi123_percent", percent(chi_snapshots["fcc"][:, 1:4].sum(axis=1) > 0), display=False)
glue("bcc_chi0_mean", round(float(chi_snapshots["bcc"][:, 0].mean()), 1), display=False)
glue("liquid_chi0_mean", round(float(chi_snapshots["liquid"][:, 0].mean()), 1), display=False)
glue("liquid_chi2_mean", round(float(chi_snapshots["liquid"][:, 2].mean()), 1), display=False)
```

Panel (a) shows the values of the table.
Panel (b) shows the mean over the atoms of each snapshot, and the error bars give the standard deviation over the atoms.
In the fcc snapshot, $\chi_0$ drops from 6 to {glue}`fcc_chi0_mean` on average, and {glue}`fcc_chi123_percent` % of the atoms have at least one angle in $\chi_1$ to $\chi_3$, which are empty in the perfect crystal.
In the bcc snapshot, $\chi_0$ drops from 7 to {glue}`bcc_chi0_mean`.
In the liquid, $\chi_0$ is {glue}`liquid_chi0_mean` on average, and many pairs fall in bins that are empty in fcc and bcc, for example {glue}`liquid_chi2_mean` in $\chi_2$.
Classifiers that test for exact counts, such as the [Ackland–Jones classification](ackland_jones), are therefore sensitive to temperature.

## Finding tetrahedral atoms

Does $A$ separate atoms in a diamond structure from atoms in other structures at finite temperature?
We add random displacements to perfect cubic and hexagonal diamond, with a standard deviation of 4 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations.
We compare them with the three MD snapshots.

```{code-cell} ipython3
a = 5.43
nearest = a * np.sqrt(3) / 4
diamonds = {
    "diamond": bulk("Si", "diamond", a=a, cubic=True).repeat(4),
    "hex. diamond": bulk("SiSi", "wurtzite", a=a / np.sqrt(2),
                         c=a / np.sqrt(2) * np.sqrt(8 / 3)).repeat((6, 6, 4)),
}

angular = {}
for name, structure in diamonds.items():
    structure.rattle(0.04 * nearest, seed=1)
    pyscal.find_neighbors(structure, method="cutoff", cutoff=0)
    angular[name] = pyscal.angular_criteria(structure)

for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=0)
    angular[name] = pyscal.angular_criteria(snapshot)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(ratio=0.42)
ax = mp[0, 0]
edges = np.logspace(-3, 1, 61)
for name, values in angular.items():
    filled = "diamond" in name
    ax.hist(values, bins=edges, density=False, weights=np.full(len(values), 1 / len(values)),
            histtype="stepfilled" if filled else "step", color=COLOURS[name],
            alpha=0.75 if filled else 1, ec=DARK if filled else COLOURS[name],
            lw=0.6 if filled else 1.8, label=name)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-3, 10)
ax.set_ylim(5e-4, 0.6)
ax.set_xlabel("$A$")
ax.set_ylabel("Fraction of atoms")
ax.legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.22));
```

```{code-cell} ipython3
:tags: [remove-cell]
largest = max(angular["diamond"].max(), angular["hex. diamond"].max())
glue("diamond_max", round(float(largest), 2), display=False)
for name in ("fcc", "bcc", "liquid"):
    glue(f"{name}_below", int(np.sum(angular[name] < largest)), display=False)
    glue(f"{name}_total", len(angular[name]), display=False)
```

Both diamond structures have $A$ below {glue}`diamond_max`.
Only {glue}`fcc_below` of the {glue}`fcc_total` atoms of the fcc snapshot and {glue}`liquid_below` of the {glue}`liquid_total` atoms of the liquid have smaller values.
In the bcc snapshot, {glue}`bcc_below` of the {glue}`bcc_total` atoms are below this value, and a few have $A$ as small as in the diamond structures.
The eight nearest neighbors in bcc sit at the corners of a cube, and four alternating corners form a regular tetrahedron.
When thermal motion brings these four atoms closer than the other four, $A$ is close to zero.
A small value of $A$ therefore shows that the four nearest neighbors are arranged as a tetrahedron, but not that the atom is in a diamond structure.
Cubic and hexagonal diamond have the same values, so $A$ does not tell them apart.
Use `diamond_structure` (see [Common neighbor analysis](cna.md#diamond-structures)) to identify and separate the two.

## Things to watch

- **Fewer than four neighbors.** An atom with fewer than four stored neighbors gets $A = 0$, the value of a perfect tetrahedron. Check the number of neighbors, for example with `np.diff(atoms.info["pyscal_bond_offsets"])`, before selecting atoms by small $A$.
- **Ties between neighbors.** In close packed crystals, the 12 nearest neighbors are at the same distance, and the four nearest are chosen among them by index. $A$ then takes several values in a perfect crystal. $A$ is meaningful only for atoms with four clearly nearest neighbors.
- **χ depends on the number of neighbors.** The χ parameters are counts, and they add up to $N(i)\,(N(i) - 1)/2$. A neighbor method that gives one extra neighbor adds $N(i)$ pairs. Compare χ parameters only when they were computed with the same neighbor method.

## References

1. M. J. Uttormark, M. O. Thompson and P. Clancy, Kinetics of crystal dissolution for a Stillinger–Weber model of silicon, *Phys. Rev. B* **47**, 15717 (1993). [doi:10.1103/PhysRevB.47.15717](https://doi.org/10.1103/PhysRevB.47.15717)
2. G. J. Ackland and A. P. Jones, Applications of local crystal structure measures in experiment and simulation, *Phys. Rev. B* **73**, 054104 (2006). [doi:10.1103/PhysRevB.73.054104](https://doi.org/10.1103/PhysRevB.73.054104)
