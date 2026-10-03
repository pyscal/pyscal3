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

# Steinhardt parameters

The bond orientational order parameters of Steinhardt, Nelson and Ronchetti [1] describe the arrangement of the neighbors of an atom with spherical harmonics.
They do not depend on the orientation of the crystal, and different crystal structures have different values.
They are used to tell crystal structures apart, to separate solid from liquid atoms, and as input to [solid–liquid classification](solid_liquid), [disorder parameters](../methods/04_disorder) and [Wigner $W_l$ parameters](../methods/10_wigner_w).

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

For an atom $i$ with $N(i)$ neighbors, the complex vector $q_{lm}(i)$ is the average of the spherical harmonics $Y_{lm}$ over the bonds to the neighbors:

$$
q_{lm}(i) = \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\mathbf{r}_{ij}), \qquad m = -l, \dots, l,
$$

where $\mathbf{r}_{ij}$ is the vector from atom $i$ to its neighbor $j$.
The Steinhardt parameter $q_l$ is the norm of this vector:

$$
q_l(i) = \left( \frac{4\pi}{2l + 1} \sum_{m=-l}^{l} \left| q_{lm}(i) \right|^2 \right)^{1/2}.
$$

The norm does not change when the crystal is rotated.
$q_l$ lies between 0 and 1.
The parameters with $l = 4$ and $l = 6$ are used most often.

With [Voronoi neighbors](../guide/neighbors), each term in the sum over $j$ is multiplied by the weight $w_{ij}$ of the neighbor, the relative area of the Voronoi face shared with $j$.
`minkowski_parameter` computes these weighted parameters directly (see [Minkowski structure metrics](../methods/09_minkowski)).

### Averaged parameters

At finite temperature, thermal vibrations broaden the distributions of $q_l$, and the distributions of different structures overlap.
Lechner and Dellago [2] proposed to average $q_{lm}$ over the atom and its neighbors before taking the norm:

$$
\bar{q}_{lm}(i) = \frac{1}{N(i) + 1} \sum_{k=0}^{N(i)} q_{lm}(k), \qquad
\bar{q}_l(i) = \left( \frac{4\pi}{2l + 1} \sum_{m=-l}^{l} \left| \bar{q}_{lm}(i) \right|^2 \right)^{1/2},
$$

where $k = 0$ is the atom $i$ itself and $k = 1, \dots, N(i)$ are its neighbors.
The averaged parameters $\bar{q}_l$ include information from the second neighbor shell, and their distributions are much narrower.

## Usage

Find the neighbors first, then compute the parameters for one or several values of $l$:

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)

q4, q6 = pyscal.steinhardt_parameter(atoms, l=[4, 6])
q4_averaged, q6_averaged = pyscal.steinhardt_parameter(atoms, l=[4, 6], averaged=True)
```

`steinhardt_parameter` returns one array per value of $l$, each with one value per atom.
It also stores the results on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_q4"]` | $(N,)$ | $q_4$ |
| `atoms.arrays["pyscal_q4_real"]`, `["pyscal_q4_imag"]` | $(N, 2l + 1)$ | real and imaginary parts of $q_{4m}$, $m = -l, \dots, l$ |
| `atoms.arrays["pyscal_avg_q4"]` | $(N,)$ | $\bar{q}_4$, with `averaged=True` |

The $q_{lm}$ are stored because the averaged parameters, [solid–liquid classification](solid_liquid) and [disorder parameters](../methods/04_disorder) need them.
With `averaged=True`, the $q_{lm}$ are computed on the way, so the plain parameters do not have to be computed first.

## Values for perfect crystals

```{code-cell} ipython3
import numpy as np
import pandas as pd
from pyscal.structures import make_crystal

crystals = {
    "fcc": make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4)),
    "bcc": make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4)),
    "hcp": make_crystal("hcp", lattice_constant=2.51, repetitions=(4, 4, 4)),
}

ls = list(range(2, 13))
perfect = {}
for name, crystal in crystals.items():
    pyscal.find_neighbors(crystal, method="cutoff", cutoff=0)
    perfect[name] = [q.mean() for q in pyscal.steinhardt_parameter(crystal, l=ls)]

pd.DataFrame(perfect, index=[f"q{l}" for l in ls]).T.round(3)
```

In a perfect crystal every atom has the same value.
The neighbors are the first shell: 12 atoms in fcc and hcp, 14 (8 + 6) in bcc.
In fcc and bcc, every neighbor at $\mathbf{r}$ has a partner at $-\mathbf{r}$, and $q_l = 0$ for all odd $l$.
In hcp the neighbors do not come in such pairs, and the odd $l$ are not zero.

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(ratio=0.42)
ax = mp[0, 0]
for (name, values), marker in zip(perfect.items(), "osD"):
    ax.plot(ls, values, marker=marker, color=COLOURS[name], mec=DARK, mew=0.7, ms=6,
            lw=1.5, label=name)
ax.set_xticks(ls)
ax.set_xlabel("$l$")
ax.set_ylabel("$q_l$")
ax.set_ylim(0, 0.75)
ax.legend(frameon=False, ncol=3, loc="upper left");
```

$q_4$ is clearly different for the three structures, and $q_6$ is about 0.5 for all of them and lower in liquids.
This is why the $(q_4, q_6)$ plane is the usual map for identifying structures.

## Crystals at finite temperature

The figure below shows the $(q_4, q_6)$ plane for three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.
The stars mark the values of the perfect crystals from the table above.

```{code-cell} ipython3
snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

maps = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff="sann")
    plain = pyscal.steinhardt_parameter(snapshot, l=[4, 6])
    averaged = pyscal.steinhardt_parameter(snapshot, l=[4, 6], averaged=True)
    maps[name] = (plain, averaged)
```

The neighbors are found with SANN, which gives the full first shell in both crystals (see [Finding neighbors](../guide/neighbors)).

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.5, wspace=0.12)
for k, title in enumerate(["(a)  $q_l$", r"(b)  $\bar{q}_l$"]):
    ax = mp[0, k]
    for name, values in maps.items():
        q4, q6 = values[k]
        ax.scatter(q4, q6, s=6, color=COLOURS[name], alpha=0.5, lw=0, label=name,
                   zorder=2 if name == "liquid" else 3)
    for name in ("fcc", "bcc", "hcp"):
        q4, q6 = perfect[name][2], perfect[name][4]
        ax.plot(q4, q6, marker="*", ms=14, color=COLOURS[name], mec=DARK, mew=0.8,
                ls="none", zorder=4)
        ax.annotate(name, (q4, q6), xytext=(6, 6), textcoords="offset points",
                    fontsize=9, color=DARK)
    label(ax, title)
    ax.set_xlim(0, 0.45)
    ax.set_ylim(0, 0.7)
    ax.set_xlabel("$q_4$" if k == 0 else r"$\bar{q}_4$")
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("$q_6$  or  $\\bar{q}_6$")
legend = mp[0, 1].legend(frameon=False, loc="lower right", markerscale=3)
for handle in legend.legend_handles:
    handle.set_alpha(1)
```

Panel (a) shows that the clouds of the fcc crystal, the bcc crystal and the liquid overlap in the plane of the plain parameters.
In panel (b), the averaged parameters form three separate clouds.
The clouds of the crystals lie below the values of the perfect crystals, because thermal motion lowers the order, but they no longer overlap.
For classifying atoms at finite temperature, use the averaged parameters.

The averaging has a cost: the value of an atom depends on its neighbors, so the resolution near an interface or a defect is about one neighbor shell coarser.

## Things to watch

- **The neighbor method changes the values.** $q_l$ depends on which atoms are counted as neighbors. Compare values only when they were computed with the same neighbor method. The values in the table above use the first shell.
- **bcc needs the second shell.** In bcc, the second shell is only 15 % farther away than the first. With only the first 8 neighbors, a perfect bcc crystal has $q_4 = 0.51$ and $q_6 = 0.63$ instead of 0.04 and 0.51.
- **Neighbor lists that some functions replace.** `minkowski_parameter`, `centrosymmetry` and `radial_distribution_function` replace the stored neighbor list. Compute the Steinhardt parameters before calling them, or call them on a copy.

## References

1. P. J. Steinhardt, D. R. Nelson and M. Ronchetti, Bond-orientational order in liquids and glasses, *Phys. Rev. B* **28**, 784 (1983). [doi:10.1103/PhysRevB.28.784](https://doi.org/10.1103/PhysRevB.28.784)
2. W. Lechner and C. Dellago, Accurate determination of crystal structures based on averaged local bond order parameters, *J. Chem. Phys.* **129**, 114707 (2008). [doi:10.1063/1.2977970](https://doi.org/10.1063/1.2977970)
