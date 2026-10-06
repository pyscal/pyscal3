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

# Minkowski structure metrics

The Minkowski structure metrics of Mickel et al. [1] are [Steinhardt parameters](steinhardt) in which each neighbor is weighted by the area of the face it shares with the atom in the [Voronoi tessellation](voronoi).
They need no cutoff or other parameter, and they change continuously when the atoms move.
Plain Steinhardt parameters jump whenever an atom enters or leaves the neighbor list.

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

The neighbors of an atom $i$ are the atoms whose Voronoi cells share a face with the cell of $i$.
Each neighbor $j$ gets a weight from the area $A_{ij}$ of the shared face:

$$
w_{ij} = \frac{A_{ij}^{\,p}}{\sum_{k} A_{ik}^{\,p}},
$$

where the sum runs over all Voronoi neighbors $k$ of $i$, and the exponent $p$ is set by `voroexp` (default 1).
The weights replace the factor $1/N(i)$ in the definition of $q_{lm}$:

$$
q_{lm}(i) = \sum_{j} w_{ij}\, Y_{lm}(\mathbf{r}_{ij}), \qquad
q_l(i) = \left( \frac{4\pi}{2l + 1} \sum_{m=-l}^{l} \left| q_{lm}(i) \right|^2 \right)^{1/2},
$$

where $Y_{lm}$ are the spherical harmonics and $\mathbf{r}_{ij}$ is the vector from atom $i$ to its neighbor $j$.
With $p = 1$, $w_{ij}$ is the fraction of the surface of the Voronoi cell of $i$ that it shares with $j$, as in the definition of Mickel et al. [1].
With $p = 0$, every Voronoi neighbor has the same weight, as in the plain Steinhardt parameters.

When an atom moves, the faces of the Voronoi cell change continuously.
A new neighbor enters the cell with a face of zero area, and a neighbor leaves when its face shrinks to zero.
With $p > 0$, $q_l$ is therefore a continuous function of the positions.

The averaged parameters $\bar{q}_l$ of Lechner and Dellago [2] are computed as on the [Steinhardt parameters](steinhardt) page, with the Voronoi neighbors.
The average over the neighbors is not weighted.
Every Voronoi neighbor counts once, whatever the area of its face.

## Usage

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")

q4, q6 = pyscal.minkowski_parameter(atoms, l=[4, 6])
q4_averaged, q6_averaged = pyscal.minkowski_parameter(atoms, l=[4, 6], averaged=True)
```

`minkowski_parameter` finds the Voronoi neighbors itself, with `find_neighbors(atoms, method="voronoi", voroexp=voroexp)`, and then calls `steinhardt_parameter`.
`find_neighbors` does not have to be called first, and a neighbor list stored on `atoms` is replaced by the Voronoi neighbors.

`minkowski_parameter` returns one array per value of $l$, each with one value per atom.
It stores the same keys as `steinhardt_parameter`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_q4"]` | $(N,)$ | $q_4$ |
| `atoms.arrays["pyscal_q4_real"]`, `["pyscal_q4_imag"]` | $(N, 2l + 1)$ | real and imaginary parts of $q_{4m}$, $m = -l, \dots, l$ |
| `atoms.arrays["pyscal_avg_q4"]` | $(N,)$ | $\bar{q}_4$, with `averaged=True` |

The Voronoi neighbors, their weights and the Voronoi volumes are also stored, as described on the [Voronoi tessellation](voronoi) page.

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
table = {}
for name, crystal in crystals.items():
    pyscal.find_neighbors(crystal, method="cutoff", cutoff=0)
    q4, q6 = pyscal.steinhardt_parameter(crystal, l=[4, 6])
    q4_minkowski, q6_minkowski = pyscal.minkowski_parameter(crystal, l=[4, 6])
    table[name] = [q4.mean(), q6.mean(), q4_minkowski.mean(), q6_minkowski.mean()]

perfect = pd.DataFrame(table, index=["q4", "q6", "q4 Minkowski", "q6 Minkowski"]).T
perfect.round(4)
```

The first two columns are the plain Steinhardt parameters with the first neighbor shell.
In fcc and hcp, all 12 faces of the Voronoi cell have the same area, and the two definitions give the same values.
In bcc, the Voronoi cell has 8 hexagonal faces shared with the first shell and 6 square faces shared with the second shell:

```{code-cell} ipython3
np.unique(crystals["bcc"].info["pyscal_bond_weight"].round(4), return_counts=True)
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
weights = np.unique(crystals["bcc"].info["pyscal_bond_weight"].round(4))
glue("ratio_hexagon_square", round(weights[1] / weights[0], 1), display=False)
glue("bcc_q4_plain", round(perfect.loc["bcc", "q4"], 3), display=False)
glue("bcc_q4_minkowski", round(perfect.loc["bcc", "q4 Minkowski"], 3), display=False)
glue("fcc_q4_minkowski", round(perfect.loc["fcc", "q4 Minkowski"], 3), display=False)
```

A hexagonal face has {glue}`ratio_hexagon_square` times the area of a square face, so the second shell counts less than in the plain parameters.
$q_4$ of bcc rises from {glue}`bcc_q4_plain` to {glue}`bcc_q4_minkowski`, close to the value of fcc, {glue}`fcc_q4_minkowski`.

## Continuity under deformation

The Bain path turns a bcc crystal into an fcc crystal by stretching it along one cube axis.
We follow it at constant volume, from the ratio $c/a = 1$ of the cubic bcc cell to $c/a = \sqrt{2}$, where the body centred tetragonal cell describes fcc.
On the way, two of the 14 neighbors of bcc move away from the atom, and the other 12 become the first shell of fcc.
We compare the Minkowski metrics with plain Steinhardt parameters for a fixed cutoff of 3 Å and for SANN neighbors.

```{code-cell} ipython3
from ase import Atoms

volume = 2.87**3 / 2           # volume per atom of the bcc crystal, in Å^3
ratios = np.linspace(1, np.sqrt(2), 85)

def bain(ratio):
    a = (2 * volume / ratio) ** (1 / 3)
    cell = [a, a, ratio * a]
    return Atoms("Fe2", scaled_positions=[[0, 0, 0], [0.5, 0.5, 0.5]], cell=cell,
                 pbc=True).repeat(4)

methods = {
    "fixed cutoff": lambda atoms: pyscal.find_neighbors(atoms, method="cutoff", cutoff=3.0),
    "SANN": lambda atoms: pyscal.find_neighbors(atoms, method="cutoff", cutoff="sann"),
}
path = {method: [] for method in [*methods, "Minkowski"]}
for ratio in ratios:
    crystal = bain(ratio)
    for method, find in methods.items():
        find(crystal)
        path[method].append([q.mean() for q in pyscal.steinhardt_parameter(crystal, l=[4, 6])])
    path["Minkowski"].append([q.mean() for q in pyscal.minkowski_parameter(crystal, l=[4, 6])])
path = {method: np.array(values) for method, values in path.items()}
```

```{code-cell} ipython3
:tags: [hide-input]
from _plotstyle import METHODS

styles = {"fixed cutoff": (*METHODS["fixed cutoff"], "--", (0, 8)),
          "SANN": (*METHODS["SANN"], ":", (4, 8)),
          "Minkowski": (*METHODS["Voronoi"], "-", (0, 8))}
mp = figure(columns=2, ratio=0.42, wspace=0.3)
for k, title in enumerate(["(a)  $q_4$", "(b)  $q_6$"]):
    ax = mp[0, k]
    for method, values in path.items():
        colour, marker, ls, every = styles[method]
        ax.plot(ratios, values[:, k], color=colour, marker=marker, markevery=every, ms=5,
                mec=DARK, mew=0.7, lw=1.6, ls=ls, label=method)
    label(ax, title)
    ax.set_xlabel("$c/a$")
    ax.set_xticks([1, 1.1, 1.2, 1.3, np.sqrt(2)], ["1\nbcc", "1.1", "1.2", "1.3", "$\\sqrt{2}$\nfcc"])
mp[0, 0].set_ylabel("$q_l$")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.12));
```

```{code-cell} ipython3
:tags: [remove-cell]
def largest_step(values):
    return round(float(np.abs(np.diff(values[:, 0])).max()), 3)

glue("step_cutoff", largest_step(path["fixed cutoff"]), display=False)
glue("step_sann", largest_step(path["SANN"]), display=False)
glue("step_minkowski", largest_step(path["Minkowski"]), display=False)
```

With the fixed cutoff and with SANN, $q_4$ and $q_6$ jump when the neighbor list changes from 14 to 12 atoms.
The largest change of $q_4$ between two consecutive points of the path is {glue}`step_cutoff` for the fixed cutoff and {glue}`step_sann` for SANN, but only {glue}`step_minkowski` for the Minkowski metrics.
The Minkowski metrics follow a smooth curve from the bcc to the fcc values.
This matters when a structure changes gradually, for example along a transformation path.

## Crystals at finite temperature

The figure below shows the $(q_4, q_6)$ plane of the Minkowski metrics for three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.
The stars mark the values of the perfect crystals from the table above.
The same plane for plain Steinhardt parameters is shown on the [Steinhardt parameters](steinhardt) page.

```{code-cell} ipython3
snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}
maps = {}
for name, snapshot in snapshots.items():
    plain = pyscal.minkowski_parameter(snapshot, l=[4, 6])
    averaged = pyscal.minkowski_parameter(snapshot, l=[4, 6], averaged=True)
    maps[name] = (plain, averaged)
```

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
        q4, q6 = perfect.loc[name, "q4 Minkowski"], perfect.loc[name, "q6 Minkowski"]
        ax.plot(q4, q6, marker="*", ms=14, color=COLOURS[name], mec=DARK, mew=0.8,
                ls="none", zorder=4)
        offset, ha = ((-8, 6), "right") if name == "fcc" else ((8, 6), "left")
        ax.annotate(name, (q4, q6), xytext=offset, textcoords="offset points", ha=ha,
                    fontsize=9, color=DARK)
    label(ax, title)
    ax.set_xlim(0, 0.5)
    ax.set_ylim(0, 0.7)
    ax.set_xlabel("$q_4$" if k == 0 else r"$\bar{q}_4$")
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("$q_6$  or  $\\bar{q}_6$")
legend = mp[0, 1].legend(frameon=False, loc="lower right", markerscale=3)
for handle in legend.legend_handles:
    handle.set_alpha(1)
```

```{code-cell} ipython3
:tags: [remove-cell]
(_, fcc_q6), (_, liquid_q6) = maps["fcc"][1], maps["liquid"][1]
glue("liquid_q6_max", round(float(liquid_q6.max()), 2), display=False)
glue("fcc_q6_min", round(float(fcc_q6.min()), 2), display=False)
glue("bcc_q6_min", round(float(maps["bcc"][1][1].min()), 2), display=False)
```

In panel (a), the clouds of the three structures overlap.
In panel (b), the averaged metrics separate the liquid from the crystals.
$\bar{q}_6$ is at most {glue:text}`liquid_q6_max:.2f` in the liquid and at least {glue:text}`fcc_q6_min:.2f` in fcc and {glue:text}`bcc_q6_min:.2f` in bcc.
The fcc and bcc clouds lie close together, because the face weights bring $q_4$ of bcc close to that of fcc (see the table above).
To tell fcc from bcc, use the averaged plain Steinhardt parameters with a neighbor method that includes the second shell of bcc, or the [Wigner parameters](wigner_w).

```{code-cell} ipython3
:tags: [remove-cell]
fcc = snapshots["fcc"]
q6_by_exponent = {e: pyscal.minkowski_parameter(fcc, l=6, voroexp=e)[0].mean() for e in (0, 1)}
glue("fcc_voronoi_neighbors", round(float(np.diff(fcc.info["pyscal_bond_offsets"]).mean()), 1),
     display=False)
glue("fcc_q6_exp0", round(float(q6_by_exponent[0]), 3), display=False)
glue("fcc_q6_exp1", round(float(q6_by_exponent[1]), 3), display=False)
glue("fcc_q6_perfect", round(float(perfect.loc["fcc", "q6 Minkowski"]), 3), display=False)
```

## Things to watch

- **The neighbor list is replaced.** `minkowski_parameter` replaces the stored neighbor list with Voronoi neighbors, and descriptors computed afterwards use them. Call it on a copy (`atoms.copy()`), or call `find_neighbors` again afterwards.
- **Values differ from plain Steinhardt parameters.** In bcc, and whenever the faces of the Voronoi cell have different areas, the Minkowski metrics differ from plain $q_l$. Do not use thresholds or reference values obtained with plain $q_l$.
- **The exponent changes the values.** In the fcc snapshot, an atom has on average {glue}`fcc_voronoi_neighbors` Voronoi neighbors instead of 12, because thermal motion creates small extra faces. With `voroexp=0` these count as much as the large faces, and the mean of $q_6$ drops to {glue:text}`fcc_q6_exp0:.3f`, compared with {glue}`fcc_q6_exp1` with `voroexp=1` and {glue}`fcc_q6_perfect` in the perfect crystal. Compare values only when they were computed with the same `voroexp`.
- **Averaging counts all Voronoi neighbors.** $\bar{q}_l$ is an unweighted average over all Voronoi neighbors, including those that share only a small face.

## References

1. W. Mickel, S. C. Kapfer, G. E. Schröder-Turk and K. Mecke, Shortcomings of the bond orientational order parameters for the analysis of disordered particulate matter, *J. Chem. Phys.* **138**, 044501 (2013). [doi:10.1063/1.4774084](https://doi.org/10.1063/1.4774084)
2. W. Lechner and C. Dellago, Accurate determination of crystal structures based on averaged local bond order parameters, *J. Chem. Phys.* **129**, 114707 (2008). [doi:10.1063/1.2977970](https://doi.org/10.1063/1.2977970)
