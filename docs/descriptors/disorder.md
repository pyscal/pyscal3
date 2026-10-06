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

# Disorder parameter

The disorder parameter of Kawasaki and Onuki [1] measures how much the orientation of the neighbor shell of an atom differs from those of its neighbors.
It is zero in a perfect crystal and large in a liquid.
It is computed from the vectors $q_{lm}$ of the [Steinhardt parameters](steinhardt), and it is used to map ordered and disordered regions, for example in glasses and polycrystals [1].
The same comparison of neighboring atoms underlies the [solid–liquid classification](solid_liquid).

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

For each atom $i$, the vector $q_{lm}(i)$, $m = -l, \dots, l$, is the average of the spherical harmonics over the bonds to its neighbors, as defined in [Steinhardt parameters](steinhardt).
The overlap of the vectors of two atoms $i$ and $j$ is

$$
S_{ij} = \frac{\mathrm{Re} \sum_{m=-l}^{l} q_{lm}(i)\, q_{lm}^*(j)}{\left( \sum_{m=-l}^{l} |q_{lm}(i)|^2 \right)^{1/2} \left( \sum_{m=-l}^{l} |q_{lm}(j)|^2 \right)^{1/2}},
$$

where the asterisk denotes the complex conjugate.
pyscal normalises the overlap, so $S_{ii} = 1$ and $S_{ij}$ lies between $-1$ and 1.
The disorder parameter of atom $i$ with $N(i)$ neighbors is

$$
D(i) = \frac{1}{N(i)} \sum_{j=1}^{N(i)} \left[ S_{ii} + S_{jj} - 2 S_{ij} \right] = 2 \left( 1 - \frac{1}{N(i)} \sum_{j=1}^{N(i)} S_{ij} \right).
$$

$D(i)$ lies between 0, when all neighbors have the same orientation as atom $i$, and 4.
Kawasaki and Onuki used $l = 6$ [1].

The averaged disorder parameter averages $D$ over the atom and its neighbors:

$$
\bar{D}(i) = \frac{1}{N(i) + 1} \sum_{k=0}^{N(i)} D(k),
$$

where $k = 0$ is the atom $i$ itself and $k = 1, \dots, N(i)$ are its neighbors.

## Usage

Find the neighbors first, then compute the disorder parameter:

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff="sann")

D = pyscal.disorder(atoms, q=6)
D_averaged = pyscal.disorder(atoms, q=6, averaged=True)
```

`q` sets the order $l$ (default 6).
`disorder` returns $D$, or $\bar{D}$ with `averaged=True`, with one value per atom.
It also stores the results on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_disorder"]` | $(N,)$ | $D$ |
| `atoms.arrays["pyscal_avg_disorder"]` | $(N,)$ | $\bar{D}$, with `averaged=True` |
| `atoms.arrays["pyscal_q6"]`, `["pyscal_q6_real"]`, `["pyscal_q6_imag"]` | $(N,)$, $(N, 13)$ | $q_6$ and the parts of $q_{6m}$, as in [Steinhardt parameters](steinhardt) |

If `atoms` already holds the $q_{lm}$ of the requested $l$, from `steinhardt_parameter` or `find_solids`, `disorder` uses them.
Otherwise it computes them from the current neighbors.

The sum over $S_{ij}$ in the definition is the average bond correlation $\langle s_{ij} \rangle$ of the [solid–liquid classification](solid_liquid), so $D(i) = 2 (1 - \langle s_{ij} \rangle)$ for the same neighbors and $l$.

## Crystal and liquid

We compute $D$ and $\bar{D}$ with $l = 6$ for three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.
The neighbors are found with SANN, which gives the full first shell in both crystals (see [Finding neighbors](../guide/neighbors)).

```{code-cell} ipython3
import numpy as np

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

values = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff="sann")
    averaged = pyscal.disorder(snapshot, q=6, averaged=True)
    values[name] = (snapshot.arrays["pyscal_disorder"], averaged)
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

# D = 2 (1 - <s_ij>) for the same neighbors and l
check = snapshots["liquid"].copy()
pyscal.find_neighbors(check, method="cutoff", cutoff="sann")
pyscal.find_solids(check, cluster=False)
assert np.allclose(pyscal.disorder(check, q=6), 2 * (1 - check.arrays["pyscal_avg_sij"]))

crystal_max = [max(values[name][k].max() for name in ("fcc", "bcc")) for k in (0, 1)]
liquid_min = [values["liquid"][k].min() for k in (0, 1)]
glue("crystal_D_max", round(crystal_max[0], 2), display=False)
glue("liquid_D_min", round(liquid_min[0], 2), display=False)
glue("crystal_Dbar_max", round(crystal_max[1], 2), display=False)
glue("liquid_Dbar_min", round(liquid_min[1], 2), display=False)
glue("overlap_D", int(sum((values[name][0] >= liquid_min[0]).sum() for name in ("fcc", "bcc"))
                      + (values["liquid"][0] <= crystal_max[0]).sum()), display=False)
glue("liquid_Dbar_mean", round(values["liquid"][1].mean(), 2), display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.25)
bins = np.linspace(0, 2.5, 51)
for k, title in enumerate(["(a)  $D$", r"(b)  $\bar{D}$"]):
    ax = mp[0, k]
    for name, (D, D_averaged) in values.items():
        data = (D, D_averaged)[k]
        ax.hist(data, bins=bins, density=True, histtype="stepfilled",
                color=COLOURS[name], alpha=0.4)
        ax.hist(data, bins=bins, density=True, histtype="step",
                color=COLOURS[name], lw=1.4, label=name)
    label(ax, title)
    ax.set_xlabel(title.split()[-1])
    ax.set_xlim(0, 2.5)
mp[0, 0].set_ylabel("Probability density")
mp[0, 1].legend(frameon=False, loc="upper right");
```

In panel (a), the largest value of $D$ in the two crystals is {glue}`crystal_D_max`, and the smallest value in the liquid is {glue}`liquid_D_min`.
The two distributions overlap only in their tails, which contain {glue}`overlap_D` atoms.
In panel (b), the averaged parameter $\bar{D}$ is at most {glue}`crystal_Dbar_max` in the crystals and at least {glue}`liquid_Dbar_min` in the liquid, and the distributions are well separated.
With $l = 6$, the fcc and the bcc crystal have almost the same distribution, so $D$ measures disorder but does not tell the crystal structures apart.

## Thermal vibrations and the choice of $l$

How does the disorder parameter grow with thermal motion, and does the choice of $l$ matter?
We give perfect fcc, bcc and hcp crystals random displacements with a standard deviation of up to 12 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations, and compute the mean of $\bar{D}$ for $l = 4$ and $l = 6$.

```{code-cell} ipython3
from ase.build import bulk

crystals = {
    "fcc": (bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6), 3.61 / np.sqrt(2)),
    "bcc": (bulk("Fe", "bcc", a=2.87, cubic=True).repeat(7), 2.87 * np.sqrt(3) / 2),
    "hcp": (bulk("Cu", "hcp", a=2.55).repeat((8, 8, 5)), 2.55),
}
noises = np.linspace(0, 0.12, 13)
ls = [4, 6]

scan = {}
for name, (perfect, nearest) in crystals.items():
    for noise in noises:
        crystal = perfect.copy()
        crystal.rattle(noise * nearest, seed=1)
        pyscal.find_neighbors(crystal, method="cutoff", cutoff="sann")
        for l in ls:
            scan[name, l, noise] = pyscal.disorder(crystal, q=l, averaged=True).mean()

liquid = snapshots["liquid"]
liquid_mean = {l: pyscal.disorder(liquid, q=l, averaged=True).mean() for l in ls}
```

```{code-cell} ipython3
:tags: [remove-cell]
glue("bcc_l4_004", round(scan["bcc", 4, noises[4]], 2), display=False)
glue("liquid_l4", round(liquid_mean[4], 2), display=False)
glue("liquid_l6", round(liquid_mean[6], 2), display=False)
glue("fcc_l6_012", round(scan["fcc", 6, noises[-1]], 2), display=False)
assert np.isclose(noises[4], 0.04)
# checks of the statements in the text
assert all(abs(scan[name, l, 0.0]) < 1e-8 for name in crystals for l in ls)
assert all(scan[name, 6, s] < liquid_mean[6] for name in crystals for s in noises)
assert all(np.all(np.diff([scan[name, 6, s] for s in noises]) > 0) for name in crystals)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.12)
for k, l in enumerate(ls):
    ax = mp[0, k]
    ax.axhline(liquid_mean[l], ls="--", color=DARK, lw=1)
    ax.annotate("liquid", (0, liquid_mean[l]), xytext=(2, 4), textcoords="offset points",
                fontsize=9, color=DARK)
    for name, marker in zip(crystals, "oDs"):
        ax.plot(100 * noises, [scan[name, l, s] for s in noises], marker=marker,
                color=COLOURS[name], mec=DARK, mew=0.7, ms=5, lw=1.6, label=name)
    label(ax, f"({'ab'[k]})  $l = {l}$")
    ax.set_xlabel("Displacement  (% of $d_1$)")
    ax.set_ylim(-0.05, 2.3)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel(r"Mean $\bar{D}$")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

The standard deviation of the displacements is given in percent of the nearest neighbor distance $d_1$.
The dashed lines mark the mean of $\bar{D}$ in the liquid snapshot.
In all three crystals, $\bar{D}$ is zero without displacements.
With $l = 6$ (b), it increases steadily with the displacements and stays below the value of the liquid over the whole range.
With $l = 4$ (a), the bcc crystal reaches {glue}`bcc_l4_004` already at a standard deviation of 4 %, close to the value of the liquid, {glue}`liquid_l4`.
The reason is that $q_4$ of a perfect bcc crystal with its full first shell is close to zero (see [Steinhardt parameters](steinhardt)).
The vectors $q_{4m}$ are then short, and small displacements change their direction completely.
Use $l = 6$ unless there is a reason for another choice.

## Things to watch

- **Stored $q_{lm}$ are reused.** `disorder` uses the $q_{lm}$ stored on `atoms` if there are any. After calling `find_neighbors` again with other settings, call `steinhardt_parameter` with the same $l$ before `disorder`, so that the $q_{lm}$ belong to the current neighbors.
- **The choice of $l$.** With $l = 4$, bcc crystals at finite temperature look as disordered as a liquid, as shown above.
- **The neighbor method changes the values.** $D$ depends on which atoms are counted as neighbors, through $q_{lm}$ and through the sum over neighbors. Compare values only when they were computed with the same neighbor method.

## References

1. T. Kawasaki and A. Onuki, Construction of a disorder variable from Steinhardt order parameters in binary mixtures at high densities in three dimensions, *J. Chem. Phys.* **135**, 174109 (2011). [doi:10.1063/1.3656762](https://doi.org/10.1063/1.3656762)
