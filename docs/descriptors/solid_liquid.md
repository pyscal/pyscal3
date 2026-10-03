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

# Solid–liquid classification

`find_solids` labels each atom as solid or liquid.
It compares the orientation of the neighbor shell of an atom with those of its neighbors, using the vectors $q_{6m}$ of the [Steinhardt parameters](steinhardt) [1, 2].
In a crystal, neighboring atoms have similar neighbor shells, and in a liquid they do not.
`find_clusters` then groups the solid atoms, or any other selection of atoms, into connected clusters.
Together they are used to follow crystal nucleation and growth, for example by measuring the size of the largest crystalline cluster in a liquid.
The [disorder parameter](disorder) is built from the same comparison, and the [entropy parameter](entropy) and the [averaged Steinhardt parameters](steinhardt) are other ways to separate solid from liquid.

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
The bond correlation between atom $i$ and its neighbor $j$ is the normalised scalar product of their vectors:

$$
s_{ij} = \frac{\mathrm{Re} \sum_{m=-l}^{l} q_{lm}(i)\, q_{lm}^*(j)}{\left( \sum_{m=-l}^{l} |q_{lm}(i)|^2 \right)^{1/2} \left( \sum_{m=-l}^{l} |q_{lm}(j)|^2 \right)^{1/2}},
$$

where the asterisk denotes the complex conjugate.
$s_{ij}$ lies between $-1$ and 1, and it is 1 when the neighbor shells of $i$ and $j$ have the same orientation.
pyscal uses $l = 6$ by default.

A bond is solid if $s_{ij}$ is larger than a threshold $s_\mathrm{c}$.
Atom $i$ is solid if both of the following hold:

1. The number $n_\mathrm{s}(i)$ of its solid bonds is larger than a minimum $n_\mathrm{c}$, or the fraction $n_\mathrm{s}(i) / N(i)$ is larger than a minimum $f_\mathrm{c}$, where $N(i)$ is the number of neighbors of $i$ [1].
2. The average bond correlation is larger than a second threshold $\bar{s}_\mathrm{c}$ [2]:

$$
\langle s_{ij} \rangle = \frac{1}{N(i)} \sum_{j=1}^{N(i)} s_{ij} > \bar{s}_\mathrm{c}.
$$

The second condition improves the classification at the boundary between solid and liquid [2].
All comparisons are strict.

Two solid atoms belong to the same cluster if they are connected by a chain of neighbor bonds between solid atoms.

## Usage

Find the neighbors first, then classify the atoms:

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)

largest = pyscal.find_solids(atoms, bonds=0.5, threshold=0.5, avgthreshold=0.6)
largest
```

The arguments are those of the definition:

| Argument | Default | Meaning |
|---|---|---|
| `q` | 6 | the order $l$ of the vectors $q_{lm}$ |
| `threshold` | 0.5 | $s_\mathrm{c}$, a bond is solid if $s_{ij} >$ `threshold` |
| `bonds` | 0.5 | an `int` sets $n_\mathrm{c}$, a `float` between 0 and 1 sets $f_\mathrm{c}$ |
| `avgthreshold` | 0.6 | $\bar{s}_\mathrm{c}$ |
| `cluster` | `True` | cluster the solid atoms and return the size of the largest cluster |
| `cutoff` | 0 | the cluster cutoff, see below |
| `right` | `True` | with `False`, the three comparisons are reversed ($<$ instead of $>$) |

With `cluster=True`, `find_solids` returns the number of atoms in the largest solid cluster, and `None` otherwise.
It stores these keys on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_solid"]` | $(N,)$ | 1.0 for solid atoms, 0.0 for liquid atoms |
| `atoms.arrays["pyscal_bonds"]` | $(N,)$ | number of solid bonds $n_\mathrm{s}(i)$ |
| `atoms.arrays["pyscal_avg_sij"]` | $(N,)$ | $\langle s_{ij} \rangle$ |
| `atoms.arrays["pyscal_sij"]` or `atoms.info["pyscal_sij"]` | $(N, k)$, or a list of lists | $s_{ij}$ for each neighbor, in the order of the neighbor list |
| `atoms.arrays["pyscal_q6"]`, `["pyscal_q6_real"]`, `["pyscal_q6_imag"]` | $(N,)$, $(N, 13)$ | $q_6$ and the parts of $q_{6m}$, as in [Steinhardt parameters](steinhardt) |
| `atoms.arrays["pyscal_cluster"]` | $(N,)$ | cluster number, starting at 1, and $-1$ for liquid atoms |
| `atoms.arrays["pyscal_largest_cluster"]` | $(N,)$ | `True` for the atoms of the largest cluster |

`pyscal_sij` is in `atoms.arrays` when every atom has the same number $k$ of neighbors, and in `atoms.info` otherwise.
The two cluster keys are written only with `cluster=True`.

`find_clusters` clusters any selection of atoms.
It takes one boolean per atom and stores the same two cluster keys:

```{code-cell} ipython3
solid = atoms.arrays["pyscal_solid"] > 0
pyscal.find_clusters(atoms, condition=solid)
```

By default, clusters follow every bond of the neighbor list.
With `cutoff` $r > 0$, in `find_solids` or `find_clusters`, only bonds shorter than or equal to $r$ are followed.

## Separating solid and liquid

How well do the default thresholds separate solid from liquid atoms?
We classify three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.

```{code-cell} ipython3
import numpy as np

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

results = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=0)
    pyscal.find_solids(snapshot, cluster=False)
    n_neighbors = np.diff(snapshot.info["pyscal_bond_offsets"])
    sij = snapshot.arrays.get("pyscal_sij", snapshot.info.get("pyscal_sij"))
    results[name] = {
        "sij": np.concatenate([np.asarray(row) for row in sij]),
        "fraction": snapshot.arrays["pyscal_bonds"] / n_neighbors,
        "average": snapshot.arrays["pyscal_avg_sij"],
        "solid": snapshot.arrays["pyscal_solid"] > 0,
    }
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

for name, values in results.items():
    glue(f"solid_{name}", int(values["solid"].sum()), display=False)
    glue(f"n_{name}", len(values["solid"]), display=False)
assert results["bcc"]["solid"].all() and not results["liquid"]["solid"].any()
glue("liquid_bonds", round(100 * np.mean(results["liquid"]["sij"] > 0.5)), display=False)
glue("liquid_fraction", round(100 * np.mean(results["liquid"]["fraction"] > 0.5), 1), display=False)
glue("liquid_average_max", round(results["liquid"]["average"].max(), 2), display=False)
glue("crystal_average_min",
     round(min(results[name]["average"].min() for name in ("fcc", "bcc")), 2), display=False)
glue("crystal_average_below",
     int(sum((results[name]["average"] <= 0.6).sum() for name in ("fcc", "bcc"))), display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
panels = [
    ("sij", "(a)  bond correlation", "$s_{ij}$", np.linspace(-0.6, 1, 41), 0.5),
    ("fraction", "(b)  solid bonds", r"$n_\mathrm{s}(i)\,/\,N(i)$", np.linspace(0, 1, 21), 0.5),
    ("average", "(c)  average", r"$\langle s_{ij} \rangle$", np.linspace(-0.2, 1, 37), 0.6),
]
mp = figure(columns=3, ratio=0.36, wspace=0.35)
for k, (key, title, xlabel, bins, threshold) in enumerate(panels):
    ax = mp[0, k]
    for name, values in results.items():
        ax.hist(values[key], bins=bins, density=True, histtype="stepfilled",
                color=COLOURS[name], alpha=0.4)
        ax.hist(values[key], bins=bins, density=True, histtype="step",
                color=COLOURS[name], lw=1.4, label=name)
    ax.axvline(threshold, ls="--", color=DARK, lw=1)
    label(ax, title)
    ax.set_xlabel(xlabel)
mp[0, 0].set_ylabel("Probability density")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

The dashed lines mark the default thresholds.
Panel (a) shows that most bonds in the crystals have $s_{ij}$ close to 1.
In the liquid, $s_{ij}$ is spread over a wide range, and {glue}`liquid_bonds` % of the bonds are solid bonds.
A single solid bond therefore says little.
Panel (b) shows the fraction of solid bonds of each atom, which is larger than 0.5 for only {glue}`liquid_fraction` % of the liquid atoms.
Panel (c) shows the average $\langle s_{ij} \rangle$.
Its largest value in the liquid is {glue}`liquid_average_max`, below the threshold of 0.6, and only {glue}`crystal_average_below` atoms of the two crystals lie below the threshold.

With both conditions, {glue}`solid_fcc` of the {glue}`n_fcc` fcc atoms, all {glue}`n_bcc` bcc atoms and none of the {glue}`n_liquid` liquid atoms are labelled solid.

## A crystal nucleus in a liquid

`find_solids` is most often used for a crystal that grows from its liquid.
As a model, we take the liquid aluminium snapshot, repeated twice in each direction, and replace its centre by a sphere of fcc aluminium with a radius of 9 Å.
The liquid atoms are removed up to 10 Å from the centre, so that liquid and crystal atoms do not overlap.
The crystal atoms get random displacements with a standard deviation of 0.15 Å in each direction, a rough model of thermal vibrations.

```{code-cell} ipython3
from ase.build import bulk

liquid = read("conf.lqd.Al.dump", format="lammps-dump-text").repeat(2)
liquid.wrap()
centre = liquid.cell.lengths() / 2
radius = 9.0

crystal = bulk("Al", "fcc", a=4.05, cubic=True).repeat(6)
crystal.positions += centre - crystal.cell.lengths() / 2
crystal = crystal[np.linalg.norm(crystal.positions - centre, axis=1) < radius]
crystal.rattle(0.15, seed=1)

liquid = liquid[np.linalg.norm(liquid.positions - centre, axis=1) > radius + 1.0]
system = liquid + crystal
inserted = np.arange(len(system)) >= len(liquid)

pyscal.find_neighbors(system, method="cutoff", cutoff=0)
pyscal.find_solids(system)
```

```{code-cell} ipython3
:tags: [remove-cell]
solid = system.arrays["pyscal_solid"] > 0
distance = np.linalg.norm(system.positions - centre, axis=1)
missed = inserted & ~solid
glue("n_system", len(system), display=False)
glue("n_inserted", int(inserted.sum()), display=False)
n_largest = int(system.arrays["pyscal_largest_cluster"].sum())
glue("n_largest", n_largest, display=False)
glue("n_missed", int(missed.sum()), display=False)
glue("missed_rmin", round(distance[missed].min(), 1), display=False)
outer = inserted & (distance > 8)
glue("outer_solid", round(100 * solid[outer].mean()), display=False)
sann = system.copy()
pyscal.find_neighbors(sann, method="cutoff", cutoff="sann")
glue("sann_largest", pyscal.find_solids(sann), display=False)
glue("sann_false", int(np.sum((sann.arrays["pyscal_solid"] > 0) & ~inserted)), display=False)

# checks of the statements in the text
assert not np.any(solid & ~inserted)
assert solid[inserted & (distance <= 8)].all()
assert not np.any((distance > 1) & (distance < 2))
assert not np.any((distance > 9) & (distance < 10))
```

The system has {glue}`n_system` atoms, {glue}`n_inserted` of them in the crystal.
The largest solid cluster has {glue}`n_largest` atoms.
No liquid atom is labelled solid, and the {glue}`n_missed` crystal atoms labelled liquid all lie more than {glue}`missed_rmin` Å from the centre, at the surface of the sphere.

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.45, wspace=0.3, width_ratios=[1, 1.3])

# (a) a slice through the centre of the sphere
ax = mp[0, 0]
in_slice = np.abs(system.positions[:, 2] - centre[2]) < 2.0
x, y = (system.positions[:, :2] - centre[:2]).T
groups = [
    (~solid & ~inserted, dict(color=COLOURS["liquid"], ec=DARK), "labelled liquid"),
    (solid, dict(color=COLOURS["fcc"], ec=DARK), "labelled solid"),
    (missed, dict(color="white", ec=COLOURS["fcc"]), "crystal atom labelled liquid"),
]
for selection, style, text in groups:
    sel = selection & in_slice
    ax.scatter(x[sel], y[sel], s=16, lw=0.8, label=text, **style)
angle = np.linspace(0, 2 * np.pi, 200)
ax.plot(radius * np.cos(angle), radius * np.sin(angle), ls="--", color=DARK, lw=1)
label(ax, "(a)  slice through the centre")
ax.set_xlabel("$x$  (Å)")
ax.set_ylabel("$y$  (Å)")
ax.set_xlim(-16, 16)
ax.set_ylim(-16, 16)
ax.set_aspect("equal")

# (b) fraction of solid atoms against the distance from the centre
ax = mp[0, 1]
edges = np.arange(0, 17, 1.0)
shell = np.digitize(distance, edges) - 1
middle = 0.5 * (edges[1:] + edges[:-1])
frac_solid = [solid[shell == k].mean() if np.any(shell == k) else np.nan
              for k in range(len(middle))]
ax.plot(middle, frac_solid, ls="-", marker="o", color=COLOURS["fcc"], mec=DARK,
        mew=0.7, ms=5, lw=1.6)
ax.axvline(radius, ls="--", color=DARK, lw=1)
ax.annotate("surface of\nthe crystal", (radius, 0.6), xytext=(8, 0),
            textcoords="offset points", fontsize=9, color=DARK, va="center")
label(ax, "(b)  profile")
ax.set_xlabel("Distance from the centre  (Å)")
ax.set_ylabel("Fraction labelled solid")
ax.set_ylim(-0.03, 1.08)
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1), markerscale=1.5);
```

Panel (a) shows the atoms within 2 Å of a plane through the centre, and the dashed circle marks the surface of the crystal.
Panel (b) shows the fraction of atoms labelled solid in shells 1 Å thick around the centre.
The shells between 1 and 2 Å and between 9 and 10 Å contain no atoms.
Every crystal atom within 8 Å of the centre is labelled solid, but only {glue}`outer_solid` % of the crystal atoms farther out.
These atoms at the surface have many liquid neighbors, and they fail the conditions.
The classification therefore places the boundary of the crystal slightly inside its true surface.

`find_clusters` works with any condition.
For example, the cluster of atoms with an [averaged Steinhardt parameter](steinhardt) $\bar{q}_6 > 0.4$:

```{code-cell} ipython3
q6_averaged = pyscal.steinhardt_parameter(system, l=6, averaged=True)[0]
pyscal.find_clusters(system, condition=q6_averaged > 0.4)
```

```{code-cell} ipython3
:tags: [remove-cell]
q6_cluster = int(system.arrays["pyscal_largest_cluster"].sum())
assert q6_cluster < n_largest
glue("q6_cluster", q6_cluster, display=False)
```

With {glue}`q6_cluster` atoms, this cluster is smaller than the one from `find_solids`.
$\bar{q}_6$ averages over the neighbors of an atom, and near the surface these include liquid atoms.

## Things to watch

- **`bonds` is a number or a fraction.** An `int` is a number of bonds and a `float` is a fraction. The comparison is strict. `bonds=6` requires at least 7 solid bonds, and `bonds=1.0` labels no atom as solid.
- **The thresholds depend on the neighbors.** $s_{ij}$ is computed from the $q_{lm}$ of the current neighbor list, and the distributions above change with the neighbor method. In the model of a nucleus above, SANN neighbors instead of the adaptive cutoff label {glue}`sann_false` liquid atoms as solid, and the largest cluster has {glue}`sann_largest` atoms instead of {glue}`n_largest`. Check the distributions, as in the figure above, before using the default thresholds for a new system.
- **Interfaces.** Atoms at the surface of a crystal are often labelled liquid. Cluster sizes are therefore smaller than the number of atoms in the crystal.
- **The cluster cutoff only removes bonds.** Clusters follow the bonds of the neighbor list. A cluster cutoff larger than the neighbor cutoff adds no bonds.
- **Stored $q_6$.** `find_solids` recomputes $q_{6m}$ from the current neighbors and overwrites `pyscal_q6`, `pyscal_q6_real` and `pyscal_q6_imag`.

## References

1. S. Auer and D. Frenkel, Numerical simulation of crystal nucleation in colloids, in *Advanced Computer Simulation*, edited by C. Holm and K. Kremer, Advances in Polymer Science, pp. 149–208 (Springer, Berlin, 2005). [doi:10.1007/b99429](https://doi.org/10.1007/b99429)
2. J. Bokeloh, G. Wilde, R. E. Rozas, R. Benjamin and J. Horbach, Nucleation barriers for the liquid-to-crystal transition in simple metals: Experiment vs. simulation, *Eur. Phys. J. Spec. Top.* **223**, 511 (2014). [doi:10.1140/epjst/e2014-02106-2](https://doi.org/10.1140/epjst/e2014-02106-2)
