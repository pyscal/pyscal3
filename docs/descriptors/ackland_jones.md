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

# Ackland–Jones classification

The method of Ackland and Jones [1] labels each atom as fcc, hcp, bcc, icosahedral or other, from the angles between the bonds to its neighbors.
It counts the bond angles in a few ranges, the [χ parameters](angular), and compares the counts with those of the perfect structures.
It gives the same kind of labels as [common neighbor analysis](cna) (CNA).
This page describes the rules that pyscal uses and compares the labels with CNA on crystals at finite temperature and on a liquid.

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

For each atom, the χ parameters $\chi_0, \dots, \chi_8$ count the pairs of neighbors whose bond angle falls in nine ranges of $\cos\theta$ (see [Angular criteria and χ parameters](angular) for the ranges).
In the perfect structures, with first shell neighbors, the counts that matter are:

| Structure | Neighbors | $\chi_0$ | $\chi_1 + \chi_2 + \chi_3$ | $\chi_2$ | $\chi_5$ |
|---|---|---|---|---|---|
| bcc | 14 | 7 | 0 | 0 | 12 |
| fcc | 12 | 6 | 0 | 0 | 12 |
| icosahedral | 12 | 6 | 0 | 0 | 0 |
| hcp | 12 | 3 | 6 | 6 | 12 |

$\chi_0$ counts pairs of neighbors on nearly opposite sides of the atom (angles above 160.9°), $\chi_1$ to $\chi_3$ count angles between 134.8° and 160.9°, and $\chi_5$ counts angles between 78.8° and 101.2°.
pyscal applies the following rules in this order:

1. **bcc** if $\chi_0 \ge 7$.
2. Otherwise, if $\chi_0 \ge 5$ and $\chi_1 + \chi_2 + \chi_3 = 0$: **icosahedral** if $\chi_5 = 0$, and **fcc** if $\chi_5 > 0$.
3. Otherwise, **hcp** if $\chi_2 > 0$.
4. Otherwise, **other**.

The method of the original paper [1] differs in two ways.
It selects the neighbors itself, from the distances to the six nearest atoms.
It also assigns atoms whose counts do not match the ideal values to the structure with the smallest deviation from them.
pyscal uses the stored neighbor list and only the rules above.

## Usage

Find the neighbors first, then classify:

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.io import read

atoms = read("conf.bcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=3.8)

labels, names = pyscal.identify_ackland_jones(atoms)
labels[:6], names[:6]
```

The cutoff of 3.8 Å is the first minimum of the radial distribution function of this bcc crystal, between the second and the third neighbor shell (see [Finding neighbors](../guide/neighbors)).
`identify_ackland_jones` returns an integer label and a name for each atom:

| Label | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| Name | other | fcc | hcp | bcc | ico |

The numbers are the same as those of `common_neighbor_analysis`.
The results are stored on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_ackland_label"]` | $(N,)$ | integer label |
| `atoms.arrays["pyscal_structure"]` | $(N,)$ | name, as a string |
| `atoms.arrays["pyscal_chiparams"]` | $(N, 9)$ | $\chi_0, \dots, \chi_8$ |

`identify_ackland_jones` uses the stored neighbor list and raises an error if `find_neighbors` has not been called.

## Crystals at finite temperature and a liquid

How well does the method label atoms in real simulations?
We classify three MD snapshots from the `examples` folder, an fcc crystal and a bcc crystal at finite temperature, and a liquid, with both methods.
The neighbors are all atoms within the first minimum of the radial distribution function, which gives 12 neighbors in fcc and 14 in bcc.

```{code-cell} ipython3
from collections import Counter

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}
cutoffs = {"fcc": 3.5, "bcc": 3.8, "liquid": 4.1}

fractions = {}
for name, snapshot in snapshots.items():
    n = len(snapshot)
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=cutoffs[name])
    _, names = pyscal.identify_ackland_jones(snapshot)
    fractions[name, "Ackland–Jones"] = {k: v / n for k, v in Counter(names).items()}
    counts = pyscal.common_neighbor_analysis(snapshot)
    fractions[name, "CNA"] = {("other" if k == "others" else k): v / n
                              for k, v in counts.items()}
```

`common_neighbor_analysis` finds its own neighbors, and it calls the label 0 `others`, which is renamed here to `other`.

```{code-cell} ipython3
:tags: [hide-input]
structures = [s for s in ("fcc", "hcp", "bcc", "ico", "other")
              if any(f.get(s, 0) > 0 for f in fractions.values())]
colours = {**COLOURS, "other": COLOURS["others"]}
rows = [(name, method) for name in snapshots for method in ("Ackland–Jones", "CNA")]
y = np.array([0, 1, 2.6, 3.6, 5.2, 6.2])[::-1]

mp = figure(ratio=0.45)
ax = mp[0, 0]
left = np.zeros(len(rows))
for structure in structures:
    width = np.array([fractions[row].get(structure, 0) for row in rows])
    ax.barh(y, width, left=left, height=0.8, color=colours[structure], ec=DARK, lw=0.6,
            label=structure)
    left += width
ax.set_yticks(y, [method for _, method in rows])
for name, centre in zip(snapshots, [y[0:2].mean(), y[2:4].mean(), y[4:6].mean()]):
    ax.text(-0.22, centre, name, ha="right", va="center", fontsize=11,
            transform=ax.get_yaxis_transform())
ax.set_xlim(0, 1)
ax.set_xlabel("Fraction of atoms")
ax.legend(frameon=False, ncol=len(structures), loc="upper center",
          bbox_to_anchor=(0.45, -0.22));
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

def percent(name, method, structure):
    return round(100 * fractions[name, method].get(structure, 0))

for name, structures_ in [("fcc", ["fcc", "hcp", "other"]), ("bcc", ["bcc", "hcp"])]:
    for structure in structures_:
        glue(f"aj_{name}_{structure}", percent(name, "Ackland–Jones", structure), display=False)
glue("aj_liquid_hcp", round(100 * fractions["liquid", "Ackland–Jones"]["hcp"], 1), display=False)
glue("cna_fcc_fcc", percent("fcc", "CNA", "fcc"), display=False)
glue("cna_fcc_other", percent("fcc", "CNA", "other"), display=False)
glue("cna_bcc_bcc", percent("bcc", "CNA", "bcc"), display=False)
glue("cna_liquid_other", percent("liquid", "CNA", "other"), display=False)
chi_liquid = snapshots["liquid"].arrays["pyscal_chiparams"]
glue("liquid_chi2", round(100 * np.mean(chi_liquid[:, 2] > 0), 1), display=False)
```

In the fcc snapshot, the Ackland–Jones method labels {glue}`aj_fcc_fcc` % of the atoms as fcc, {glue}`aj_fcc_hcp` % as hcp and {glue}`aj_fcc_other` % as other.
CNA labels {glue}`cna_fcc_fcc` % as fcc and {glue}`cna_fcc_other` % as other.
In the bcc snapshot, the method labels {glue}`aj_bcc_bcc` % as bcc and {glue}`aj_bcc_hcp` % as hcp, and CNA labels {glue}`cna_bcc_bcc` % as bcc.
In the liquid, the method labels {glue}`aj_liquid_hcp` % of the atoms as hcp, and CNA labels {glue}`cna_liquid_other` % as other.

The hcp labels come from rule 3.
Thermal motion moves angles out of $\chi_0$ or into $\chi_1$ to $\chi_3$, so that an fcc atom fails the fcc test.
Any atom that fails the tests for bcc, fcc and icosahedral and has an angle in $\chi_2$ is labelled hcp.
In the liquid, {glue}`liquid_chi2` % of the atoms have such an angle.
CNA labels atoms that it cannot assign as other, so a missed atom does not look like a different crystal structure.

## Random displacements

To see how the labels change with the size of the thermal displacements, we add random displacements to perfect fcc, hcp and bcc crystals, with a standard deviation $\sigma$ up to 10 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations.
The neighbors are all atoms within a cutoff between the first and second neighbor shell (fcc, hcp), or between the second and third (bcc).
We compare with adaptive CNA.

```{code-cell} ipython3
from ase.build import bulk

def perfect(name):
    """Crystal, nearest neighbor distance and neighbor cutoff."""
    if name == "hcp":
        return bulk("Cu", "hcp", a=2.55).repeat((8, 8, 5)), 2.55, 0.5 * (1 + np.sqrt(2)) * 2.55
    if name == "fcc":
        a = 3.61
        return bulk("Cu", "fcc", a=a, cubic=True).repeat(6), a / np.sqrt(2), 0.5 * (1 / np.sqrt(2) + 1) * a
    a = 2.87
    return bulk("Fe", "bcc", a=a, cubic=True).repeat(6), a * np.sqrt(3) / 2, 0.5 * (1 + np.sqrt(2)) * a

sigmas = np.linspace(0, 0.10, 11)
scan = {}
for name in ("fcc", "hcp", "bcc"):
    crystal, nearest, cutoff = perfect(name)
    for sigma in sigmas:
        noisy = crystal.copy()
        noisy.rattle(sigma * nearest, seed=1)
        pyscal.find_neighbors(noisy, method="cutoff", cutoff=cutoff)
        _, names = pyscal.identify_ackland_jones(noisy)
        scan[name, sigma, "Ackland–Jones"] = Counter(names)
        scan[name, sigma, "CNA"] = pyscal.common_neighbor_analysis(noisy)
        scan[name, sigma, "n"] = len(noisy)
```

```{code-cell} ipython3
:tags: [hide-input]
from matplotlib.lines import Line2D

mp = figure(columns=3, ratio=0.36, wspace=0.12)
for k, name in enumerate(("fcc", "hcp", "bcc")):
    ax = mp[0, k]
    n = scan[name, 0, "n"]
    for structure, marker in zip(("fcc", "hcp", "bcc", "other"), "o^sD"):
        values = [scan[name, s, "Ackland–Jones"].get(structure, 0) / n for s in sigmas]
        colour = COLOURS["others"] if structure == "other" else COLOURS[structure]
        ax.plot(100 * sigmas, values, marker=marker, color=colour, mec=DARK, mew=0.7,
                ms=4.5, lw=1.5)
    cna = [scan[name, s, "CNA"][name] / n for s in sigmas]
    ax.plot(100 * sigmas, cna, ls=":", marker="o", color=COLOURS[name], mfc="white",
            mec=COLOURS[name], mew=1, ms=4.5, lw=1.5)
    label(ax, f"({'abc'[k]})  {name}")
    ax.set_ylim(-0.03, 1.08)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("Fraction of atoms")
mp[0, 1].set_xlabel(r"$\sigma$  (% of nearest neighbor distance)")
handles = [Line2D([], [], marker=m, color=COLOURS[s] if s != "other" else COLOURS["others"],
                  mec=DARK, mew=0.7, ms=5, lw=1.5, label=f"{s} (Ackland–Jones)")
           for s, m in zip(("fcc", "hcp", "bcc", "other"), "o^sD")]
handles.append(Line2D([], [], ls=":", marker="o", color=DARK, mfc="white", mec=DARK,
                      ms=5, lw=1.5, label="correct label (CNA)"))
mp.fig.legend(handles=handles, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
def fraction(name, sigma, method, structure):
    return round(100 * scan[name, sigma, method].get(structure, 0) / scan[name, sigma, "n"])

glue("lowest_4", min(fraction(name, s, method, name) for name in ("fcc", "hcp", "bcc")
                     for method in ("Ackland–Jones", "CNA") for s in sigmas[:5]), display=False)
s6 = sigmas[6]
glue("aj_fcc_6", fraction("fcc", s6, "Ackland–Jones", "fcc"), display=False)
glue("cna_fcc_6", fraction("fcc", s6, "CNA", "fcc"), display=False)
glue("aj_bcc_6", fraction("bcc", s6, "Ackland–Jones", "bcc"), display=False)
glue("cna_bcc_6", fraction("bcc", s6, "CNA", "bcc"), display=False)
s10 = sigmas[10]
glue("aj_fcc_hcp_10", fraction("fcc", s10, "Ackland–Jones", "hcp"), display=False)
glue("aj_bcc_hcp_10", fraction("bcc", s10, "Ackland–Jones", "hcp"), display=False)
glue("aj_hcp_10", fraction("hcp", s10, "Ackland–Jones", "hcp"), display=False)
```

The solid lines show the fraction of atoms that the Ackland–Jones method assigns to each structure, and the dotted lines the fraction that CNA labels correctly.
For $\sigma$ up to 4 %, both methods label at least {glue}`lowest_4` % of the atoms correctly.
At $\sigma = 6$ %, the method labels {glue}`aj_fcc_6` % of the fcc atoms as fcc and {glue}`aj_bcc_6` % of the bcc atoms as bcc, where CNA labels {glue}`cna_fcc_6` % and {glue}`cna_bcc_6` %.
At $\sigma = 10$ %, the method labels {glue}`aj_fcc_hcp_10` % of the fcc atoms and {glue}`aj_bcc_hcp_10` % of the bcc atoms as hcp.
The hcp crystal stays at {glue}`aj_hcp_10` % hcp, because hcp is the label for atoms that fail the other tests.
A high fraction of hcp atoms is therefore not evidence for an hcp phase.
Check it with CNA or with the [averaged Steinhardt parameters](steinhardt).

```{code-cell} ipython3
:tags: [remove-cell]
# checks for the statements below
def share(atoms, structure):
    names = np.array(pyscal.identify_ackland_jones(atoms)[1])
    return round(100 * np.mean(names == structure))

pyscal.find_neighbors(snapshots["bcc"], method="cutoff", cutoff=0)
glue("adaptive_bcc", share(snapshots["bcc"], "bcc"), display=False)
pyscal.find_neighbors(snapshots["fcc"], method="voronoi")
glue("voronoi_fcc", share(snapshots["fcc"], "fcc"), display=False)

bcc = bulk("Fe", "bcc", a=2.87, cubic=True).repeat(4)
pyscal.find_neighbors(bcc, method="number", nmax=8)
assert share(bcc, "other") == 100
assert np.all(bcc.arrays["pyscal_chiparams"][:, 0] == 4)

diamond = bulk("Si", "diamond", a=5.43, cubic=True).repeat(3)
pyscal.find_neighbors(diamond, method="number", nmax=4)
assert share(diamond, "other") == 100
```

## Things to watch

- **The neighbor list decides the result.** The rules expect the 12 first neighbors in fcc, hcp and icosahedral environments, and 14 (8 + 6) in bcc. With only the 8 first neighbors, a bcc atom has $\chi_0 = 4$ and is labelled other. Check the number of neighbors per atom before classifying. Methods that add or drop neighbors give fewer correct labels (see [Finding neighbors](../guide/neighbors)). With the adaptive cutoff and its defaults, only {glue}`adaptive_bcc` % of the atoms of the bcc snapshot are labelled bcc. With Voronoi neighbors, only {glue}`voronoi_fcc` % of the atoms of the fcc snapshot are labelled fcc.
- **hcp is the fallback label.** Every atom with an angle in $\chi_2$ that is not labelled bcc, fcc or icosahedral is labelled hcp, including the atoms of a liquid.
- **Diamond structures are not covered.** With their four first neighbors, the atoms of a perfect diamond crystal are labelled other. Use `diamond_structure` (see [Common neighbor analysis](cna.md#diamond-structures)).
- **`pyscal_structure` is shared with CNA.** `common_neighbor_analysis` stores integer labels under the same key, and `identify_ackland_jones` stores names. The function called last overwrites the result of the other. The integer labels of the Ackland–Jones method are also in `pyscal_ackland_label`.

## References

1. G. J. Ackland and A. P. Jones, Applications of local crystal structure measures in experiment and simulation, *Phys. Rev. B* **73**, 054104 (2006). [doi:10.1103/PhysRevB.73.054104](https://doi.org/10.1103/PhysRevB.73.054104)
