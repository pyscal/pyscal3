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

The method of Ackland and Jones [1] labels each atom as fcc, hcp, bcc, icosahedral or unknown, from the angles between the bonds to its nearest neighbors.
It counts the bond angles in eight ranges and compares the counts with those of the perfect structures.
It gives the same kind of labels as [common neighbor analysis](cna) (CNA), and it chooses its neighbors from the distances to the six nearest atoms, so it needs no cutoff.
This page describes the method and compares its labels with CNA on crystals at finite temperature and on a liquid.

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

For each atom $i$, let $r_0^2$ be the mean of the squared distances to its six nearest atoms.
The method uses two groups of neighbors:
the $N_0$ atoms with $r_{ij}^2 < 1.45\, r_0^2$, and the $N_1$ atoms with $r_{ij}^2 < 1.55\, r_0^2$.
For every pair $j, k$ of the $N_0$ neighbors, the cosine of the bond angle $\theta_{jik}$ is counted in one of eight ranges:

| Count | $\chi_0$ | $\chi_1$ | $\chi_2$ | $\chi_3$ | $\chi_4$ | $\chi_5$ | $\chi_6$ | $\chi_7$ |
|---|---|---|---|---|---|---|---|---|
| $\cos\theta$ from | −1 | −0.945 | −0.915 | −0.755 | −0.195 | 0.195 | 0.245 | 0.795 |
| to | −0.945 | −0.915 | −0.755 | −0.195 | 0.195 | 0.245 | 0.795 | 1 |

In the perfect structures the counts are:

| Structure | $N_0$ | $\chi_0$ | $\chi_1$ | $\chi_2$ | $\chi_3$ | $\chi_4$ | $\chi_5$ | $\chi_6$ | $\chi_7$ |
|---|---|---|---|---|---|---|---|---|---|
| fcc | 12 | 6 | 0 | 0 | 24 | 12 | 0 | 24 | 0 |
| hcp | 12 | 3 | 0 | 6 | 21 | 12 | 0 | 24 | 0 |
| bcc | 14 | 7 | 0 | 0 | 36 | 12 | 0 | 36 | 0 |
| icosahedral | 12 | 6 | 0 | 0 | 30 | 0 | 0 | 30 | 0 |

From the counts, the method computes the deviations from a bcc, a close packed, an fcc and an hcp environment:

$$
\delta_\mathrm{bcc} = \frac{0.35\, \chi_4}{\chi_5 + \chi_6 - \chi_4}, \quad
\delta_\mathrm{cp} = \left| 1 - \frac{\chi_6}{24} \right|, \quad
\delta_\mathrm{fcc} = \frac{0.61}{6} \left( |\chi_0 + \chi_1 - 6| + \chi_2 \right), \quad
\delta_\mathrm{hcp} = \frac{1}{12} \left( |\chi_0 - 3| + |\chi_0 + \chi_1 + \chi_2 + \chi_3 - 9| \right).
$$

$\delta_\mathrm{bcc}$ is set to 0 if $\chi_0 = 7$, $\delta_\mathrm{fcc}$ if $\chi_0 = 6$, and $\delta_\mathrm{hcp}$ if $\chi_0 \le 3$.
The atom is then labelled by the first rule that applies:

1. **unknown** if $\chi_7 > 0$, that is, if two neighbors are closer than about 37° as seen from the atom.
2. **icosahedral** if $\chi_4 < 3$ and $11 \le N_1 \le 13$, and **unknown** if $\chi_4 < 3$ otherwise.
3. **bcc** if $\delta_\mathrm{bcc} \le \delta_\mathrm{cp}$ and $N_1 \ge 11$, and **unknown** if $\delta_\mathrm{bcc} \le \delta_\mathrm{cp}$ otherwise.
4. **unknown** if $N_1$ is not 11 or 12.
5. **fcc** if $\delta_\mathrm{fcc} < \delta_\mathrm{hcp}$, and **hcp** otherwise.

pyscal follows the original implementation of `compute ackland/atom` in LAMMPS.

## Usage

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.io import read

atoms = read("conf.bcc.dump", format="lammps-dump-text")
labels, names = pyscal.identify_ackland_jones(atoms)
labels[:6], names[:6]
```

`identify_ackland_jones` finds its own neighbors, so `find_neighbors` is not needed, and a neighbor list stored on `atoms` is not used and not changed.
It returns an integer label and a name for each atom.
The labels are those of `common_neighbor_analysis`:

| Label | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| Name | other | fcc | hcp | bcc | ico |

The results are stored on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_ackland_label"]` | $(N,)$ | integer label |
| `atoms.arrays["pyscal_structure"]` | $(N,)$ | integer label, the key that `common_neighbor_analysis` also uses |
| `atoms.arrays["pyscal_ackland_chi"]` | $(N, 8)$ | $\chi_0, \dots, \chi_7$ |

## Crystals at finite temperature and a liquid

How well does the method label atoms in real simulations?
We classify three MD snapshots from the `examples` folder, an fcc crystal and a bcc crystal at finite temperature, and a liquid, with the Ackland–Jones method and with adaptive CNA.

```{code-cell} ipython3
from collections import Counter

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

fractions = {}
for name, snapshot in snapshots.items():
    n = len(snapshot)
    _, names = pyscal.identify_ackland_jones(snapshot)
    fractions[name, "Ackland–Jones"] = {k: v / n for k, v in Counter(names).items()}
    counts = pyscal.common_neighbor_analysis(snapshot)
    fractions[name, "CNA"] = {("other" if k == "others" else k): v / n
                              for k, v in counts.items()}
```

`common_neighbor_analysis` calls the label 0 `others`, which is renamed here to `other`.

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

for name in ("fcc", "bcc", "liquid"):
    for structure in ("fcc", "hcp", "bcc", "ico", "other"):
        glue(f"aj_{name}_{structure}", percent(name, "Ackland–Jones", structure), display=False)
glue("cna_fcc_fcc", percent("fcc", "CNA", "fcc"), display=False)
glue("cna_bcc_bcc", percent("bcc", "CNA", "bcc"), display=False)
glue("cna_liquid_other", percent("liquid", "CNA", "other"), display=False)
assert percent("fcc", "Ackland–Jones", "fcc") > percent("fcc", "CNA", "fcc")
assert percent("bcc", "Ackland–Jones", "bcc") < percent("bcc", "CNA", "bcc")
assert percent("liquid", "Ackland–Jones", "other") > 70
```

In the fcc snapshot, the Ackland–Jones method labels {glue}`aj_fcc_fcc` % of the atoms as fcc, more than CNA with {glue}`cna_fcc_fcc` %.
The method assigns the structure with the smallest deviation, so an atom whose angles are moved slightly by thermal motion can still be labelled fcc, where CNA requires the exact signatures.
In the bcc snapshot, the method labels {glue}`aj_bcc_bcc` % of the atoms as bcc, fewer than CNA with {glue}`cna_bcc_bcc` %, and {glue}`aj_bcc_fcc` % as fcc.
In bcc, the second shell is only 15 % farther away than the first, and thermal motion moves second shell atoms across the limit of $1.45\, r_0^2$, which changes $\chi_0$ and $N_1$.
In the liquid, {glue}`aj_liquid_other` % of the atoms are unknown, and CNA labels {glue}`cna_liquid_other` % as other.

## Random displacements

To see how the labels change with the size of the thermal displacements, we add random displacements to perfect fcc, hcp and bcc crystals, with a standard deviation $\sigma$ up to 10 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations.
We compare with adaptive CNA.

```{code-cell} ipython3
from ase.build import bulk

def perfect(name):
    """Crystal and nearest neighbor distance."""
    if name == "hcp":
        return bulk("Cu", "hcp", a=2.55).repeat((8, 8, 5)), 2.55
    if name == "fcc":
        return bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6), 3.61 / np.sqrt(2)
    return bulk("Fe", "bcc", a=2.87, cubic=True).repeat(6), 2.87 * np.sqrt(3) / 2

sigmas = np.linspace(0, 0.10, 11)
scan = {}
for name in ("fcc", "hcp", "bcc"):
    crystal, nearest = perfect(name)
    for sigma in sigmas:
        noisy = crystal.copy()
        noisy.rattle(sigma * nearest, seed=1)
        _, names = pyscal.identify_ackland_jones(noisy)
        scan[name, sigma, "Ackland–Jones"] = Counter(names)
        scan[name, sigma, "CNA"] = pyscal.common_neighbor_analysis(noisy)
        scan[name, sigma, "n"] = len(noisy)
```

```{code-cell} ipython3
:tags: [hide-input]
from matplotlib.lines import Line2D

kinds = ("fcc", "hcp", "bcc", "ico", "other")
markers = dict(zip(kinds, "o^sDv"))
colour_of = lambda s: COLOURS["others"] if s == "other" else COLOURS[s]
mp = figure(columns=3, ratio=0.36, wspace=0.12)
for k, name in enumerate(("fcc", "hcp", "bcc")):
    ax = mp[0, k]
    n = scan[name, 0, "n"]
    for structure in kinds:
        values = [scan[name, s, "Ackland–Jones"].get(structure, 0) / n for s in sigmas]
        ax.plot(100 * sigmas, values, marker=markers[structure], color=colour_of(structure),
                mec=DARK, mew=0.7, ms=4.5, lw=1.5)
    cna = [scan[name, s, "CNA"][name] / n for s in sigmas]
    ax.plot(100 * sigmas, cna, ls=":", marker="o", color=COLOURS[name], mfc="white",
            mec=COLOURS[name], mew=1, ms=4.5, lw=1.5)
    label(ax, f"({'abc'[k]})  {name}")
    ax.set_ylim(-0.03, 1.08)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("Fraction of atoms")
mp[0, 1].set_xlabel(r"$\sigma$  (% of nearest neighbor distance)")
handles = [Line2D([], [], marker=markers[s], color=colour_of(s), mec=DARK, mew=0.7, ms=5,
                  lw=1.5, label=f"{s} (Ackland–Jones)") for s in kinds]
handles.append(Line2D([], [], ls=":", marker="o", color=DARK, mfc="white", mec=DARK,
                      ms=5, lw=1.5, label="correct label (CNA)"))
mp.fig.legend(handles=handles, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
def fraction(name, sigma, method, structure):
    return round(100 * scan[name, sigma, method].get(structure, 0) / scan[name, sigma, "n"])

s4, s6, s8 = sigmas[4], sigmas[6], sigmas[8]
for name in ("fcc", "hcp"):
    glue(f"aj_{name}_8", fraction(name, s8, "Ackland–Jones", name), display=False)
    glue(f"cna_{name}_8", fraction(name, s8, "CNA", name), display=False)
    # up to 6 %, the two methods agree to within a few percent
    assert all(abs(fraction(name, s, "Ackland–Jones", name) - fraction(name, s, "CNA", name)) <= 10
               for s in sigmas[:7])
    assert fraction(name, s8, "Ackland–Jones", name) > fraction(name, s8, "CNA", name)
glue("aj_bcc_4", fraction("bcc", s4, "Ackland–Jones", "bcc"), display=False)
glue("cna_bcc_4", fraction("bcc", s4, "CNA", "bcc"), display=False)
assert fraction("bcc", s4, "Ackland–Jones", "bcc") < fraction("bcc", s4, "CNA", "bcc")
```

The solid lines show the fraction of atoms that the Ackland–Jones method assigns to each label, and the dotted lines the fraction that CNA labels correctly.
For fcc and hcp, the two methods label about the same fraction of atoms correctly up to $\sigma = 6$ %.
At larger displacements the Ackland–Jones method keeps more atoms: at $\sigma = 8$ % it labels {glue}`aj_fcc_8` % of the fcc atoms and {glue}`aj_hcp_8` % of the hcp atoms correctly, against {glue}`cna_fcc_8` % and {glue}`cna_hcp_8` % with CNA.
For bcc it loses atoms much earlier: at $\sigma = 4$ % it labels {glue}`aj_bcc_4` % of the atoms as bcc, where CNA still labels {glue}`cna_bcc_4` %.

## Things to watch

- **bcc at high temperature.** Thermal motion moves atoms of the second bcc shell across the limits of $N_0$ and $N_1$, and bcc atoms are labelled fcc or unknown. Check bcc fractions with CNA.
- **Diamond structures are not covered.** The atoms of a perfect diamond crystal are labelled unknown. Use `diamond_structure` (see [Common neighbor analysis](cna.md#diamond-structures)).
- **Surfaces.** Atoms at a free surface have too few neighbors in the angle counts and are labelled unknown.
- **`pyscal_structure` is shared with CNA.** `common_neighbor_analysis` stores its labels under the same key, with the same codes. The function called last overwrites the result of the other. The Ackland–Jones labels are also in `pyscal_ackland_label`.

## References

1. G. J. Ackland and A. P. Jones, Applications of local crystal structure measures in experiment and simulation, *Phys. Rev. B* **73**, 054104 (2006). [doi:10.1103/PhysRevB.73.054104](https://doi.org/10.1103/PhysRevB.73.054104)
