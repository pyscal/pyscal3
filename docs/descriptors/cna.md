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

# Common neighbor analysis

Common neighbor analysis (CNA) [1, 2] assigns a crystal structure to each atom from the way its neighbors are bonded to each other.
It labels atoms as fcc, hcp, bcc, icosahedral or other, and, in its extended form, as cubic or hexagonal diamond.
CNA is the usual choice for finding defects in a crystal: stacking faults, grain boundaries, dislocation cores and surfaces show up as atoms that are not labelled with the structure of the crystal.

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

CNA looks at each bond between an atom $i$ and one of its neighbors $j$ and records three numbers:

1. the number of neighbors that $i$ and $j$ have in common,
2. the number of bonds between these common neighbors,
3. the number of bonds in the longest chain formed by these bonds.

Each bond gets a signature of three digits, for example 421.
The signatures of the bonds of an atom identify its structure:

| Structure | Neighbors | Signatures |
|---|---|---|
| fcc | 12 | 12 × 421 |
| hcp | 12 | 6 × 421, 6 × 422 |
| bcc | 14 | 8 × 666, 6 × 444 |
| icosahedral | 12 | 12 × 555 |

An atom whose signatures match none of these is labelled *others*.

### Adaptive and conventional CNA

The neighbors in CNA are the atoms within a cutoff $r_c$ between the first and the second neighbor shell, or, for bcc, between the second and the third.
**Conventional CNA** computes this cutoff from a given lattice constant $a$:
$r_c = \tfrac{1}{2}(1/\sqrt{2} + 1)\,a \approx 0.854\,a$ for fcc and hcp, and $r_c = \tfrac{1}{2}(1 + \sqrt{2})\,a \approx 1.207\,a$ for bcc.

**Adaptive CNA** [3] computes the cutoff for each atom from its own nearest atoms.
For fcc and hcp it uses the 12 nearest atoms $j$, at distances $r_{ij}$ sorted in increasing order:

$$
r_c^\mathrm{fcc}(i) = \frac{1 + \sqrt{2}}{2} \cdot \frac{1}{12} \sum_{j=1}^{12} r_{ij}.
$$

For bcc it uses the 14 nearest atoms:

$$
r_c^\mathrm{bcc}(i) = \frac{1 + \sqrt{2}}{2} \cdot \frac{1}{14} \left( \frac{2}{\sqrt{3}} \sum_{j=1}^{8} r_{ij} + \sum_{j=9}^{14} r_{ij} \right).
$$

An atom is first tested for fcc, hcp and icosahedral with $r_c^\mathrm{fcc}$, and then for bcc with $r_c^\mathrm{bcc}$.
Adaptive CNA needs no lattice constant, and it works when the lattice constant changes from place to place, for example under strain, at interfaces between two phases, or with thermal expansion.

## Usage

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.bcc.dump", format="lammps-dump-text")
pyscal.common_neighbor_analysis(atoms)
```

Without `lattice_constant`, `common_neighbor_analysis` uses adaptive CNA.
With `lattice_constant=a`, it uses conventional CNA.
It returns the number of atoms of each structure and stores the label of each atom in `atoms.arrays["pyscal_structure"]`:

| Label | 0 | 1 | 2 | 3 | 4 |
|---|---|---|---|---|---|
| Structure | others | fcc | hcp | bcc | icosahedral |

`common_neighbor_analysis` finds its own neighbors, so `find_neighbors` is not needed, and a neighbor list stored on `atoms` is left unchanged.

## Strain and temperature

To compare adaptive and conventional CNA, we strain perfect fcc and bcc crystals by a linear strain $\varepsilon$ and add random displacements with a standard deviation of 4 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations.
Conventional CNA gets the lattice constant of the unstrained crystal.

```{code-cell} ipython3
import numpy as np
from ase.build import bulk

lattice_constants = {"fcc": 3.61, "bcc": 2.87}
nearest = {"fcc": 3.61 / np.sqrt(2), "bcc": 2.87 * np.sqrt(3) / 2}
strains = np.linspace(-0.16, 0.24, 21)

def fraction_labelled(name, strain, noise, adaptive):
    a = lattice_constants[name]
    crystal = bulk("Cu", name, a=a * (1 + strain), cubic=True).repeat(7)
    crystal.rattle(noise * nearest[name] * (1 + strain), seed=2)
    counts = pyscal.common_neighbor_analysis(
        crystal, lattice_constant=None if adaptive else a)
    return counts[name] / len(crystal)

strain_scan = {
    (name, adaptive): [fraction_labelled(name, e, 0.04, adaptive) for e in strains]
    for name in lattice_constants for adaptive in (True, False)
}
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.12)
for k, name in enumerate(lattice_constants):
    ax = mp[0, k]
    ax.axvline(0, ls=":", color=DARK, lw=0.8)
    for adaptive, ls, marker, text in [(True, "-", "o", "adaptive CNA"),
                                       (False, ":", "s", "conventional CNA")]:
        ax.plot(100 * strains, strain_scan[name, adaptive], ls=ls, marker=marker,
                color=COLOURS[name], mfc=COLOURS[name] if adaptive else "white",
                mec=DARK if adaptive else COLOURS[name], mew=0.9, ms=5, lw=1.6,
                label=text)
    label(ax, f"({'ab'[k]})  {name}")
    ax.set_xlabel(r"Strain  $\varepsilon$  (%)")
    ax.set_ylim(-0.03, 1.12)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel(f"Fraction labelled correctly")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.02));
```

Adaptive CNA labels every atom correctly over the whole range.
Conventional CNA loses atoms once the strain exceeds a few percent, because the thermal displacements move atoms across its fixed cutoff.
Use conventional CNA only when the lattice constant is known and uniform.

Large thermal displacements limit both variants.
Below, we label perfect fcc, bcc and hcp crystals with increasing random displacements, using adaptive CNA.

```{code-cell} ipython3
noises = np.linspace(0, 0.14, 15)

def noisy(name, noise):
    if name == "hcp":
        crystal, d = bulk("Cu", "hcp", a=2.55).repeat((10, 10, 6)), 2.55
    else:
        crystal, d = bulk("Cu", name, a=lattice_constants[name], cubic=True).repeat(7), nearest[name]
    crystal.rattle(noise * d, seed=1)
    return pyscal.common_neighbor_analysis(crystal)[name] / len(crystal)

noise_scan = {name: [noisy(name, s) for s in noises] for name in ("fcc", "hcp", "bcc")}
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(ratio=0.42)
ax = mp[0, 0]
for (name, values), marker in zip(noise_scan.items(), "osD"):
    ax.plot(100 * noises, values, marker=marker, color=COLOURS[name], mec=DARK,
            mew=0.7, ms=5, lw=1.6, label=name)
ax.set_xlabel("Standard deviation of the displacements  (% of nearest neighbor distance)")
ax.set_ylabel("Fraction labelled correctly")
ax.legend(frameon=False);
```

Above a standard deviation of about 5 % of the nearest neighbor distance, the fraction of correctly labelled atoms drops quickly.
By the Lindemann criterion, the root mean square displacement at melting is roughly 10 to 15 % of the nearest neighbor distance, or 6 to 9 % in each direction.
A crystal close to its melting point therefore has many atoms labelled *others*, as in [A first analysis](../tour).
At high temperature, the [averaged Steinhardt parameters](steinhardt) still separate crystal from liquid atoms where CNA labels many crystal atoms as *others*.

## Example: stacking faults in fcc

An intrinsic stacking fault in fcc is a missing close-packed layer.
The stacking sequence ABCABC becomes ABCACABC, and the two layers next to the fault have hcp surroundings.
We build an fcc crystal with two such faults by removing two (111) layers.

```{code-cell} ipython3
from ase.build import fcc111

perfect = fcc111("Cu", size=(6, 6, 21), a=3.61, orthogonal=True, periodic=True)
spacing = perfect.cell[2, 2] / 21
layer = np.round(perfect.positions[:, 2] / spacing).astype(int)

removed = [6, 15]
keep = ~np.isin(layer, removed)
faulted = perfect[keep]
layer = layer[keep]

# close the gaps left by the removed layers
below = np.array([sum(l > r for r in removed) for l in layer])
faulted.positions[:, 2] -= below * spacing
faulted.cell[2, 2] -= len(removed) * spacing
layer -= below

pyscal.common_neighbor_analysis(faulted)
```

```{code-cell} ipython3
:tags: [hide-input]
names = {0: "others", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico"}
labels = faulted.arrays["pyscal_structure"]
mp = figure(ratio=0.3)
ax = mp[0, 0]
for code in np.unique(labels):
    sel = labels == code
    ax.scatter(faulted.positions[sel, 2], faulted.positions[sel, 1], s=22,
               color=COLOURS[names[code]], ec=DARK, lw=0.4, label=names[code])
ax.set_xlabel("$z$  along [111]  (Å)")
ax.set_ylabel("$y$  (Å)")
ax.set_aspect("equal")
ax.legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.28));
```

The figure shows the atoms projected on the $y$–$z$ plane, with $z$ along [111] and the close-packed layers vertical.
Each stacking fault appears as two layers of hcp atoms in the fcc crystal.

## Diamond structures

Cubic and hexagonal diamond have only four nearest neighbors, and the signatures above do not apply.
`diamond_structure` follows Maras et al. [4]: the 12 second neighbors of an atom (the neighbors of its neighbors) are tested with the fcc and hcp signatures.
An fcc arrangement means cubic diamond, an hcp arrangement hexagonal diamond.

```{code-cell} ipython3
a = 5.43
cubic = bulk("Si", "diamond", a=a, cubic=True).repeat(4)
hexagonal = bulk("SiSi", "wurtzite", a=a / np.sqrt(2), c=a / np.sqrt(2) * np.sqrt(8 / 3)).repeat((6, 6, 4))

{name: {k: v for k, v in pyscal.diamond_structure(s).items() if v}
 for name, s in [("cubic", cubic), ("hexagonal", hexagonal)]}
```

The labels in `atoms.arrays["pyscal_structure"]` are:

| Label | Structure |
|---|---|
| 0 | others |
| 1 | cubic diamond |
| 2 | first neighbor of a cubic diamond atom |
| 3 | second neighbor of a cubic diamond atom |
| 4 | hexagonal diamond |
| 5 | first neighbor of a hexagonal diamond atom |
| 6 | second neighbor of a hexagonal diamond atom |

Labels 2, 3, 5 and 6 mark atoms that are not diamond themselves but bonded to diamond atoms, such as atoms at a surface or a defect.

## Things to watch

- **Thermal displacements.** At high temperature, many atoms of a crystal are labelled *others*. Use averaged Steinhardt parameters, or average the labels over time, when this matters.
- **Surfaces.** Atoms at a free surface have too few neighbors and are labelled *others*.
- **Few atoms.** CNA needs 14 candidate neighbors for every atom (4 for `diamond_structure`). In very small or very dilute systems, it warns and labels atoms with too few candidates as *others*.

## References

1. J. D. Honeycutt and H. C. Andersen, Molecular dynamics study of melting and freezing of small Lennard-Jones clusters, *J. Phys. Chem.* **91**, 4950 (1987). [doi:10.1021/j100303a014](https://doi.org/10.1021/j100303a014)
2. D. Faken and H. Jónsson, Systematic analysis of local atomic structure combined with 3D computer graphics, *Comput. Mater. Sci.* **2**, 279 (1994). [doi:10.1016/0927-0256(94)90109-0](https://doi.org/10.1016/0927-0256(94)90109-0)
3. A. Stukowski, Structure identification methods for atomistic simulations of crystalline materials, *Modelling Simul. Mater. Sci. Eng.* **20**, 045021 (2012). [doi:10.1088/0965-0393/20/4/045021](https://doi.org/10.1088/0965-0393/20/4/045021)
4. E. Maras, O. Trushin, A. Stukowski, T. Ala-Nissila and H. Jónsson, Global transition path search for dislocation formation in Ge on Si(001), *Comput. Phys. Commun.* **205**, 13 (2016). [doi:10.1016/j.cpc.2016.04.001](https://doi.org/10.1016/j.cpc.2016.04.001)
