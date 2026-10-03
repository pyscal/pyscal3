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

# Centrosymmetry parameter

The centrosymmetry parameter of Kelchner, Plimpton and Hamilton [1] measures how far the neighbors of an atom are from inversion symmetry.
In a perfect fcc or bcc crystal, every neighbor of an atom has a partner on the opposite side, and the parameter is zero.
Defects such as surfaces, stacking faults, dislocation cores and vacancies break this symmetry and give positive values.
The parameter is used to find defects in fcc and bcc metals.
It gives a continuous value rather than a structure label, unlike [common neighbor analysis](cna).
For defects defined relative to a reference structure, see [Atomic deformation](deformation) and [Wigner–Seitz analysis](wigner_seitz).

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

For an atom $i$, take its $N$ nearest neighbors, where $N$ is even.
For each of the $N(N - 1)/2$ pairs of neighbors $j$ and $k$, compute

$$
w_{jk} = \left| \mathbf{r}_{ij} + \mathbf{r}_{ik} \right|^2,
$$

where $\mathbf{r}_{ij}$ is the vector from atom $i$ to its neighbor $j$.
$w_{jk} = 0$ when $j$ and $k$ are on exactly opposite sides of $i$ at the same distance.
The centrosymmetry parameter is the sum of the $N/2$ smallest values of $w_{jk}$:

$$
P(i) = \sum_{m=1}^{N/2} w_{(m)},
$$

where $w_{(1)} \le w_{(2)} \le \dots$ are the values of $w_{jk}$ sorted in increasing order.
$P$ has the unit of a length squared, Å² for positions in Å.

The selected pairs are not required to be disjoint, so one neighbor can belong to more than one of the $N/2$ pairs.
This is the scheme used in LAMMPS, called greedy edge selection by Larsen [2], who also discusses a variant in which every neighbor belongs to exactly one pair.
pyscal implements the greedy scheme.

$N$ is usually the number of first shell neighbors of the crystal: 12 for fcc and 8 for bcc.
For bcc, $N = 14$ includes the second shell, which also has inversion symmetry.

## Usage

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
csp = pyscal.centrosymmetry(atoms, nmax=12)
csp[:5].round(2)
```

`centrosymmetry` returns $P$ for each atom, with `nmax` $= N$ (default 12).
`nmax` must be a positive even number, otherwise the function raises a `ValueError`.
The result is stored in `atoms.arrays["pyscal_centrosymmetry"]`, with shape $(N_\mathrm{atoms},)$.

`centrosymmetry` finds its own neighbors, so `find_neighbors` is not needed.
It calls `find_neighbors(atoms, method="number", nmax=nmax)`, which replaces the neighbor list stored on `atoms`:

```{code-cell} ipython3
np.unique(np.diff(atoms.info["pyscal_bond_offsets"]))
```

Every atom now has the 12 neighbors used for $P$.
Compute other descriptors before calling `centrosymmetry`, or call it on a copy (`atoms.copy()`).

## Stacking faults, vacancies and surfaces

Which values does $P$ take at common defects in fcc?
We build an fcc crystal with close packed (111) layers, as in [Common neighbor analysis](cna), with three kinds of defects:
a stacking fault, made by removing one (111) layer and closing the gap,
a vacancy, made by removing one atom,
and two free surfaces, made by switching off the periodic boundary conditions along [111].

```{code-cell} ipython3
from ase.build import fcc111
from ase.geometry import get_distances

a = 3.61
crystal = fcc111("Cu", size=(6, 6, 18), a=a, orthogonal=True, periodic=True)
spacing = crystal.cell[2, 2] / 18
layer = np.round(crystal.positions[:, 2] / spacing).astype(int)

# stacking fault: remove layer 6 and close the gap
keep = layer != 6
crystal, layer = crystal[keep], layer[keep]
crystal.positions[layer > 6, 2] -= spacing

# free surfaces at the bottom and the top
crystal.pbc = [True, True, False]

# vacancy: remove the atom of layer 12 closest to the centre of the layer
centre = crystal.cell.diagonal()[:2] / 2
in_layer = np.flatnonzero(layer == 12)
vacancy = in_layer[np.argmin(np.linalg.norm(crystal.positions[in_layer, :2] - centre, axis=1))]
_, distance = get_distances(crystal.positions[vacancy], crystal.positions,
                            cell=crystal.cell, pbc=crystal.pbc)
site = np.full(len(crystal), "bulk", dtype=object)
site[np.isin(layer, [5, 7])] = "stacking fault"
site[np.isin(layer, [layer.min(), layer.max()])] = "surface"
site[distance[0] < 0.85 * a] = "vacancy neighbor"
del crystal[vacancy]
site = np.delete(site, vacancy)

perfect_csp = pyscal.centrosymmetry(crystal, nmax=12)
```

`site` records which defect each atom belongs to: the two layers next to the stacking fault, the outermost layer on each side, and the 12 neighbors of the vacancy.
We also add random displacements with a standard deviation of 4 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations.

```{code-cell} ipython3
vibrating = crystal.copy()
vibrating.rattle(0.04 * a / np.sqrt(2), seed=3)
vibrating_csp = pyscal.centrosymmetry(vibrating, nmax=12)
```

```{code-cell} ipython3
:tags: [hide-input]
from _plotstyle import RED

colours = {"bulk": COLOURS["fcc"], "stacking fault": COLOURS["hcp"],
           "vacancy neighbor": RED, "surface": COLOURS["others"]}
markers = {"bulk": "o", "stacking fault": "^", "vacancy neighbor": "s", "surface": "D"}

mp = figure(columns=2, ratio=0.42, wspace=0.12)
for k, (structure, values, title) in enumerate([
        (crystal, perfect_csp, "(a)  perfect positions"),
        (vibrating, vibrating_csp, "(b)  random displacements")]):
    ax = mp[0, k]
    for name in colours:
        sel = site == name
        ax.scatter(structure.positions[sel, 2], values[sel], s=16, marker=markers[name],
                   color=colours[name], ec=DARK, lw=0.4, label=name, zorder=3)
    label(ax, title)
    ax.set_xlabel("$z$  along [111]  (Å)")
    ax.set_ylim(-1, 22)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("$P$  (Å$^2$)")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.1), markerscale=1.5);
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

def values(csp, name):
    return csp[site == name]

for name in ("stacking fault", "vacancy neighbor", "surface"):
    v = values(perfect_csp, name)
    assert np.ptp(v) < 1e-8
    glue(f"perfect_{name.split()[0]}", round(float(v[0]), 2), display=False)
assert np.all(values(perfect_csp, "bulk") < 1e-8)
glue("a_squared_half", round(a**2 / 2, 2), display=False)
glue("bulk_max", round(float(values(vibrating_csp, "bulk").max()), 1), display=False)
glue("defect_min", round(float(min(values(vibrating_csp, n).min()
                                   for n in ("stacking fault", "vacancy neighbor", "surface"))), 1),
     display=False)
```

With perfect positions (a), $P$ is zero in the bulk.
The atoms next to the stacking fault and the neighbors of the vacancy both have $P =$ {glue}`perfect_stacking` Å², which is $a^2/2$ for the lattice constant $a = 3.61$ Å.
The surface atoms have $P =$ {glue}`perfect_surface` Å².
With random displacements (b), the bulk atoms reach up to {glue}`bulk_max` Å², and the lowest value of a defect atom is {glue}`defect_min` Å².
A threshold between these two values finds all defect atoms in this crystal, but the margin is small.
$P$ does not tell a stacking fault from a vacancy, because both have the same value.
To label the stacking fault, use [common neighbor analysis](cna).

## Thermal vibrations in MD snapshots

At higher temperature, the thermal values of $P$ approach those of defects.
We compute $P$ for the fcc and bcc MD snapshots from the `examples` folder and compare it with the value of a vacancy neighbor in a perfect crystal with the same lattice constant.
For bcc, we compare $N = 8$ and $N = 14$.

```{code-cell} ipython3
from ase.build import bulk

def vacancy_value(structure, a, nmax):
    """P of the neighbors of a vacancy in a perfect crystal."""
    perfect = bulk("Cu", structure, a=a, cubic=True).repeat(4)
    del perfect[0]
    return pyscal.centrosymmetry(perfect, nmax=nmax).max()

fcc = read("conf.fcc.dump", format="lammps-dump-text")
bcc = read("conf.bcc.dump", format="lammps-dump-text")
a_fcc = (4 * fcc.get_volume() / len(fcc)) ** (1 / 3)
a_bcc = (2 * bcc.get_volume() / len(bcc)) ** (1 / 3)

snapshot_csp = {
    ("fcc", 12): pyscal.centrosymmetry(fcc, nmax=12),
    ("bcc", 8): pyscal.centrosymmetry(bcc, nmax=8),
    ("bcc", 14): pyscal.centrosymmetry(bcc, nmax=14),
}
reference = {
    ("fcc", 12): vacancy_value("fcc", a_fcc, 12),
    ("bcc", 8): vacancy_value("bcc", a_bcc, 8),
    ("bcc", 14): vacancy_value("bcc", a_bcc, 14),
}
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.12)
bins = np.linspace(0, 20, 51)
panels = [[("fcc", 12)], [("bcc", 8), ("bcc", 14)]]
styles = {("fcc", 12): (COLOURS["fcc"], "-"), ("bcc", 8): (COLOURS["bcc"], ":"),
          ("bcc", 14): (COLOURS["bcc"], "-")}
for k, keys in enumerate(panels):
    ax = mp[0, k]
    for key in keys:
        colour, ls = styles[key]
        ax.hist(snapshot_csp[key], bins=bins, density=True, histtype="step", color=colour,
                ls=ls, lw=1.8, label=f"$N = {key[1]}$")
    ax.axvline(reference[keys[0]], ls="--", color=DARK, lw=1)
    note(ax, f"vacancy neighbor\n{reference[keys[0]]:.1f} Å$^2$", loc="upper right")
    ax.set_xlabel("$P$  (Å$^2$)")
    ax.set_xlim(0, 20)
    ax.set_ylim(0, 0.42)
    ax.legend(frameon=False, loc="center right")
    if k:
        ax.tick_params(labelleft=False)
label(mp[0, 0], "(a)  fcc snapshot")
label(mp[0, 1], "(b)  bcc snapshot")
mp[0, 0].set_ylabel("Probability density  (Å$^{-2}$)");
```

```{code-cell} ipython3
:tags: [remove-cell]
assert np.isclose(reference["bcc", 8], reference["bcc", 14])
for (name, nmax), csp in snapshot_csp.items():
    glue(f"above_{name}_{nmax}", round(100 * float(np.mean(csp > reference[name, nmax]))),
         display=False)
    glue(f"median_{name}_{nmax}", round(float(np.median(csp)), 1), display=False)
glue("a_fcc", round(float(a_fcc), 2), display=False)

# atoms whose 8 nearest neighbors include an atom of the second shell (along <100>)
nearest8 = bcc.copy()
pyscal.find_neighbors(nearest8, method="number", nmax=8)
offsets = nearest8.info["pyscal_bond_offsets"]
second = np.abs(nearest8.info["pyscal_bond_vector"]).max(axis=1) > 0.75 * a_bcc
glue("second_shell_8", round(100 * float(np.mean(np.add.reduceat(second, offsets[:-1]) > 0))),
     display=False)
glue("a_bcc", round(float(a_bcc), 2), display=False)
```

The dashed lines mark the value of a vacancy neighbor.
The lattice constants of the perfect crystals, {glue}`a_fcc` Å for fcc and {glue}`a_bcc` Å for bcc, are computed from the volume per atom of the snapshots.
In the fcc snapshot (a), the median of $P$ is {glue}`median_fcc_12` Å², and {glue}`above_fcc_12` % of the atoms lie above the value of a vacancy neighbor, although the snapshot has as many atoms as lattice sites.
A threshold for defects has to be placed above most of the thermal distribution, and some defect atoms can then be missed.

In the bcc snapshot (b), the choice of $N$ matters.
With $N = 8$, the median is {glue}`median_bcc_8` Å², and {glue}`above_bcc_8` % of the atoms lie above the value of a vacancy neighbor.
With $N = 14$, the median is {glue}`median_bcc_14` Å², and {glue}`above_bcc_14` % of the atoms lie above it.
The second shell in bcc is only 15 % farther away than the first.
At this temperature, the 8 nearest atoms of {glue}`second_shell_8` % of the atoms include at least one atom of the second shell, which has no partner among the 8.
The separate peaks of the distribution for $N = 8$ belong to atoms with none, one and two such atoms among their 8 nearest neighbors.
With $N = 14$, both shells are included, and the pairs stay complete.

```{code-cell} ipython3
:tags: [remove-cell]
# check: P = d^2 in perfect hcp
from pyscal.structures import make_crystal
d = 2.51
hcp = make_crystal("hcp", lattice_constant=d, repetitions=(4, 4, 4))
assert np.allclose(pyscal.centrosymmetry(hcp, nmax=12), d**2)
```

## Things to watch

- **The neighbor list is replaced.** `centrosymmetry` stores its own list of the `nmax` nearest neighbors. Compute descriptors that need another neighbor list first, or call `centrosymmetry` on a copy.
- **Choose `nmax` for the crystal.** Use 12 for fcc. For bcc, use 14 at finite temperature, as shown above, or 8 at low temperature. Values computed with different `nmax` cannot be compared.
- **hcp is not centrosymmetric.** In a perfect hcp crystal, every atom has $P = d^2$ with $N = 12$, where $d$ is the nearest neighbor distance. $P$ is therefore not suited for finding defects in hcp crystals.
- **Values scale with the lattice constant.** $P$ grows with the square of the interatomic distances. To compare materials, divide by $d^2$.

## References

1. C. L. Kelchner, S. J. Plimpton and J. C. Hamilton, Dislocation nucleation and defect structure during surface indentation, *Phys. Rev. B* **58**, 11085 (1998). [doi:10.1103/PhysRevB.58.11085](https://doi.org/10.1103/PhysRevB.58.11085)
2. P. M. Larsen, Revisiting the common neighbour analysis and the centrosymmetry parameter, arXiv:2003.08879 (2020). [arXiv:2003.08879](https://arxiv.org/abs/2003.08879)
