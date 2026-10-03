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

# Coordination numbers

The coordination number of an atom is the number of its neighbors.
It is the simplest local descriptor, and it changes at surfaces, vacancies and interfaces.
pyscal computes four measures from the neighbor list: the coordination number, the effective coordination number of Hoppe [1], the generalized coordination number of Calle-Vallejo et al. [2], and a local density.
All four depend on the neighbor list, so read [Finding neighbors](../guide/neighbors) first.
The running coordination number obtained from $g(r)$ is described in [Distribution functions](distributions).

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

### Coordination number

The coordination number $\mathrm{CN}(i)$ of atom $i$ is the number $N(i)$ of its neighbors.
It is an integer, and it changes in steps when a neighbor crosses the cutoff.

### Effective coordination number

The effective coordination number (ECoN) of Hoppe [1] gives each neighbor a weight that decreases smoothly with its distance:

$$
\mathrm{ECoN}(i) = \sum_{j=1}^{N(i)} \exp\left[ 1 - \left( \frac{r_{ij}}{\bar{r}(i)} \right)^6 \right],
$$

where $r_{ij}$ is the distance between atom $i$ and its neighbor $j$, and $\bar{r}(i)$ is a weighted mean bond length, defined with the same weights:

$$
\bar{r}(i) = \frac{\sum_{j} r_{ij}\, w_{ij}}{\sum_{j} w_{ij}}, \qquad
w_{ij} = \exp\left[ 1 - \left( \frac{r_{ij}}{\bar{r}(i)} \right)^6 \right].
$$

pyscal solves this equation for $\bar{r}(i)$ by iteration, starting from the shortest bond of atom $i$, until $\bar{r}(i)$ changes by less than one part in $10^{12}$.

A neighbor at the distance $\bar{r}(i)$ counts 1.
A neighbor 10 % farther away counts 0.46, and one 20 % farther away counts 0.14.
Neighbors closer than $\bar{r}(i)$ count more than 1.
When all neighbors are at the same distance, ECoN equals the coordination number.

### Generalized coordination number

The generalized coordination number (GCN) of Calle-Vallejo et al. [2] counts each neighbor with its own coordination number:

$$
\mathrm{GCN}(i) = \frac{1}{\mathrm{CN}_\mathrm{max}} \sum_{j=1}^{N(i)} \mathrm{CN}(j),
$$

where $\mathrm{CN}_\mathrm{max}$ is the coordination number of an atom in the bulk, 12 for the first shell of fcc.
A neighbor that has lost neighbors itself counts less than 1.
GCN therefore distinguishes surface sites that have the same coordination number but different surroundings.
It is used as a descriptor of the adsorption energy on metal surfaces and nanoparticles [2].

### Local density

The local density of atom $i$ is its coordination number divided by the volume of a sphere with the mean bond length $\langle r \rangle_i$ as radius:

$$
\rho(i) = \frac{\mathrm{CN}(i)}{\frac{4}{3} \pi \langle r \rangle_i^3}, \qquad
\langle r \rangle_i = \frac{1}{N(i)} \sum_{j=1}^{N(i)} r_{ij}.
$$

It is larger for atoms with more or closer neighbors.

## Usage

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=3.5)

cn = pyscal.coordination_number(atoms)
econ = pyscal.effective_coordination_number(atoms)
gcn = pyscal.generalized_coordination_number(atoms, cn_max=12)
density = pyscal.local_density(atoms)
```

The cutoff of 3.5 Å is close to the first minimum of $g(r)$ of this fcc crystal (see [Finding neighbors](../guide/neighbors)).
Each function returns one value per atom and stores it on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_cn"]` | $(N,)$, integer | $\mathrm{CN}$ |
| `atoms.arrays["pyscal_econ"]` | $(N,)$ | $\mathrm{ECoN}$ |
| `atoms.arrays["pyscal_gcn"]` | $(N,)$ | $\mathrm{GCN}$ |
| `atoms.arrays["pyscal_local_density"]` | $(N,)$ | $\rho$, in atoms per Å$^3$ |

All four functions read the stored neighbor list, so `find_neighbors` must be called first.
`generalized_coordination_number` computes the coordination numbers on the way and also stores `pyscal_cn`.
Without `cn_max`, it uses the largest coordination number in the structure.

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

perfect = {}
for name, crystal in crystals.items():
    pyscal.find_neighbors(crystal, method="cutoff", cutoff=0)
    perfect[name] = {
        "CN": pyscal.coordination_number(crystal).mean(),
        "ECoN": pyscal.effective_coordination_number(crystal).mean(),
        "GCN": pyscal.generalized_coordination_number(crystal).mean(),
        "rho / (N/V)": pyscal.local_density(crystal).mean() * crystal.get_volume() / len(crystal),
    }

pd.DataFrame(perfect).T.round(2)
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
glue("econ_bcc", round(perfect["bcc"]["ECoN"], 2), display=False)
glue("rho_fcc", round(perfect["fcc"]["rho / (N/V)"], 2), display=False)
```

The adaptive cutoff takes the first shell in fcc and hcp (12 neighbors), and the first two shells in bcc (8 + 6).
In fcc and hcp all neighbors are at the same distance, and ECoN equals the coordination number.
In bcc, the 6 atoms of the second shell are 15 % farther away than the 8 of the first, and ECoN is {glue}`econ_bcc`.
With the default `cn_max`, GCN equals the coordination number in a perfect crystal, because every neighbor is itself fully coordinated.
The last column compares the local density with the number density $N/V$ of the crystal.
For fcc it is {glue}`rho_fcc` times larger, so $\rho$ is a relative measure and not the number density.

## Dependence on the cutoff

The coordination number depends on where the cutoff is placed.
We compute the mean CN and ECoN of the three MD snapshots from the `examples` folder, an fcc crystal and a bcc crystal at finite temperature and a liquid, for cutoffs from 2.8 to 5.6 Å.

```{code-cell} ipython3
snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

cutoffs = np.arange(2.8, 5.61, 0.1)
scan = {}
for name, snapshot in snapshots.items():
    values = []
    for cutoff in cutoffs:
        pyscal.find_neighbors(snapshot, method="cutoff", cutoff=cutoff)
        values.append((pyscal.coordination_number(snapshot).mean(),
                       pyscal.effective_coordination_number(snapshot).mean()))
    scan[name] = np.array(values)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.3)
for k, (title, ylabel) in enumerate([("(a)  CN", "Mean CN"), ("(b)  ECoN", "Mean ECoN")]):
    ax = mp[0, k]
    for (name, values), marker in zip(scan.items(), "osD"):
        ax.plot(cutoffs, values[:, k], marker=marker, color=COLOURS[name], mec=DARK,
                mew=0.7, ms=4, lw=1.5, label=name)
    label(ax, title)
    ax.set_xlabel("Cutoff  $r_c$  (Å)")
    ax.set_ylabel(ylabel)
for n in (12, 14):
    mp[0, 0].axhline(n, ls=":", color=DARK, lw=0.8)
mp[0, 1].set_ylim(0, 12)
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
plateau = {name: values[-1, 1] for name, values in scan.items()}
glue("econ_plateau_fcc", round(plateau["fcc"], 1), display=False)
glue("econ_plateau_bcc", round(plateau["bcc"], 1), display=False)
glue("econ_plateau_liquid", round(plateau["liquid"], 1), display=False)
```

The coordination number (a) grows with the cutoff.
In the crystals it rises slowly between the first and the second shell, near 12 for fcc and 14 for bcc, and in the liquid it rises steadily.
ECoN (b) stops changing once the cutoff passes the first shell, because distant neighbors have weights close to zero.
A cutoff that is too large does not change ECoN, and a cutoff after the first minimum of $g(r)$ is enough.

At finite temperature, ECoN is much lower than the coordination number: {glue}`econ_plateau_fcc` in fcc, {glue}`econ_plateau_bcc` in bcc and {glue}`econ_plateau_liquid` in the liquid.
Thermal motion spreads the bond lengths.
The weights favour short bonds, so $\bar{r}(i)$ is shorter than the mean bond length, and the neighbors beyond $\bar{r}(i)$ count less than 1.
Compare ECoN values only between structures at similar temperatures.

## Surface sites of a nanoparticle

On a metal nanoparticle, atoms on facets, edges and corners have different coordination numbers.
Atoms with the same coordination number can still differ in the coordination of their neighbors, for example a facet atom next to an edge and one in the middle of the facet.
GCN tells them apart.
We build a truncated octahedron of {glue}`n_particle` Pt atoms with ASE and use a fixed cutoff of 3.3 Å, between the first (2.77 Å) and the second (3.92 Å) neighbor shell.

```{code-cell} ipython3
from ase.cluster import Octahedron

particle = Octahedron("Pt", length=9, cutoff=3, latticeconstant=3.92)
pyscal.find_neighbors(particle, method="cutoff", cutoff=3.3)
cn = pyscal.coordination_number(particle)
gcn = pyscal.generalized_coordination_number(particle, cn_max=12)
```

The particle has no cell, so `find_neighbors` treats it as an isolated cluster.
For comparison, we compute the coordination numbers of the top layer of flat (111) and (100) surfaces of Pt.

```{code-cell} ipython3
from ase.build import fcc100, fcc111

surfaces = {}
for name, build in [("(111)", fcc111), ("(100)", fcc100)]:
    slab = build("Pt", size=(4, 4, 6), a=3.92, vacuum=8, periodic=True)
    pyscal.find_neighbors(slab, method="cutoff", cutoff=3.3)
    top = slab.positions[:, 2] > slab.positions[:, 2].max() - 0.1
    surfaces[name] = (pyscal.coordination_number(slab)[top][0],
                      pyscal.generalized_coordination_number(slab, cn_max=12)[top][0])
surfaces
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.48, wspace=0.25, width_ratios=[1, 1.3])

# (a) surface atoms on the front half of the particle, seen along [111]
ax = mp[0, 0]
view = np.array([1, 1, 1]) / np.sqrt(3)
across = np.array([1, -1, 0]) / np.sqrt(2)
up = np.cross(view, across)
positions = particle.positions - particle.positions.mean(axis=0)
front = np.flatnonzero((cn < 12) & (positions @ view > -2))
front = front[np.argsort(positions[front] @ view)]
points = ax.scatter(positions[front] @ across, positions[front] @ up, c=gcn[front],
                    cmap="viridis", vmin=4, vmax=8, s=38, ec=DARK, lw=0.4)
ax.set_aspect("equal")
ax.set_axis_off()
label(ax, "(a)  surface, seen along [111]")
colourbar = mp.fig.colorbar(points, ax=ax, orientation="horizontal", fraction=0.06,
                            pad=0.04, aspect=25)
colourbar.set_label("GCN")

# (b) GCN against CN, marker area proportional to the number of atoms
ax = mp[0, 1]
pairs, counts = np.unique(np.column_stack([cn, gcn.round(3)]), axis=0, return_counts=True)
ax.scatter(pairs[:, 0], pairs[:, 1], s=10 + 2.5 * counts, color=COLOURS["fcc"], ec=DARK,
           lw=0.7, zorder=3)
sites = {6: "corners", 7: "edges", 8: "(100) facets", 9: "(111) facets",
         12: "inside"}
for c, text in sites.items():
    g = pairs[pairs[:, 0] == c, 1]
    left = c == 8
    ax.annotate(text, (c, g.min()), xytext=(-8, 6) if left else (7, -7),
                textcoords="offset points", ha="right" if left else "left",
                va="bottom" if left else "top", fontsize=9, color=DARK)
ax.set_xlabel("CN")
ax.set_ylabel("GCN")
ax.set_xlim(5, 13.6)
ax.set_ylim(2.5, 12.8)
label(ax, "(b)  all atoms")
note(ax, "marker area: number of atoms");
```

```{code-cell} ipython3
:tags: [remove-cell]
facet = cn == 9
glue("n_particle", len(particle), display=False)
glue("gcn_facet_min", round(gcn[facet].min(), 2), display=False)
glue("gcn_facet_max", round(gcn[facet].max(), 2), display=False)
glue("gcn_111", round(surfaces["(111)"][1], 2), display=False)
glue("gcn_100", round(surfaces["(100)"][1], 2), display=False)
glue("gcn_bulk_min", round(gcn[cn == 12].min(), 2), display=False)
```

Panel (a) shows the surface atoms on the front half of the particle, seen along [111].
A hexagonal (111) facet faces the viewer, and the other facets are seen at an angle.
Atoms in the middle of a facet have the highest GCN of the surface, and atoms at edges and corners the lowest.

Panel (b) shows that GCN splits each coordination number into several values.
The atoms with CN = 9 sit on (111) facets.
Their GCN is {glue}`gcn_facet_max` when all their neighbors are facet atoms or atoms inside the particle, and drops to {glue}`gcn_facet_min` next to a corner.
The flat Pt(111) surface has GCN = {glue}`gcn_111` and Pt(100) has {glue}`gcn_100`, as in [2].
Atoms right below the surface have 12 neighbors like bulk atoms, but GCN as low as {glue}`gcn_bulk_min`, because some of their neighbors are surface atoms.

```{code-cell} ipython3
:tags: [remove-cell]
pyscal.find_neighbors(atoms, method="cutoff", cutoff=3.5)
glue("cn_max_fcc", int(pyscal.coordination_number(atoms).max()), display=False)
glue("gcn_default_fcc", round(pyscal.generalized_coordination_number(atoms).mean(), 2), display=False)
glue("gcn_twelve_fcc", round(pyscal.generalized_coordination_number(atoms, cn_max=12).mean(), 2), display=False)
```

## Things to watch

- **The neighbor list sets the values.** All four measures count the atoms in the neighbor list. The adaptive cutoff gives 14 neighbors in perfect bcc, a fixed cutoff between the shells gives 8. Compare values only for the same neighbor method.
- **ECoN at finite temperature.** Thermal motion spreads the bond lengths, and ECoN drops well below the coordination number, as shown above.
- **The default of `cn_max`.** Without `cn_max`, `generalized_coordination_number` divides by the largest coordination number in the structure. At finite temperature some atoms have more neighbors than in the perfect crystal, and this default is larger than the bulk value. In the fcc snapshot with a cutoff of 3.5 Å, the largest coordination number is {glue}`cn_max_fcc`, and the mean GCN without `cn_max` is {glue}`gcn_default_fcc` instead of {glue}`gcn_twelve_fcc`. Pass `cn_max` explicitly.
- **The local density is not $N/V$.** The sphere of radius $\langle r \rangle_i$ is smaller than the volume per atom times $\mathrm{CN}(i)$, so $\rho$ is about twice the number density in close packed crystals. Use it to compare atoms within one structure, for example to find compressed or expanded regions.

## References

1. R. Hoppe, Effective coordination numbers (ECoN) and mean fictive ionic radii (MEFIR), *Z. Kristallogr.* **150**, 23 (1979). [doi:10.1524/zkri.1979.150.14.23](https://doi.org/10.1524/zkri.1979.150.14.23)
2. F. Calle-Vallejo, J. Tymoczko, V. Colic, Q. H. Vu, M. D. Pohl, K. Morgenstern, D. Loffreda, P. Sautet, W. Schuhmann and A. S. Bandarenka, Finding optimal surface sites on heterogeneous catalysts by counting nearest neighbors, *Science* **350**, 185 (2015). [doi:10.1126/science.aab3501](https://doi.org/10.1126/science.aab3501)
