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

# Wigner–Seitz defect analysis

Wigner–Seitz analysis finds point defects by comparing a structure with a perfect reference lattice [1].
Each atom is assigned to the nearest lattice site, and the number of atoms on each site gives the vacancies, the interstitials and, in alloys, the antisites.
It is used to count the defects left by collision cascades in radiation damage simulations and to follow vacancies and interstitials as they diffuse.
[Common neighbor analysis](cna) and the [centrosymmetry parameter](centrosymmetry.md) label atoms by the structure around them, and thermal vibrations at high temperature make many atoms look defective.
Wigner–Seitz analysis counts atoms per site instead, so it gives the number of point defects also at high temperature.
The [deformation descriptors](deformation) also use a reference structure, but compare the neighbors of each atom.

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

The reference structure has $M$ lattice sites at positions $\mathbf{R}_s$, $s = 1, \dots, M$.
The Wigner–Seitz cell of a site is the region of space closer to this site than to any other.
Each atom $i$ of the analysed structure, at position $\mathbf{r}_i$, is assigned to the site whose Wigner–Seitz cell contains it:

$$
s(i) = \arg\min_{s} \left| \mathbf{r}_i - \mathbf{R}_s \right|,
$$

using the nearest periodic image.
The occupancy $n_s$ of site $s$ is the number of atoms assigned to it.
An empty site ($n_s = 0$) is a vacancy, and every atom beyond the first on a site ($n_s \ge 2$) is an interstitial:

$$
N_\mathrm{vac} = \left| \{ s : n_s = 0 \} \right|, \qquad
N_\mathrm{int} = \sum_{s=1}^{M} \max(n_s - 1, 0).
$$

Since every atom is assigned to exactly one site, $N_\mathrm{int} - N_\mathrm{vac} = N - M$, where $N$ is the number of atoms.
In an alloy, a site with one atom whose chemical species differs from that of the site in the reference is an antisite.

The analysis has no parameters.
It needs only the reference lattice, which must have the same cell, orientation and origin as the analysed structure.

## Usage

We remove one atom from a perfect fcc crystal:

```{code-cell} ipython3
import pyscal
from ase.build import bulk

reference = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6)
atoms = reference.copy()
del atoms[100]

result = pyscal.wigner_seitz_analysis(atoms, reference)
result["vacancy_count"], result["interstitial_count"], result["vacancy_indices"]
```

The analysis finds one vacancy, at site 100 of the reference.
`wigner_seitz_analysis` returns a dictionary:

| Key | Shape | Content |
|---|---|---|
| `"vacancy_count"` | integer | $N_\mathrm{vac}$ |
| `"interstitial_count"` | integer | $N_\mathrm{int}$ |
| `"occupancy"` | $(M,)$ | $n_s$ of each reference site |
| `"site_index"` | $(N,)$ | $s(i)$ of each atom |
| `"vacancy_indices"` | $(N_\mathrm{vac},)$ | the empty sites |
| `"interstitial_sites"` | array | the sites with $n_s \ge 2$ |
| `"occupancy_by_type"` | dictionary of $(M,)$ | $n_s$ for each chemical species, with `per_type_occupancies=True` |

It also stores results on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_ws_site_index"]` | $(N,)$ | $s(i)$ |
| `atoms.arrays["pyscal_ws_occupancy"]` | $(N,)$ | $n_{s(i)}$, the occupancy of the site of each atom |
| `atoms.info["pyscal_ws_vacancy_count"]` | integer | $N_\mathrm{vac}$ |
| `atoms.info["pyscal_ws_interstitial_count"]` | integer | $N_\mathrm{int}$ |

`identify_defect_atoms` calls `wigner_seitz_analysis` with `per_type_occupancies=True` and adds per atom masks and a summary:

```{code-cell} ipython3
defects = pyscal.identify_defect_atoms(atoms, reference)
defects["defect_summary"], defects["vacancy_positions"]
```

| Key | Shape | Content |
|---|---|---|
| `"perfect_mask"` | $(N,)$ | `True` for atoms alone on their site |
| `"interstitial_mask"` | $(N,)$ | `True` for atoms on a site with $n_s \ge 2$ |
| `"vacancy_positions"` | $(N_\mathrm{vac}, 3)$ | positions of the empty sites in the reference |
| `"defect_summary"` | string | the numbers of vacancies, interstitials and antisites |

The dictionary also contains all keys returned by `wigner_seitz_analysis`.

Neither function uses neighbors, so `find_neighbors` is not needed.
The structure and the reference can have different numbers of atoms, and the atoms can be in any order.
Atoms outside the cell are wrapped back into it along the periodic directions of the reference.

## Point defects at high temperature

The `examples` folder contains an MD snapshot of an fcc crystal at high temperature.
As the reference, we build a perfect fcc lattice with the same cell as the snapshot.

```{code-cell} ipython3
import numpy as np
from ase.io import read

snapshot = read("conf.fcc.dump", format="lammps-dump-text")
a = snapshot.cell.lengths()[1] / 6
lattice = bulk(snapshot.get_chemical_symbols()[0], "fcc", a=a, cubic=True).repeat((7, 6, 6))
lattice.set_cell(snapshot.cell, scale_atoms=True)

cna = pyscal.common_neighbor_analysis(snapshot.copy())
result = pyscal.wigner_seitz_analysis(snapshot, lattice)
cna["others"], result["vacancy_count"], result["interstitial_count"]
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
from ase.geometry import find_mic
displacement, _ = find_mic(snapshot.positions - lattice.positions[result["site_index"]], snapshot.cell)
nearest = a / np.sqrt(2)
sigma_snapshot = displacement.std(axis=0).mean()
glue("n_snapshot", len(snapshot), display=False)
glue("cna_others", cna["others"], display=False)
glue("cna_others_pct", round(100 * cna["others"] / len(snapshot)), display=False)
glue("sigma_snapshot", round(sigma_snapshot, 2), display=False)
glue("sigma_snapshot_pct", round(100 * sigma_snapshot / nearest, 1), display=False)
assert result["vacancy_count"] == 0 and result["interstitial_count"] == 0
```

CNA labels {glue}`cna_others` of the {glue}`n_snapshot` atoms ({glue}`cna_others_pct` %) as *others*, because of the large thermal displacements.
The standard deviation of the displacements from the lattice sites is {glue}`sigma_snapshot` Å in each direction, {glue}`sigma_snapshot_pct` % of the nearest neighbor distance.
The Wigner–Seitz analysis assigns exactly one atom to every site and finds no defect.

We now add defects to the snapshot.
Two atoms are removed, which leaves two vacancies.
Two more atoms are taken from their sites and inserted as split interstitials at two other sites.
At each of these sites, the inserted atom and the atom already there are displaced by $\pm 0.3\,a$ along $x$ from the site.
This gives two Frenkel pairs, so the structure has four vacancies and two interstitials.
All six sites lie in a slab of the crystal, which is shown in the figure below.

```{code-cell} ipython3
rng = np.random.default_rng(4)
z = lattice.positions[:, 2]
in_slab = (z > 0.25 * a) & (z < 1.25 * a)
chosen = rng.choice(np.where(in_slab)[0], size=6, replace=False)
removed_sites, emptied_sites, split_sites = chosen[:2], chosen[2:4], chosen[4:]

on_site = {s: np.where(result["site_index"] == s)[0][0] for s in chosen}
defected = snapshot.copy()
shift = np.array([0.3 * a, 0.0, 0.0])
for empty, split in zip(emptied_sites, split_sites):
    centre = defected.positions[on_site[split]].copy()
    defected.positions[on_site[split]] = centre - shift
    defected.positions[on_site[empty]] = centre + shift
del defected[[on_site[s] for s in removed_sites]]

defects = pyscal.identify_defect_atoms(defected, lattice)
defects["defect_summary"]
```

```{code-cell} ipython3
:tags: [remove-cell]
assert set(defects["vacancy_indices"]) == set(removed_sites) | set(emptied_sites)
assert set(defects["interstitial_sites"]) == set(split_sites)
glue("n_interstitial_mask", int(defects["interstitial_mask"].sum()), display=False)
glue("n_interstitial", defects["interstitial_count"], display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
names = {0: "others", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico"}
pyscal.common_neighbor_analysis(defected)
structure = defected.arrays["pyscal_structure"]
zd = defected.positions[:, 2]
slab = (zd > 0.25 * a - 0.5) & (zd < 1.25 * a + 0.5)
xy = defected.positions[:, :2]

mp = figure(columns=2, ratio=0.5, wspace=0.08)
ax = mp[0, 0]
for code in (1, 0):
    sel = slab & (structure == code)
    ax.scatter(*xy[sel].T, s=14, color=COLOURS[names[code]], ec=DARK, lw=0.3,
               label=names[code])
label(ax, "(a)  CNA")

ax = mp[0, 1]
ordinary = slab & defects["perfect_mask"]
ax.scatter(*xy[ordinary].T, s=14, color=COLOURS["liquid"], ec=DARK, lw=0.3, alpha=0.35,
           label="one atom on its site")
ax.scatter(*xy[slab & defects["interstitial_mask"]].T, s=34, marker="D", color=COLOURS["bcc"],
           ec=DARK, lw=0.8, label="atoms on a site with two atoms", zorder=3)
ax.scatter(*defects["vacancy_positions"][:, :2].T, s=70, marker="s", facecolor="none",
           ec=COLOURS["diamond"], lw=1.6, label="vacancy", zorder=3)
label(ax, "(b)  Wigner–Seitz")
ax.tick_params(labelleft=False)
for ax in (mp[0, 0], mp[0, 1]):
    ax.set_aspect("equal")
    ax.set_xlim(-1.5, snapshot.cell[0, 0] + 1.5)
    ax.set_ylim(-1.5, snapshot.cell[1, 1] + 1.5)
    ax.set_xlabel("$x$  (Å)")
mp[0, 0].set_ylabel("$y$  (Å)")
mp[0, 0].legend(frameon=False, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.2),
                markerscale=1.5)
mp[0, 1].legend(frameon=False, ncol=1, loc="upper center", bbox_to_anchor=(0.5, -0.2));
```

The figure shows the atoms in the slab, projected on the $x$–$y$ plane.
In panel (a), CNA labels atoms as *others* throughout the slab, and the defects cannot be told apart from the thermal disorder.
In panel (b), the Wigner–Seitz analysis finds the four empty sites (squares) and the two sites with two atoms (diamonds).
`interstitial_mask` marks both atoms of a split interstitial, so it is `True` for {glue}`n_interstitial_mask` atoms while there are {glue}`n_interstitial` interstitials.

## How large can the displacements be?

An atom is assigned to the wrong site when it moves outside the Wigner–Seitz cell of its own site, about half the nearest neighbor distance $d$ away.
This creates a spurious vacancy and a spurious interstitial.
To find out when this happens, we displace the atoms of a perfect fcc crystal randomly, with a standard deviation $\sigma$ in each direction, and count the vacancies.

```{code-cell} ipython3
d = 3.61 / np.sqrt(2)
perfect = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(8)
fractions = np.arange(0, 0.25, 0.02)
spurious_noise = []
for f in fractions:
    rattled = perfect.copy()
    rattled.rattle(f * d, seed=1)
    spurious_noise.append(pyscal.wigner_seitz_analysis(rattled, perfect)["vacancy_count"])
```

A second source of spurious defects is a change of the cell, for example by thermal expansion in a simulation at constant pressure.
We expand the crystal uniformly by a linear strain $\varepsilon$, add random displacements with $\sigma = 0.05\,d$, and analyse it with and without `affine_mapping="to_reference"`.
This option maps the positions to the cell of the reference before the analysis:
$\mathbf{r}_i' = \mathbf{H}_\mathrm{ref} \mathbf{H}^{-1} \mathbf{r}_i$, where the columns of $\mathbf{H}$ and $\mathbf{H}_\mathrm{ref}$ are the cell vectors of the structure and of the reference.

```{code-cell} ipython3
strains = np.arange(0, 0.0501, 0.005)
spurious_strain = {}
for repeat in (6, 12):
    perfect = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(repeat)
    for mapping in ("none", "to_reference"):
        counts = []
        for e in strains:
            expanded = perfect.copy()
            expanded.set_cell(perfect.cell * (1 + e), scale_atoms=True)
            expanded.rattle(0.05 * d, seed=1)
            r = pyscal.wigner_seitz_analysis(expanded, perfect, affine_mapping=mapping)
            counts.append(r["vacancy_count"] / len(perfect))
        spurious_strain[repeat, mapping] = counts
```

```{code-cell} ipython3
:tags: [remove-cell]
from _plotstyle import RED, TEAL, PURPLE
n_noise = len(bulk("Cu", "fcc", cubic=True).repeat(8))
noise_fraction = np.array(spurious_noise) / n_noise
first = fractions[np.argmax(noise_fraction > 0)]
glue("n_noise", n_noise, display=False)
glue("d_nearest", round(d, 2), display=False)
glue("noise_last_clean_pct", round(100 * fractions[np.argmax(noise_fraction > 0) - 1]), display=False)
glue("noise_first_pct", round(100 * first), display=False)
glue("noise_max_pct", round(100 * fractions[-1]), display=False)
glue("noise_max_vac_pct", round(100 * noise_fraction[-1], 1), display=False)
for repeat in (6, 12):
    counts = np.array(spurious_strain[repeat, "none"])
    e_first = strains[np.argmax(counts > 0)]
    length = 3.61 * repeat
    glue(f"L_{repeat}", round(length, 1), display=False)
    glue(f"e_first_{repeat}", round(100 * e_first, 1), display=False)
    glue(f"shift_first_{repeat}", round(e_first * length, 2), display=False)
    assert max(spurious_strain[repeat, "to_reference"]) == 0
assert all(c == 0 for c in spurious_noise[:np.argmax(noise_fraction > 0)])
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.3)
ax = mp[0, 0]
ax.plot(100 * fractions, 100 * noise_fraction, marker="o", color=COLOURS["fcc"], mec=DARK,
        mew=0.7, ms=5, lw=1.6)
ax.axvline(100 * sigma_snapshot / nearest, ls="--", color=DARK, lw=1)
ax.annotate("fcc snapshot", (100 * sigma_snapshot / nearest, 10), xytext=(5, 0),
            textcoords="offset points", ha="left", fontsize=9, color=DARK)
label(ax, "(a)  random displacements")
ax.set_xlabel(r"$\sigma / d$  (%)")
ax.set_ylabel("Spurious vacancies  (% of sites)")
ax.set_ylim(-0.5, 15)

ax = mp[0, 1]
styles = {(6, "none"): (RED, "o", "-", "6 × 6 × 6 cells, no mapping"),
          (12, "none"): (PURPLE, "s", "-", "12 × 12 × 12 cells, no mapping"),
          (12, "to_reference"): (TEAL, "D", ":", "12 × 12 × 12 cells, to_reference")}
for key, (colour, marker, ls, text) in styles.items():
    ax.plot(100 * strains, 100 * np.array(spurious_strain[key]), marker=marker, ls=ls,
            color=colour, mec=DARK, mew=0.7, ms=5, lw=1.6, label=text)
label(ax, "(b)  uniform expansion")
ax.set_xlabel(r"Linear strain  $\varepsilon$  (%)")
ax.set_ylabel("Spurious vacancies  (% of sites)")
ax.set_ylim(-0.5, 15)
mp.fig.legend(*ax.get_legend_handles_labels(), frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

Panel (a) shows that the analysis finds no spurious defect in the {glue}`n_noise` atoms up to $\sigma$ = {glue}`noise_last_clean_pct` % of $d$.
The first spurious defects appear at {glue}`noise_first_pct` %, and at {glue}`noise_max_pct` % of $d$ {glue}`noise_max_vac_pct` % of the sites are wrongly counted as vacant.
The dashed line marks the fcc snapshot above.
By the Lindemann criterion, crystals melt when the root mean square displacement is roughly 10 to 15 % of $d$, or 6 to 9 % in each direction.
The thermal displacements in a crystal below its melting point therefore stay in the range where this test finds no spurious defects.
The number of spurious defects at a given $\sigma$ grows with the number of atoms.

Panel (b) shows the effect of a uniform expansion.
Without mapping, an atom at a distance $x$ from the origin of the cell is shifted by $\varepsilon x$ from its site.
In the crystal with an edge length $L$ = {glue}`L_6` Å, spurious defects appear at $\varepsilon$ = {glue}`e_first_6` %, and in the crystal with $L$ = {glue}`L_12` Å already at {glue}`e_first_12` %.
At these strains, the largest shift $\varepsilon L$ is {glue}`shift_first_6` Å and {glue}`shift_first_12` Å, about a third of $d$ = {glue}`d_nearest` Å.
With `affine_mapping="to_reference"`, no spurious defect appears for either size, and the line for 6 × 6 × 6 cells (not shown) is also zero.

## Antisites in an ordered alloy

In an L1₂ ordered alloy such as Cu₃Au, the Au atoms occupy the corners of the cubic cell and the Cu atoms the face centres.
We swap three Au atoms with one of their Cu neighbors each, which gives six antisites, and add random displacements with a standard deviation of 0.1 Å as a rough model of thermal vibrations.

```{code-cell} ipython3
unit = bulk("Cu", "fcc", a=3.75, cubic=True)
unit.symbols[0] = "Au"
ordered = unit.repeat(6)

alloy = ordered.copy()
rng = np.random.default_rng(1)
for i in rng.choice(np.where(alloy.symbols == "Au")[0], size=3, replace=False):
    distances = alloy.get_distances(i, range(len(alloy)), mic=True)
    j = next(j for j in np.argsort(distances) if alloy.symbols[j] == "Cu")
    alloy.symbols[i], alloy.symbols[j] = "Cu", "Au"
alloy.rattle(0.1, seed=2)

defects = pyscal.identify_defect_atoms(alloy, ordered)
defects["defect_summary"]
```

`identify_defect_atoms` reports the number of antisites only in `defect_summary`.
The antisite atoms are those whose species differs from the species of their site in the reference:

```{code-cell} ipython3
species = np.array(alloy.get_chemical_symbols())
site_species = np.array(ordered.get_chemical_symbols())[defects["site_index"]]
antisite = species != site_species
np.where(antisite)[0], species[antisite]
```

The occupancies of each species, in `defects["occupancy_by_type"]`, show which sites hold a Cu atom and which an Au atom.

## Things to watch

- **The reference must match.** The reference lattice must have the same cell, orientation and origin as the structure. If the whole crystal drifts during a simulation, for example because the total momentum was not set to zero, subtract the drift of the centre of mass before the analysis.
- **Changes of the cell.** Use `affine_mapping="to_reference"` when the cell of the structure differs from that of the reference. The mapping removes uniform strain only. A crystal that deforms inhomogeneously can still give spurious defects.
- **Large displacements.** An atom that moves about half the nearest neighbor distance away from its site is assigned to a neighboring site, and a spurious vacancy and interstitial appear. In the test above, this happens from $\sigma$ = {glue}`noise_first_pct` % of $d$. The analysis is meaningless for liquids or amorphous regions.
- **Interstitial atoms.** `interstitial_mask` marks every atom on a site with two or more atoms, not only the extra ones. A split interstitial gives two marked atoms.
- **Antisites.** `identify_defect_atoms` counts antisites only on sites with exactly one atom, and only when the structure contains more than one species. It does not return an antisite mask, but the mask follows from `site_index` as shown above.

## References

1. K. Nordlund, M. Ghaly, R. S. Averback, M. Caturla, T. Diaz de la Rubia and J. Tarus, Defect production in collision cascades in elemental semiconductors and fcc metals, *Phys. Rev. B* **57**, 7556 (1998). [doi:10.1103/PhysRevB.57.7556](https://doi.org/10.1103/PhysRevB.57.7556)
