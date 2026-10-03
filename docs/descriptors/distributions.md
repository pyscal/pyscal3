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

# Distribution functions

Distribution functions describe a whole structure rather than single atoms.
The radial distribution function $g(r)$ is the density of atoms at a distance $r$ from an atom, relative to the mean density.
The angular distribution function (ADF) is the distribution of the angles between the bonds of an atom, and the bond length distribution function (BLDF) the distribution of the bond lengths in the neighbor list.
They are used to compare structures with experiment and with each other, to choose a cutoff for the neighbors, and to check which bonds a neighbor method includes.
[Finding neighbors](../guide/neighbors) uses $g(r)$ to choose a cutoff.
This page explains the normalisation of $g(r)$ and shows the ADF and the BLDF of crystals and a liquid.

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

### Radial distribution function

For $N$ atoms in a volume $V$, with number density $\rho = N/V$, the radial distribution function is

$$
g(r) = \frac{1}{N \rho} \sum_{i=1}^{N} \sum_{j \neq i} \frac{\delta(r - r_{ij})}{4 \pi r^2},
$$

where $r_{ij}$ is the distance between atoms $i$ and $j$, using the nearest periodic image, and $\delta$ is the Dirac delta function.
pyscal computes $g(r)$ as a histogram with bins of width $\Delta r$.
For the bin that starts at $r_k$,

$$
g(r_k) = \frac{n_k}{N \rho\, \Delta V_k}, \qquad
\Delta V_k = \frac{4}{3} \pi \left[ (r_k + \Delta r)^3 - r_k^3 \right],
$$

where $n_k$ is the number of pairs $(i, j)$ with $r_k \le r_{ij} < r_k + \Delta r$, each pair counted once for $i$ and once for $j$.
With this normalisation, $g(r) = 1$ for an ideal gas, and $g(r) \to 1$ at large $r$ in a liquid.
The running coordination number

$$
n(r) = 4 \pi \rho \int_0^r g(r') \, r'^2 \, dr'
$$

is the mean number of atoms within the distance $r$ of an atom [1].

### Angular distribution function

Each pair of neighbors $j$ and $k$ of an atom $i$ forms a bond angle $\theta_{jik}$:

$$
\cos \theta_{jik} = \frac{\mathbf{r}_{ij} \cdot \mathbf{r}_{ik}}{r_{ij}\, r_{ik}},
$$

where $\mathbf{r}_{ij}$ is the vector from atom $i$ to its neighbor $j$.
An atom with $N(i)$ neighbors has $N(i)(N(i) - 1)/2$ bond angles.
The ADF $P(\theta)$ is the histogram of the bond angles of all atoms between 0° and 180°, normalised so that $\int_0^{180°} P(\theta)\, d\theta = 1$, with $\theta$ in degrees.

If the bond directions were uncorrelated, $P(\theta)$ would be proportional to $\sin \theta$:

$$
P_\mathrm{random}(\theta) = \frac{\pi}{360} \sin \theta.
$$

### Bond length distribution function

The BLDF $P(r)$ is the histogram of the distances $r_{ij}$ of all bonds in the neighbor list, normalised so that $\int P(r)\, dr = 1$.
Unlike $g(r)$, it is not divided by the volume of a shell, and it contains only the bonds in the neighbor list.

## Usage

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.lqd.Al.dump", format="lammps-dump-text")
g, r = pyscal.radial_distribution_function(atoms, rmin=0, rmax=8.0, bins=160)

pyscal.find_neighbors(atoms, method="cutoff", cutoff=4.1)
adf, angles = pyscal.angular_distribution_function(atoms, bins=180)
bldf, lengths = pyscal.bond_length_distribution(atoms, bins=100)
```

All three functions return the histogram first and the **left edges** of the bins second.

- `radial_distribution_function(atoms, rmin, rmax, bins)` returns $g(r)$ and the left edges $r_k$, in Å. It finds its own neighbors with a fixed cutoff of `rmax`, and **replaces the neighbor list stored on `atoms`**. This is why `find_neighbors` is called after it above. It does not store $g(r)$ on `atoms`. It needs a cell that is periodic in all three directions, because it uses $\rho = N/V$.
- `angular_distribution_function(atoms, bins)` returns $P(\theta)$, in 1/degree, and the left edges in degrees. The bins divide the range from 0° to 180°.
- `bond_length_distribution(atoms, bins, rmin, rmax)` returns $P(r)$, in 1/Å, and the left edges in Å. Without `rmin` and `rmax`, the range runs from 0.9 times the shortest to 1.1 times the longest bond.

The ADF and the BLDF use the stored neighbor list, so `find_neighbors` must be called first.
They store their results in `atoms.info`:

| Key | Shape | Content |
|---|---|---|
| `atoms.info["pyscal_adf"]` | (`bins`,) | $P(\theta)$ |
| `atoms.info["pyscal_adf_angles"]` | (`bins`,) | left bin edges, in degrees |
| `atoms.info["pyscal_bldf"]` | (`bins`,) | $P(r)$ |
| `atoms.info["pyscal_bldf_r"]` | (`bins`,) | left bin edges, in Å |

`angular_distribution_function` computes the bond angles with `chi_params`, so it also stores the cosines of the bond angles of each atom in `atoms.info["pyscal_cosines"]` and the $\chi$ parameters of the [Ackland–Jones method](ackland_jones) in `atoms.arrays["pyscal_chiparams"]`.

## Normalisation of g(r)

We compare $g(r)$ of the liquid with an ideal gas: 4000 atoms at random positions in the same cell.
Then we compute the running coordination number $n(r)$ of the three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and the liquid.

```{code-cell} ipython3
import numpy as np
from ase import Atoms

liquid = read("conf.lqd.Al.dump", format="lammps-dump-text")
rng = np.random.default_rng(1)
gas = Atoms("Al4000", cell=liquid.cell, pbc=True)
gas.positions = rng.random((len(gas), 3)) @ gas.cell.array

pair = {}
for name, structure in [("liquid", liquid), ("ideal gas", gas)]:
    g, r = pyscal.radial_distribution_function(structure, rmax=10.0, bins=100)
    pair[name] = (r + 0.5 * (r[1] - r[0]), g)          # bin centres
```

```{code-cell} ipython3
snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": liquid,
}

running = {}
first_minimum = {}
for name, snapshot in snapshots.items():
    g, r = pyscal.radial_distribution_function(snapshot, rmax=6.0, bins=120)
    dr = r[1] - r[0]
    rho = len(snapshot) / snapshot.get_volume()
    n = rho * np.cumsum(g * 4 / 3 * np.pi * ((r + dr) ** 3 - r ** 3))
    peak = np.argmax(g)
    minimum = peak + np.argmin(g[peak:peak + 40])
    running[name] = (r + dr, n)                         # n(r) at the right bin edges
    first_minimum[name] = (r[minimum] + dr, n[minimum])
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.3)

ax = mp[0, 0]
ax.axhline(1, ls=":", color=DARK, lw=0.8)
r, g = pair["ideal gas"]
ax.plot(r, g, ls="--", color=DARK, lw=1.2, label="ideal gas")
r, g = pair["liquid"]
ax.plot(r, g, color=COLOURS["liquid"], lw=1.8, label="liquid")
label(ax, "(a)  $g(r)$")
ax.set_xlabel("$r$  (Å)")
ax.set_ylabel("$g(r)$")
ax.set_xlim(0, 10)
ax.set_ylim(0, 3.2)
ax.legend(frameon=False, loc="upper right")

ax = mp[0, 1]
for (name, (r, n)), ls in zip(running.items(), ["-", "--", "-"]):
    ax.plot(r, n, color=COLOURS[name], lw=1.8, ls=ls, label=name)
    rmin, nmin = first_minimum[name]
    ax.plot(rmin, nmin, marker="o", ms=6, color=COLOURS[name], mec=DARK, mew=0.8,
            zorder=4)
for value in (12, 14):
    ax.axhline(value, ls=":", color=DARK, lw=0.8)
label(ax, "(b)  $n(r)$")
ax.set_xlabel("$r$  (Å)")
ax.set_ylabel("$n(r)$")
ax.set_xlim(2, 6)
ax.set_ylim(0, 60)
ax.legend(frameon=False, loc="upper left")
note(ax, "dots: first minimum of $g(r)$", loc="lower right");
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

far = pair["liquid"][0] > 6
glue("g_far_liquid", round(pair["liquid"][1][far].mean(), 2), display=False)
for name in snapshots:
    glue(f"n_min_{name}", round(first_minimum[name][1], 1), display=False)
    glue(f"r_min_{name}", round(first_minimum[name][0], 2), display=False)

# n(r) is the mean number of neighbors within r
check = snapshots["liquid"].copy()
pyscal.find_neighbors(check, method="cutoff", cutoff=first_minimum["liquid"][0])
glue("cn_check_liquid", round(np.diff(check.info["pyscal_bond_offsets"]).mean(), 1), display=False)
```

In panel (a), $g(r)$ of the ideal gas is 1 at all distances, apart from statistical noise.
The noise is largest at small $r$, where the shells contain few atoms.
$g(r)$ of the liquid is zero below about 2 Å, where atoms do not overlap, has peaks at the neighbor shells, and tends to 1 at large $r$: its mean beyond 6 Å is {glue}`g_far_liquid`.

Panel (b) shows $n(r)$.
At the first minimum of $g(r)$ (dots), it is {glue}`n_min_fcc` in fcc, {glue}`n_min_bcc` in bcc and {glue}`n_min_liquid` in the liquid.
These are the numbers of neighbors with a fixed cutoff at that distance.
For the liquid, `find_neighbors` with this cutoff gives {glue}`cn_check_liquid` neighbors per atom.
In the crystals, $n(r)$ is nearly flat around the first minimum, at 12 in fcc and 14 in bcc.
In the liquid, it rises steadily, and the number of neighbors depends on where exactly the cutoff is placed.

## Bond angles

We find the neighbors of each snapshot with a fixed cutoff at the first minimum of $g(r)$, as above, and compute the ADF.

```{code-cell} ipython3
adfs = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=first_minimum[name][0])
    adf, angles = pyscal.angular_distribution_function(snapshot, bins=90)
    adfs[name] = (angles + 1, adf)                      # bin centres
```

```{code-cell} ipython3
:tags: [hide-input]
ideal = {
    "fcc": [60, 90, 120, 180],
    "bcc": [54.74, 70.53, 90, 109.47, 125.26, 180],
}
theta = np.linspace(0, 180, 181)
mp = figure(columns=3, ratio=0.33, wspace=0.12)
for k, (name, (angles, adf)) in enumerate(adfs.items()):
    ax = mp[0, k]
    for angle in ideal.get(name, []):
        ax.axvline(angle, ls=":", color=DARK, lw=0.8)
    ax.plot(theta, np.pi / 360 * np.sin(np.radians(theta)), ls="--", color=DARK, lw=1.2,
            label=r"uncorrelated, $\propto \sin\theta$")
    ax.plot(angles, adf, color=COLOURS[name], lw=1.8, label=name)
    label(ax, f"({'abc'[k]})  {name}")
    ax.set_xlabel(r"$\theta$  (degree)")
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 60, 120, 180])
    ax.set_ylim(0, 0.028)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel(r"$P(\theta)$  (1/degree)")
handles, labels = mp[0, 2].get_legend_handles_labels()
mp.fig.legend(handles[:1], labels[:1], frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
from scipy.signal import find_peaks

def peaks(name):
    angles, adf = adfs[name]
    found, _ = find_peaks(adf, prominence=0.001)
    return angles[found]

glue("liquid_peak_1", int(peaks("liquid")[0]), display=False)
glue("liquid_peak_2", int(peaks("liquid")[1]), display=False)
glue("fcc_last_peak", int(peaks("fcc")[-1]), display=False)
```

The dotted lines mark the angles in the perfect crystals: 60°, 90°, 120° and 180° between the 12 neighbors in fcc, and 54.7°, 70.5°, 90°, 109.5°, 125.3° and 180° between the 14 neighbors in bcc.
At finite temperature the peaks are broad.
In fcc, they lie near the angles of the perfect crystal.
In bcc, the peaks at 54.7°, 90° and 125.3° are clear, and those at 70.5° and 109.5°, between the bonds of the first shell, are only shoulders.
In the liquid there are two broad peaks, near {glue}`liquid_peak_1`° and {glue}`liquid_peak_2`°.
The dashed line is the ADF of uncorrelated bond directions.
It is zero at 0° and 180°, because few directions lie close to a given axis.
The same factor $\sin \theta$ moves the peak of the crystals at 180° to {glue}`fcc_last_peak`° in fcc.

## Bond lengths and the neighbor method

The BLDF shows which bonds a neighbor method includes.
We compute it for the bcc crystal and the liquid with the four methods compared in [Finding neighbors](../guide/neighbors).
We multiply each BLDF by the mean number of neighbors, so that the area under each curve is the number of neighbors per atom.

```{code-cell} ipython3
methods = {
    "fixed cutoff": lambda name: dict(method="cutoff", cutoff=first_minimum[name][0]),
    "adaptive": lambda name: dict(method="cutoff", cutoff=0),
    "SANN": lambda name: dict(method="cutoff", cutoff="sann"),
    "Voronoi": lambda name: dict(method="voronoi"),
}

bldfs = {}
for name in ("bcc", "liquid"):
    snapshot = snapshots[name]
    for method, options in methods.items():
        pyscal.find_neighbors(snapshot, **options(name))
        bldf, lengths = pyscal.bond_length_distribution(snapshot, bins=90, rmin=2.0, rmax=6.5)
        neighbors = len(snapshot.info["pyscal_bond_neighbors"]) / len(snapshot)
        longest = snapshot.info["pyscal_bond_distance"].max()
        bldfs[name, method] = (lengths + 0.025, bldf * neighbors, neighbors, longest)
```

```{code-cell} ipython3
:tags: [hide-input]
from _plotstyle import METHODS

mp = figure(columns=2, ratio=0.42, wspace=0.12)
for k, name in enumerate(("bcc", "liquid")):
    ax = mp[0, k]
    # the methods with more neighbors are drawn first, so that every band shows
    # the bonds that one method adds to the next
    for method in ("Voronoi", "fixed cutoff", "SANN", "adaptive"):
        lengths, values = bldfs[name, method][:2]
        ax.fill_between(lengths, values, step="mid", color=METHODS[method][0], ec=DARK,
                        lw=0.6, label=method)
    label(ax, f"({'ab'[k]})  {name}")
    ax.set_xlabel("Bond length  $r_{ij}$  (Å)")
    ax.set_xlim(2, 6.5)
    ax.set_ylim(0, 22)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel("Neighbors per atom per Å")
a = (2 * snapshots["bcc"].get_volume() / len(snapshots["bcc"])) ** (1 / 3)
for shell in (np.sqrt(3) / 2 * a, a):
    mp[0, 0].axvline(shell, ls=":", color=DARK, lw=0.8)
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles[::-1], labels[::-1], frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
glue("voronoi_longest", round(bldfs["liquid", "Voronoi"][3], 1), display=False)
for name in ("bcc", "liquid"):
    for method, key in [("adaptive", "adaptive"), ("Voronoi", "voronoi"), ("fixed cutoff", "fixed")]:
        glue(f"nn_{name}_{key}", round(bldfs[name, method][2], 1), display=False)
```

In bcc (a), the dotted lines mark the first shell, at $\sqrt{3}a/2$, and the second shell, at $a$, of the perfect crystal with the same density.
At finite temperature the two shells overlap.
The fixed cutoff includes both shells, {glue}`nn_bcc_fixed` neighbors per atom.
The adaptive cutoff stops inside the second shell, with {glue}`nn_bcc_adaptive` neighbors.
The Voronoi neighbors are nearly the same as with the fixed cutoff, plus a few longer bonds from the small faces created by thermal motion.

In the liquid (b), the methods differ more.
The adaptive cutoff gives {glue}`nn_liquid_adaptive` neighbors per atom, the fixed cutoff {glue}`nn_liquid_fixed` and the Voronoi method {glue}`nn_liquid_voronoi`, with bonds up to {glue}`voronoi_longest` Å.
All methods agree on the short bonds, and they differ in where they cut the tail of the first shell.
Descriptors computed from these neighbor lists differ in the same way.

## Things to watch

- **`radial_distribution_function` replaces the neighbor list.** It calls `find_neighbors` with a fixed cutoff of `rmax`. Compute $g(r)$ first, or on a copy (`atoms.copy()`), when the neighbor list is needed afterwards.
- **The bins are given by their left edges.** Plot $g(r)$, the ADF and the BLDF at the bin centres, $r + \Delta r / 2$, as done above. Plotting at the left edges shifts the curves by half a bin.
- **Total, not partial, $g(r)$.** `radial_distribution_function` counts all pairs, regardless of the species. For partial $g_{AB}(r)$ of an alloy, select the pairs from the neighbor list yourself.
- **The ADF is not divided by $\sin \theta$.** Uncorrelated bonds give a distribution proportional to $\sin \theta$, so angles near 90° are more frequent than angles near 0° or 180° for purely geometric reasons. Divide by $\sin \theta$ to compare the frequency of angles directly.
- **`rmin` and `rmax` of the BLDF.** Bonds outside the range are left out, and the rest is normalised to 1. Choose a range that contains all bonds when comparing neighbor methods, as above.

## References

1. M. P. Allen and D. J. Tildesley, *Computer Simulation of Liquids*, 2nd ed. (Oxford University Press, Oxford, 2017).
