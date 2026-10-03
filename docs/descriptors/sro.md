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

# Chemical short range order

The Warren–Cowley parameter [1, 2] measures whether the atoms of an alloy have more or fewer neighbors of another species than in a random mixture.
It is zero for a random solid solution, negative when unlike neighbors are preferred, as in ordered compounds, and positive when like atoms cluster.
It is used to follow ordering and clustering in simulations of alloys, for example high entropy alloys, and it is measured in experiments by diffuse X-ray or neutron scattering.
Like every descriptor here, it depends on the neighbor list, see [Finding neighbors](../guide/neighbors).

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

For an atom $i$ of species A, let $N(i)$ be its number of neighbors and $N_B(i)$ the number of its neighbors of species B.
The Warren–Cowley parameter of atom $i$ is

$$
\alpha_{AB}(i) = 1 - \frac{p_{AB}(i)}{c_B}, \qquad p_{AB}(i) = \frac{N_B(i)}{N(i)},
$$

where $c_B$ is the fraction of atoms of species B in the whole structure.
pyscal averages it over the $N_A$ atoms of species A:

$$
\bar{\alpha}_{AB} = \frac{1}{N_A} \sum_{i \in A} \alpha_{AB}(i).
$$

- $\bar{\alpha}_{AB} = 0$: the A atoms have as many B neighbors as in a random mixture.
- $\bar{\alpha}_{AB} < 0$: they have more B neighbors, the alloy orders.
- $\bar{\alpha}_{AB} > 0$: they have fewer B neighbors, like atoms cluster.

The value of one atom lies between $1 - 1/c_B$, when all neighbors are B, and 1, when none is.
When all atoms have the same number of neighbors, $\bar{\alpha}_{AB} = \bar{\alpha}_{BA}$, and $\bar{\alpha}_{AB}$ is the pair parameter of Cowley [1].
If B is the minority species and all atoms have the same number of neighbors, the most negative value is $-c_B/c_A$, reached when B atoms have only A neighbors.

Cowley defines one parameter $\alpha_n$ for each neighbor shell $n$ of the lattice [1].
In pyscal, the shell is set by the neighbor list: a cutoff between the first and the second shell gives $\alpha_1$, and the `shell_thickness` option of `find_neighbors` selects a single outer shell.

### Like pairs

With B = A, the same formula gives $\alpha_{AA}(i) = 1 - p_{AA}(i)/c_A$.
Its sign has the opposite meaning: $\alpha_{AA} > 0$ means fewer like neighbors than in a random mixture, so ordering.
In a binary alloy in which all atoms have the same number of neighbors, $\bar{\alpha}_{AA} = -(c_B/c_A)\, \bar{\alpha}_{AB}$, and the two carry the same information.

## Usage

```{code-cell} ipython3
import pyscal
from pyscal.structures import make_crystal

# Cu3Au in the L1_2 structure: Au on the cube corners, Cu on the face centres
atoms = make_crystal("l12", lattice_constant=3.75, repetitions=(6, 6, 6),
                     element=["Au", "Cu"])
pyscal.find_neighbors(atoms, method="cutoff", cutoff=3.2)

alpha = pyscal.short_range_order(atoms, reference_type="Cu", compare_type="Au")
per_atom = pyscal.short_range_order(atoms, "Cu", "Au", average=False)
alpha
```

The cutoff of 3.2 Å lies between the first (2.65 Å) and the second (3.75 Å) neighbor shell, so this is $\alpha_1$.
`short_range_order` returns $\bar{\alpha}_{AB}$, or with `average=False` the value of each atom.
It stores the values of the atoms in `atoms.arrays`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_sro"]` | $(N,)$ | $\alpha_{AB}(i)$ for atoms of species A, `NaN` for all other atoms and for atoms without neighbors |

The species are given as chemical symbols or atomic numbers.
Without `reference_type`, A is the most abundant species, and without `compare_type`, B is the most abundant of the others.
For equal amounts, the species with the lower atomic number comes first.
`short_range_order` reads the stored neighbor list, so `find_neighbors` must be called first.

### Values for ordered alloys

```{code-cell} ipython3
import numpy as np
import pandas as pd

nial = make_crystal("b2", lattice_constant=2.88, repetitions=(6, 6, 6), element=["Ni", "Al"])
values = {}

for first, second in [("Cu", "Au"), ("Au", "Cu"), ("Cu", "Cu"), ("Au", "Au")]:
    values["Cu3Au, L1_2, 12 neighbors", f"alpha_{first}{second}"] = \
        pyscal.short_range_order(atoms, first, second)

for name, options in [("8 neighbors", dict(cutoff=2.7)), ("14 neighbors", dict(cutoff=0))]:
    pyscal.find_neighbors(nial, method="cutoff", **options)
    for first, second in [("Ni", "Al"), ("Ni", "Ni")]:
        values[f"NiAl, B2, {name}", f"alpha_{first}{second}"] = \
            pyscal.short_range_order(nial, first, second)

pd.Series(values).round(3)
```

In Cu$_3$Au, each Cu atom has 4 Au and 8 Cu neighbors, so $p_\mathrm{CuAu} = 1/3$, $c_\mathrm{Au} = 1/4$ and $\alpha_\mathrm{CuAu} = -1/3$.
This is the most negative value possible at this composition, because each Au atom has only Cu neighbors.
In B2 NiAl, the 8 nearest neighbors of a Ni atom are all Al, and $\alpha_\mathrm{NiAl} = -1$.
The adaptive cutoff (`cutoff=0`) also includes the 6 Ni atoms of the second shell, and $\alpha_\mathrm{NiAl}$ becomes $1 - (8/14)/(1/2) = -1/7$.

## Long range order and short range order

Between the fully ordered and the random alloy lie partially ordered states.
We start from Cu$_3$Au and exchange randomly chosen Au atoms on the cube corners with randomly chosen Cu atoms on the face centres.
We measure the order with the Bragg–Williams long range order parameter [2]

$$
S = \frac{r_\mathrm{Au} - c_\mathrm{Au}}{1 - c_\mathrm{Au}},
$$

where $r_\mathrm{Au}$ is the fraction of corner sites occupied by Au.
$S$ is 1 in the ordered alloy and 0 in the random one.

```{code-cell} ipython3
ordered = make_crystal("l12", lattice_constant=3.75, repetitions=(8, 8, 8),
                       element=["Au", "Cu"])
corner = ordered.numbers == 79                 # the Au sublattice of the ordered alloy

def partially_ordered(order, seed=1):
    """Exchange Au and Cu atoms between the sublattices until S = order."""
    rng = np.random.default_rng(seed)
    alloy = ordered.copy()
    n = int(round(0.75 * (1 - order) * corner.sum()))
    gold = rng.choice(np.flatnonzero(corner), n, replace=False)
    copper = rng.choice(np.flatnonzero(~corner), n, replace=False)
    alloy.numbers[gold], alloy.numbers[copper] = 29, 79
    return alloy

def long_range_order(alloy):
    r = np.mean(alloy.numbers[corner] == 79)
    return (r - 0.25) / 0.75
```

We compute $\alpha_n$ for the first five shells of fcc, at $a/\sqrt{2}$, $a$, $\sqrt{3/2}\,a$, $\sqrt{2}\,a$ and $\sqrt{5/2}\,a$.
Each shell is selected with a fixed cutoff and `shell_thickness`, with the limits halfway between neighboring shells.

```{code-cell} ipython3
a = 3.75
shells = a * np.sqrt([0.5, 1, 1.5, 2, 2.5, 3])
limits = np.concatenate([[1.0], 0.5 * (shells[1:] + shells[:-1])])

def warren_cowley(alloy, n):
    pyscal.find_neighbors(alloy, method="cutoff", cutoff=limits[n - 1],
                          shell_thickness=limits[n] - limits[n - 1])
    return pyscal.short_range_order(alloy, "Cu", "Au")

orders = np.linspace(0, 1, 11)
scan = []
for order in orders:
    alloy = partially_ordered(order)
    scan.append((long_range_order(alloy), warren_cowley(alloy, 1), warren_cowley(alloy, 2)))
scan = np.array(scan)

states = {"ordered, S = 1": 1.0, "partially ordered, S = 0.5": 0.5, "random, S = 0": 0.0}
by_shell = {name: [warren_cowley(partially_ordered(order), n) for n in range(1, 6)]
            for name, order in states.items()}
```

```{code-cell} ipython3
:tags: [hide-input]
from _plotstyle import BLUE, PURPLE, TEAL, GREY

STATE_COLOURS = dict(zip(states, [PURPLE, TEAL, GREY]))

mp = figure(columns=2, ratio=0.42, wspace=0.3)

ax = mp[0, 0]
s = np.linspace(0, 1, 101)
ax.axhline(0, ls="--", color=DARK, lw=1)
ax.plot(s, -s**2 / 3, ls=":", color=DARK, lw=1.2)
ax.plot(s, s**2, ls=":", color=DARK, lw=1.2, label="random within sublattices")
ax.plot(scan[:, 0], scan[:, 1], "o", color=BLUE, mec=DARK, mew=0.8, ms=6,
        label=r"$\alpha_1$, first shell")
ax.plot(scan[:, 0], scan[:, 2], "s", color="white", mec=BLUE, mew=1.2, ms=6,
        label=r"$\alpha_2$, second shell")
label(ax, "(a)  Cu$_3$Au")
ax.set_xlabel("Long range order  $S$")
ax.set_ylabel(r"$\bar{\alpha}_\mathrm{CuAu}$")
ax.set_ylim(-0.5, 1.12)
ax.legend(frameon=False, loc="upper left", fontsize=9)

ax = mp[0, 1]
ax.axhline(0, ls="--", color=DARK, lw=1)
for (name, values), marker in zip(by_shell.items(), "osD"):
    ax.plot(range(1, 6), values, marker=marker, color=STATE_COLOURS[name], mec=DARK,
            mew=0.7, ms=6, lw=1.5, label=name)
label(ax, "(b)  shells")
ax.set_xlabel("Neighbor shell  $n$")
ax.set_ylabel(r"$\bar{\alpha}_n$")
ax.set_xticks(range(1, 6))
ax.set_ylim(-0.5, 1.12)
handles, labels = ax.get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

random_alloy = by_shell["random, S = 0"]
glue("alpha_random_max", round(max(abs(v) for v in random_alloy), 3), display=False)
half = by_shell["partially ordered, S = 0.5"]
glue("alpha1_half", round(half[0], 3), display=False)
glue("alpha2_half", round(half[1], 3), display=False)
glue("n_alloy", len(ordered), display=False)
```

The alloy has {glue}`n_alloy` atoms.
In panel (a), $\alpha_1$ goes from $-1/3$ in the ordered alloy to 0 in the random one, and $\alpha_2$ from 1 to 0.
The dotted lines are the values for atoms placed at random within each sublattice, $\alpha_1 = -S^2/3$ and $\alpha_2 = S^2$.
The computed values follow these lines, because the exchanges above place the atoms at random within each sublattice.
At $S = 0.5$, $\alpha_1$ is only {glue}`alpha1_half` and $\alpha_2$ is {glue}`alpha2_half`: the parameters fall with the square of the long range order.

Panel (b) shows the first five shells.
In the ordered alloy, the sign alternates.
The vectors to the second and fourth shells, of length $a$ and $\sqrt{2}\,a$, are lattice vectors of the cubic cell, which map every site onto a site of the same species.
The neighbors of a Cu atom in these shells are therefore all Cu, and $\alpha_2 = \alpha_4 = 1$.
In the first, third and fifth shells, a Cu atom has one Au neighbor for every two Cu neighbors, the most Au neighbors possible, and $\alpha_n = -1/3$.
In the random alloy, all $|\alpha_n|$ are at most {glue}`alpha_random_max`.
Above its ordering temperature, an alloy has short range order without long range order, and $\alpha_n$ decays with $n$ [1].
In the partially ordered states here, the order is long range, and $|\alpha_n|$ does not decay.

## Per-atom values

The average $\bar{\alpha}_{AB}$ over many atoms is precise, but the value of a single atom depends on only 12 neighbors.
The figure shows the distribution of $\alpha_\mathrm{CuAu}(i)$ over the Cu atoms in the three states, for the first shell.

```{code-cell} ipython3
distributions = {}
for name, order in states.items():
    alloy = partially_ordered(order)
    pyscal.find_neighbors(alloy, method="cutoff", cutoff=limits[1])
    values = pyscal.short_range_order(alloy, "Cu", "Au", average=False)
    distributions[name] = values[~np.isnan(values)]
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(ratio=0.42)
ax = mp[0, 0]
grid = 1 - np.arange(0, 13) / 3              # alpha for 0, 1, ..., 12 Au neighbors
width = 0.27 / 3
for k, (name, values) in enumerate(distributions.items()):
    fraction = [np.mean(np.isclose(values, g)) for g in grid]
    ax.bar(grid + (k - 1) * width, fraction, width=width, color=STATE_COLOURS[name],
           ec=DARK, lw=0.8, label=name, zorder=3)
ax.axvline(0, ls="--", color=DARK, lw=1)
ax.set_xlim(-1.9, 1.2)
ax.set_xlabel(r"$\alpha_\mathrm{CuAu}(i)$")
ax.set_ylabel("Fraction of Cu atoms")
secondary = ax.secondary_xaxis("top", functions=(lambda x: 3 * (1 - x), lambda n: 1 - n / 3))
secondary.set_xlabel("Au neighbors of a Cu atom")
secondary.set_xticks(range(0, 9))
ax.legend(frameon=False, loc="upper left");
```

```{code-cell} ipython3
:tags: [remove-cell]
spread = distributions["random, S = 0"]
glue("random_mean", round(spread.mean(), 3), display=False)
glue("random_std", round(spread.std(), 2), display=False)
glue("random_min", round(spread.min(), 2), display=False)
glue("random_max", round(spread.max(), 2), display=False)
```

With 12 neighbors, a Cu atom can only take the values $1 - k/3$, where $k$ is its number of Au neighbors.
In the ordered alloy every Cu atom has $k = 4$.
In the random alloy, $k$ follows a binomial distribution.
The mean of $\alpha_\mathrm{CuAu}(i)$ is {glue}`random_mean`, but the values of single atoms range from {glue}`random_min` to {glue}`random_max`, with a standard deviation of {glue}`random_std`.
A single atom with a negative value is therefore no sign of ordering.
Average the values over many atoms, over a region, or over time before interpreting them.

```{code-cell} ipython3
:tags: [remove-cell]
rattled = ordered.copy()
rattled.rattle(0.04 * a / np.sqrt(2), seed=2)
glue("alpha5_rattled", round(warren_cowley(rattled, 5), 2), display=False)
glue("alpha1_rattled", round(warren_cowley(rattled, 1), 3), display=False)
warren_cowley(rattled, 5)
glue("fifth_complete", round(100 * np.mean(pyscal.coordination_number(rattled) == 24)), display=False)
```

## Things to watch

- **The neighbor list sets the shell.** The adaptive cutoff and SANN take 14 neighbors in a perfect bcc lattice, and $\alpha_\mathrm{NiAl}$ of B2 NiAl is $-1/7$ instead of $-1$, as in the table above. Use a fixed cutoff between the shells, or `shell_thickness`, to get the parameter of a given shell.
- **Thermal displacements mix the outer shells.** At finite temperature, the distances to the outer shells overlap, and some neighbors are counted in the wrong shell. In the ordered Cu$_3$Au with random displacements of standard deviation 4 % of the nearest neighbor distance, a rough model of thermal vibrations, only {glue}`fifth_complete` % of the atoms have the 24 neighbors of the fifth shell, and $\alpha_5$ is {glue}`alpha5_rattled` instead of $-1/3$, while $\alpha_1$ is {glue}`alpha1_rattled`. When the lattice sites are known, find the neighbors on the ideal positions.
- **Like pairs have the opposite sign.** For B = A, positive values mean ordering, see [Like pairs](#like-pairs).
- **Different numbers of neighbors.** pyscal averages $\alpha_{AB}(i)$ over atoms. When the atoms have different numbers of neighbors, for example with the adaptive cutoff or Voronoi neighbors, this differs from the ratio of pair counts of Cowley, and $\bar{\alpha}_{AB} \neq \bar{\alpha}_{BA}$ in general.

## References

1. J. M. Cowley, An approximate theory of order in alloys, *Phys. Rev.* **77**, 669 (1950). [doi:10.1103/PhysRev.77.669](https://doi.org/10.1103/PhysRev.77.669)
2. B. E. Warren, *X-Ray Diffraction* (Addison-Wesley, Reading, 1969).
