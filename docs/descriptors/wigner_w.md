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

# Wigner $W_l$ parameters

The Wigner parameters $W_l$ of Steinhardt, Nelson and Ronchetti [1] are third order invariants of the same vectors $q_{lm}$ as the [Steinhardt parameters](steinhardt).
Unlike $q_l$, they can be negative, and their sign differs between structures with similar values of $q_l$.
They are used together with the Steinhardt parameters to tell fcc, hcp, bcc and icosahedral environments apart.

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

For an atom $i$ with $N(i)$ neighbors, $q_{lm}(i)$ is the average of the spherical harmonics over the bonds to the neighbors, as on the [Steinhardt parameters](steinhardt) page:

$$
q_{lm}(i) = \frac{1}{N(i)} \sum_{j=1}^{N(i)} Y_{lm}(\mathbf{r}_{ij}), \qquad m = -l, \dots, l,
$$

where $\mathbf{r}_{ij}$ is the vector from atom $i$ to its neighbor $j$.
With Voronoi neighbors, each term is multiplied by the face weight $w_{ij}$ (see [Voronoi tessellation](voronoi)), as in the [Minkowski structure metrics](minkowski).
The Wigner parameter $W_l$ is the third order invariant of this vector:

$$
W_l(i) = \sum_{\substack{m_1, m_2, m_3 \\ m_1 + m_2 + m_3 = 0}}
\begin{pmatrix} l & l & l \\ m_1 & m_2 & m_3 \end{pmatrix}
q_{lm_1}(i)\, q_{lm_2}(i)\, q_{lm_3}(i),
$$

where the term in brackets is the Wigner $3j$ symbol, and the sum runs over all $m_1$, $m_2$, $m_3$ between $-l$ and $l$ that add up to zero.
$W_l$ does not change when the crystal is rotated.
It is usually normalised by the norm of $q_{lm}$:

$$
\hat{W}_l(i) = \frac{W_l(i)}{\left( \sum_{m=-l}^{l} \left| q_{lm}(i) \right|^2 \right)^{3/2}}.
$$

$\hat{W}_l$ depends only on the shape of the arrangement of the neighbors and not on the magnitude of $q_{lm}$.
For odd $l$, the $3j$ symbol changes sign when two of its columns are exchanged, while the product of the $q_{lm}$ does not, so the terms cancel and $W_l = 0$.
The parameters with $l = 4$ and $l = 6$ are used most often.

### Averaged parameters

As for the Steinhardt parameters, Lechner and Dellago [2] proposed to average $q_{lm}$ over the atom and its neighbors before computing the invariant:

$$
\bar{q}_{lm}(i) = \frac{1}{N(i) + 1} \sum_{k=0}^{N(i)} q_{lm}(k),
$$

where $k = 0$ is the atom $i$ itself and $k = 1, \dots, N(i)$ are its neighbors.
$\bar{W}_l$ and $\hat{\bar{W}}_l$ follow from the equations above with $\bar{q}_{lm}$ in place of $q_{lm}$.

## Usage

Find the neighbors first, then compute the parameters for one or several values of $l$:

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)

w4, w6 = pyscal.wigner_w_parameter(atoms, l=[4, 6])
w4_averaged, w6_averaged = pyscal.wigner_w_parameter(atoms, l=[4, 6], averaged=True)
w6_raw = pyscal.wigner_w_parameter(atoms, l=6, normalized=False)[0]
```

`wigner_w_parameter` returns one array per value of $l$, each with one value per atom.
By default it returns the normalised $\hat{W}_l$, and with `normalized=False` the raw $W_l$.
It computes the $q_{lm}$ itself, so `steinhardt_parameter` does not have to be called first.
It stores the results on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_w4"]` | $(N,)$ | $W_4$ |
| `atoms.arrays["pyscal_what4"]` | $(N,)$ | $\hat{W}_4$ |
| `atoms.arrays["pyscal_avg_w4"]` | $(N,)$ | $\bar{W}_4$, with `averaged=True` |
| `atoms.arrays["pyscal_avg_what4"]` | $(N,)$ | $\hat{\bar{W}}_4$, with `averaged=True` |
| `atoms.arrays["pyscal_q4"]`, `["pyscal_q4_real"]`, `["pyscal_q4_imag"]` | $(N,)$, $(N, 2l + 1)$ | $q_4$ and the parts of $q_{4m}$, as stored by `steinhardt_parameter` |

Both the raw and the normalised values are stored, whatever the value of `normalized`.
With `averaged=True`, only the averaged values are stored.

## Values for perfect crystals

We compute $q_l$ and $\hat{W}_l$ for perfect fcc, bcc and hcp crystals, and for the central atom of a 13 atom icosahedral cluster.

```{code-cell} ipython3
import numpy as np
import pandas as pd
from ase.cluster import Icosahedron
from pyscal.structures import make_crystal

crystals = {
    "fcc": make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4)),
    "bcc": make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4)),
    "hcp": make_crystal("hcp", lattice_constant=2.51, repetitions=(4, 4, 4)),
}
perfect = {}
for name, crystal in crystals.items():
    pyscal.find_neighbors(crystal, method="cutoff", cutoff=0)
    q = pyscal.steinhardt_parameter(crystal, l=[4, 6])
    w = pyscal.wigner_w_parameter(crystal, l=[4, 6])
    perfect[name] = [v.mean() for v in q + w]

cluster = Icosahedron("Cu", noshells=2)
cluster.center(vacuum=8)
cluster.pbc = True
pyscal.find_neighbors(cluster, method="number", nmax=12)
centre = np.argmin(np.linalg.norm(cluster.positions - cluster.get_center_of_mass(), axis=1))
q = pyscal.steinhardt_parameter(cluster, l=[4, 6])
w = pyscal.wigner_w_parameter(cluster, l=[4, 6])
perfect["ico"] = [v[centre] for v in q + w]

pd.DataFrame(perfect, index=["q4", "q6", "w4_hat", "w6_hat"]).T.round(4)
```

The neighbors are the first shell: 12 atoms in fcc, hcp and the icosahedron, 14 (8 + 6) in bcc.
fcc and bcc have values of $\hat{W}_4$ and $\hat{W}_6$ of the same magnitude and opposite sign.
fcc and hcp have opposite signs of $\hat{W}_4$ and the same sign of $\hat{W}_6$.
The icosahedron has $q_4 = 0$, and pyscal returns $\hat{W}_4 = 0$.
Steinhardt et al. [1] used this large negative $\hat{W}_6$ to look for icosahedral order in supercooled liquids.

## Crystals at finite temperature

We compute $\hat{W}_4$ and $\hat{W}_6$ for three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.
There is no hcp snapshot, so we build a perfect hcp crystal and add random displacements with a standard deviation of 6 % of the nearest neighbor distance in each direction, a rough model of thermal vibrations.
The neighbors are all atoms closer than the first minimum of the radial distribution function $g(r)$ (see [Finding neighbors](../guide/neighbors)).

```{code-cell} ipython3
from ase.build import bulk

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "hcp": bulk("Al", "hcp", a=2.86).repeat((8, 8, 5)),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}
snapshots["hcp"].rattle(0.06 * 2.86, seed=1)

def first_minimum(atoms):
    # radial_distribution_function replaces the neighbor list, so use a copy
    g, r = pyscal.radial_distribution_function(atoms.copy(), rmax=6.0, bins=120)
    r = r + 0.5 * (r[1] - r[0])
    peak = np.argmax(g)
    return r[peak + np.argmin(g[peak:peak + 40])]

maps = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=first_minimum(snapshot))
    plain = pyscal.wigner_w_parameter(snapshot, l=[4, 6])
    averaged = pyscal.wigner_w_parameter(snapshot, l=[4, 6], averaged=True)
    maps[name] = (plain, averaged)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.55, wspace=0.12)
for k, title in enumerate([r"(a)  $\hat{W}_l$", r"(b)  $\hat{\bar{W}}_l$"]):
    ax = mp[0, k]
    ax.axhline(0, ls=":", color=DARK, lw=0.8, zorder=1)
    ax.axvline(0, ls=":", color=DARK, lw=0.8, zorder=1)
    for name, values in maps.items():
        w4, w6 = values[k]
        ax.scatter(w4, w6, s=6, color=COLOURS[name], alpha=0.5, lw=0, label=name,
                   zorder=2 if name == "liquid" else 3)
    for name in ("fcc", "hcp", "bcc"):
        w4, w6 = perfect[name][2], perfect[name][3]
        ax.plot(w4, w6, marker="*", ms=14, color=COLOURS[name], mec=DARK, mew=0.8,
                ls="none", zorder=4)
    label(ax, title)
    ax.set_xlim(-0.2, 0.2)
    ax.set_ylim(-0.16, 0.16)
    ax.set_xticks([-0.15, 0, 0.15])
    ax.set_xlabel(r"$\hat{W}_4$" if k == 0 else r"$\hat{\bar{W}}_4$")
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel(r"$\hat{W}_6$  or  $\hat{\bar{W}}_6$")
handles, labels = mp[0, 0].get_legend_handles_labels()
legend = mp.fig.legend(handles, labels, frameon=False, ncol=4, loc="upper center",
                       bbox_to_anchor=(0.5, -0.1), markerscale=3)
for handle in legend.legend_handles:
    handle.set_alpha(1);
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

def percent(name, l, sign, averaged=True):
    values = maps[name][1 if averaged else 0][0 if l == 4 else 1]
    return round(100 * np.mean(sign * values > 0))

glue("fcc_w4_negative", percent("fcc", 4, -1), display=False)
glue("hcp_w4_positive", percent("hcp", 4, +1), display=False)
glue("bcc_w4_positive", percent("bcc", 4, +1), display=False)
glue("bcc_w6_positive", percent("bcc", 6, +1), display=False)
glue("fcc_w6_negative", percent("fcc", 6, -1), display=False)
glue("hcp_w6_negative", percent("hcp", 6, -1), display=False)
glue("liquid_w4_positive", percent("liquid", 4, +1), display=False)
glue("liquid_w6_positive", percent("liquid", 6, +1), display=False)
glue("bcc_w4_positive_plain", percent("bcc", 4, +1, averaged=False), display=False)
```

The stars mark the values of the perfect crystals from the table above.
In panel (a), the clouds of $\hat{W}_l$ are broad, and the crystals overlap with each other and with the liquid.
For example, only {glue}`bcc_w4_positive_plain` % of the bcc atoms have $\hat{W}_4 > 0$.
In panel (b), the averaged parameters of the three crystals lie in three different quadrants:

- $\hat{\bar{W}}_4 < 0$ for {glue}`fcc_w4_negative` % of the fcc atoms, and $\hat{\bar{W}}_4 > 0$ for {glue}`hcp_w4_positive` % of the hcp atoms and {glue}`bcc_w4_positive` % of the bcc atoms.
- $\hat{\bar{W}}_6 > 0$ for {glue}`bcc_w6_positive` % of the bcc atoms, and $\hat{\bar{W}}_6 < 0$ for {glue}`fcc_w6_negative` % of the fcc atoms and {glue}`hcp_w6_negative` % of the hcp atoms.

The liquid is spread over all quadrants.
Of its atoms, {glue}`liquid_w4_positive` % have $\hat{\bar{W}}_4 > 0$ and {glue}`liquid_w6_positive` % have $\hat{\bar{W}}_6 > 0$.
$\hat{W}_l$ measures the shape of the environment and not the degree of order, so it does not separate a liquid from a crystal.
Separate solid from liquid atoms first, for example with $\bar{q}_6$ (see [Steinhardt parameters](steinhardt)), and then use the signs of $\hat{\bar{W}}_4$ and $\hat{\bar{W}}_6$ to assign the crystal structure.

## bcc and the choice of neighbors

The sign of $\hat{W}_4$ in bcc depends on the second neighbor shell.
For the 8 atoms of the first shell alone, which sit at the corners of a cube, $\hat{W}_4$ has the same value as in fcc.
With the full first and second shell, it has the opposite sign:

```{code-cell} ipython3
bcc = crystals["bcc"]
shells = {}
for n in (8, 14):
    pyscal.find_neighbors(bcc, method="number", nmax=n)
    shells[n] = pyscal.wigner_w_parameter(bcc, l=[4, 6])[0].mean()
{n: round(w4, 4) for n, w4 in shells.items()}
```

At finite temperature, the result therefore depends on how many atoms of the second shell the neighbor method includes, and how it weights them.
We compute $\hat{\bar{W}}_4$ and $\hat{\bar{W}}_6$ for the bcc snapshot with four neighbor methods.

```{code-cell} ipython3
bcc = snapshots["bcc"]
methods = {
    "fixed cutoff": dict(method="cutoff", cutoff=first_minimum(bcc)),
    "adaptive": dict(method="cutoff", cutoff=0),
    "SANN": dict(method="cutoff", cutoff="sann"),
    "Voronoi": dict(method="voronoi"),
}
by_method = {}
for method, options in methods.items():
    pyscal.find_neighbors(bcc, **options)
    neighbors = np.diff(bcc.info["pyscal_bond_offsets"])
    by_method[method] = (neighbors, pyscal.wigner_w_parameter(bcc, l=[4, 6], averaged=True))
```

```{code-cell} ipython3
:tags: [hide-input]
from _plotstyle import METHODS

mp = figure(columns=2, ratio=0.42, wspace=0.12)
for k, (title, bins) in enumerate([(r"(a)  $\hat{\bar{W}}_4$", np.linspace(-0.2, 0.2, 61)),
                                   (r"(b)  $\hat{\bar{W}}_6$", np.linspace(-0.01, 0.02, 61))]):
    ax = mp[0, k]
    for (method, (neighbors, values)), ls in zip(by_method.items(), ["-", "--", ":", "-."]):
        colour, marker = METHODS[method]
        ax.hist(values[k], bins=bins, density=True, histtype="step", color=colour, lw=1.6,
                ls=ls, label=method)
    ax.axvline(0, ls=":", color=DARK, lw=0.8)
    label(ax, title)
    ax.set_xlabel(title.split()[-1])
    ax.set_yticks([])
mp[0, 0].set_ylabel("Probability density")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=4, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
for method, key in [("fixed cutoff", "fixed"), ("adaptive", "adaptive"),
                    ("SANN", "sann"), ("Voronoi", "voronoi")]:
    neighbors, (w4, w6) = by_method[method]
    glue(f"bcc_{key}_neighbors", round(neighbors.mean(), 1), display=False)
    glue(f"bcc_{key}_w4_positive", round(100 * np.mean(w4 > 0)), display=False)
    glue(f"bcc_{key}_w6_positive", round(100 * np.mean(w6 > 0)), display=False)
```

| Method | Mean neighbors per atom | Atoms with $\hat{\bar{W}}_4 > 0$ (%) | Atoms with $\hat{\bar{W}}_6 > 0$ (%) |
|---|---|---|---|
| fixed cutoff | {glue}`bcc_fixed_neighbors` | {glue}`bcc_fixed_w4_positive` | {glue}`bcc_fixed_w6_positive` |
| adaptive | {glue}`bcc_adaptive_neighbors` | {glue}`bcc_adaptive_w4_positive` | {glue}`bcc_adaptive_w6_positive` |
| SANN | {glue}`bcc_sann_neighbors` | {glue}`bcc_sann_w4_positive` | {glue}`bcc_sann_w6_positive` |
| Voronoi | {glue}`bcc_voronoi_neighbors` | {glue}`bcc_voronoi_w4_positive` | {glue}`bcc_voronoi_w6_positive` |

Only the fixed cutoff at the first minimum of $g(r)$ includes the full second shell with the same weight as the first, and gives $\hat{\bar{W}}_4 > 0$ for almost all bcc atoms.
The adaptive cutoff with its default parameters and SANN miss part of the second shell.
Voronoi neighbors include the second shell, but the faces shared with the second shell are smaller than those shared with the first, so their weights are smaller.
With these three methods, many or all bcc atoms have $\hat{\bar{W}}_4 < 0$, as in fcc.
The sign of $\hat{\bar{W}}_6$ does not have this problem, and it is positive for the bcc atoms with every method (last column).

## Things to watch

- **The neighbor method changes the values.** $\hat{W}_4$ of bcc changes sign depending on how the second shell is counted, as shown above. Compare values only when they were computed with the same neighbor method, and check $\hat{W}_4$ in a perfect crystal with the same method before using it as a criterion.
- **Use the averaged parameters at finite temperature.** The plain $\hat{W}_l$ of the crystals overlap at finite temperature. The averaged $\hat{\bar{W}}_l$ separate them.
- **$\hat{W}_l$ does not measure order.** A liquid atom can have any value of $\hat{W}_l$. Combine $\hat{W}_l$ with $q_l$ or $\bar{q}_l$.
- **Odd $l$.** $W_l$ is zero for odd $l$, and `wigner_w_parameter` returns zeros.
- **Undefined normalisation.** When $\sum_m |q_{lm}|^2$ is zero, as for $l = 4$ in an icosahedron, $\hat{W}_l$ is undefined and `wigner_w_parameter` returns 0.

## References

1. P. J. Steinhardt, D. R. Nelson and M. Ronchetti, Bond-orientational order in liquids and glasses, *Phys. Rev. B* **28**, 784 (1983). [doi:10.1103/PhysRevB.28.784](https://doi.org/10.1103/PhysRevB.28.784)
2. W. Lechner and C. Dellago, Accurate determination of crystal structures based on averaged local bond order parameters, *J. Chem. Phys.* **129**, 114707 (2008). [doi:10.1063/1.2977970](https://doi.org/10.1063/1.2977970)
