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

# Entropy parameter

The entropy parameter of Piaggi and Parrinello [1] is a fingerprint of the local order around an atom.
It is computed from the radial distribution function centred on the atom, in the form of the pair contribution to the entropy.
It is lower in a crystal than in a liquid, and it is used to separate solid from liquid atoms and to find defects [1].
Unlike the [Steinhardt parameters](steinhardt), the [disorder parameter](disorder.md) and the [solid–liquid classification](solid_liquid), it uses only the distances to the neighbors and no bond angles.

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

The radial distribution function centred on atom $i$, smoothed with Gaussians of width $\sigma$, is

$$
g_m^i(r) = \frac{1}{4 \pi \rho r^2} \sum_{j} \frac{1}{\sqrt{2 \pi \sigma^2}} \exp\left( -\frac{(r - r_{ij})^2}{2 \sigma^2} \right),
$$

where the sum runs over the neighbors $j$ of atom $i$, $r_{ij}$ is the distance between $i$ and $j$, and $\rho = N / V$ is the number density of the system, with $N$ atoms in the volume $V$.
Piaggi and Parrinello [1] define the entropy parameter as

$$
s_s^i = -2 \pi \rho k_\mathrm{B} \int_0^{r_m} \left[ g_m^i(r) \ln g_m^i(r) - g_m^i(r) + 1 \right] r^2 \, \mathrm{d}r,
$$

where $r_m$ is the upper limit of the integral and $k_\mathrm{B}$ is the Boltzmann constant.
pyscal leaves out the factor $2 \pi k_\mathrm{B}$ and computes the dimensionless value

$$
s(i) = -\rho \int_{r_0}^{r_m} \left[ g_m^i(r) \ln g_m^i(r) - g_m^i(r) + 1 \right] r^2 \, \mathrm{d}r = \frac{s_s^i}{2 \pi k_\mathrm{B}}.
$$

The integral is evaluated with the trapezoidal rule from $r_0$ to $r_m$ in steps of $h$.
The term in square brackets is zero for $g_m^i = 1$ and positive otherwise, so $s(i) \le 0$.
The more structured the surroundings of atom $i$, the more negative $s(i)$.

The averaged entropy parameter averages $s$ over the atom and its neighbors:

$$
\bar{s}(i) = \frac{1}{N(i) + 1} \sum_{k=0}^{N(i)} s(k),
$$

where $N(i)$ is the number of neighbors of $i$, $k = 0$ is the atom $i$ itself and $k = 1, \dots, N(i)$ are its neighbors.

## Usage

Find the neighbors with a fixed cutoff of at least $r_m + 3\sigma$ (see [Things to watch](#things-to-watch)), then compute the entropy parameter:

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
rm, sigma = 5.7, 0.2
pyscal.find_neighbors(atoms, method="cutoff", cutoff=rm + 3 * sigma)

s = pyscal.entropy(atoms, rm=rm, sigma=sigma)
s_averaged = pyscal.entropy(atoms, rm=rm, sigma=sigma, averaged=True)
```

The arguments are those of the definition, with lengths in the units of the positions, usually Å:

| Argument | Default | Meaning |
|---|---|---|
| `rm` | required | upper limit $r_m$ of the integral |
| `sigma` | 0.2 | width $\sigma$ of the Gaussians |
| `rstart` | 0.001 | lower limit $r_0$ of the integral |
| `h` | 0.001 | step $h$ of the integration |
| `averaged` | `False` | return $\bar{s}$ instead of $s$ (`average` does the same) |
| `local` | `False` | use a local density for each atom instead of $\rho$, see below |

`entropy` returns $s$, or $\bar{s}$ with `averaged=True`, with one value per atom.
It also stores the results on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `atoms.arrays["pyscal_entropy"]` | $(N,)$ | $s$ |
| `atoms.arrays["pyscal_average_entropy"]` | $(N,)$ | $\bar{s}$, with `averaged=True` |

With `local=True`, the density of each atom is $\rho_i = N(i) / (\tfrac{4}{3} \pi r_\mathrm{c}(i)^3)$, where $r_\mathrm{c}(i)$ is its neighbor cutoff, instead of the global density $\rho = N / V$.
`entropy` needs a fully periodic cell.

## Solid and liquid

We compute $s$ and $\bar{s}$ for three MD snapshots from the `examples` folder: an fcc crystal and a bcc crystal at finite temperature, and a liquid.
The snapshots are of different elements, with different densities.
If all positions, $r_m$ and $\sigma$ are scaled by the same factor, $s$ does not change.
To compare the snapshots, we therefore give $r_m$ and $\sigma$ in units of $\ell = (V / N)^{1/3}$, the cube root of the volume per atom: $r_m = 2.2\,\ell$ and $\sigma = 0.08\,\ell$.

```{code-cell} ipython3
import numpy as np

snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}

def entropies(snapshot, rm_scaled, sigma_scaled):
    """s and s_averaged with rm and sigma given in units of (V/N)^(1/3)"""
    ell = (snapshot.get_volume() / len(snapshot)) ** (1 / 3)
    rm, sigma = rm_scaled * ell, sigma_scaled * ell
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=rm + 3 * sigma)
    s_averaged = pyscal.entropy(snapshot, rm=rm, sigma=sigma, averaged=True)
    return snapshot.arrays["pyscal_entropy"], s_averaged

values = {name: entropies(snapshot, 2.2, 0.08) for name, snapshot in snapshots.items()}
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

ell_liquid = (snapshots["liquid"].get_volume() / len(snapshots["liquid"])) ** (1 / 3)
glue("rm_liquid", round(2.2 * ell_liquid, 1), display=False)
glue("sigma_liquid", round(0.08 * ell_liquid, 2), display=False)

crystal_max = [max(values[name][k].max() for name in ("fcc", "bcc")) for k in (0, 1)]
liquid_min = [values["liquid"][k].min() for k in (0, 1)]
glue("crystal_s_max", round(crystal_max[0], 2), display=False)
glue("liquid_s_min", round(liquid_min[0], 2), display=False)
glue("liquid_s_overlap",
     round(100 * np.mean(values["liquid"][0] < crystal_max[0])), display=False)
glue("crystal_sbar_max", round(crystal_max[1], 2), display=False)
glue("liquid_sbar_min", round(liquid_min[1], 2), display=False)
assert crystal_max[1] < liquid_min[1]

# s does not change when all lengths are scaled by the same factor
scaled = snapshots["liquid"].copy()
scaled.set_cell(1.3 * scaled.cell, scale_atoms=True)
assert np.allclose(entropies(scaled, 2.2, 0.08)[0], values["liquid"][0], atol=1e-3)
```

For the liquid, these are $r_m =$ {glue}`rm_liquid` Å and $\sigma =$ {glue}`sigma_liquid` Å.

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.25)
bins = np.linspace(-1.4, 0, 57)
for k, title in enumerate(["(a)  $s$", r"(b)  $\bar{s}$"]):
    ax = mp[0, k]
    for name, data in values.items():
        ax.hist(data[k], bins=bins, density=True, histtype="stepfilled",
                color=COLOURS[name], alpha=0.4)
        ax.hist(data[k], bins=bins, density=True, histtype="step",
                color=COLOURS[name], lw=1.4, label=name)
    label(ax, title)
    ax.set_xlabel(title.split()[-1])
    ax.set_xlim(-1.4, 0)
mp[0, 0].set_ylabel("Probability density")
mp[0, 0].legend(frameon=False, loc="upper left");
```

In panel (a), the distributions of $s$ overlap.
The largest value of $s$ in the two crystals is {glue}`crystal_s_max`, and {glue}`liquid_s_overlap` % of the liquid atoms have a lower value.
In panel (b), the averaged parameter $\bar{s}$ is at most {glue}`crystal_sbar_max` in the crystals and at least {glue}`liquid_sbar_min` in the liquid, so the distributions are separated.
As for the [Steinhardt parameters](steinhardt), use the averaged values to classify atoms.

## The choice of $r_m$ and $\sigma$

The values of $s$ depend strongly on $r_m$ and $\sigma$.
We vary each of them, with the other one fixed at the value used above, and plot the mean of $\bar{s}$ in each snapshot.

```{code-cell} ipython3
rm_values = np.arange(1.4, 3.21, 0.1)
sigma_values = np.arange(0.03, 0.201, 0.01)

rm_scan, sigma_scan = {}, {}
for name, snapshot in snapshots.items():
    rm_scan[name] = [entropies(snapshot, rm, 0.08)[1] for rm in rm_values]
    sigma_scan[name] = [entropies(snapshot, 2.2, sigma)[1] for sigma in sigma_values]
```

```{code-cell} ipython3
:tags: [remove-cell]
def separated(scan):
    """True if every value of s_averaged in the crystals is below every value in the liquid"""
    return all(max(scan["fcc"][k].max(), scan["bcc"][k].max()) < scan["liquid"][k].min()
               for k in range(len(scan["liquid"])))

assert separated(rm_scan) and separated(sigma_scan)
glue("sigma_min", round(sigma_values[0], 2), display=False)
glue("sigma_max", round(sigma_values[-1], 2), display=False)
glue("liquid_sigma_min", round(sigma_scan["liquid"][0].mean(), 2), display=False)
glue("fcc_sigma_008", round(values["fcc"][1].mean(), 2), display=False)
assert sigma_scan["liquid"][0].mean() < values["fcc"][1].mean()
for name in snapshots:
    assert np.all(np.diff([v.mean() for v in rm_scan[name]]) < 0)
    assert np.all(np.diff([v.mean() for v in sigma_scan[name]]) > 0)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.12)
scans = [(rm_values, rm_scan, r"$r_m\,/\,\ell$", 2.2), (sigma_values, sigma_scan, r"$\sigma\,/\,\ell$", 0.08)]
for k, (x, scan, xlabel, used) in enumerate(scans):
    ax = mp[0, k]
    ax.axvline(used, ls=":", color=DARK, lw=1)
    for name, marker in zip(snapshots, "oDs"):
        ax.plot(x, [v.mean() for v in scan[name]], marker=marker, color=COLOURS[name],
                mec=DARK, mew=0.7, ms=4.5, lw=1.6, label=name, zorder=3)
        ax.fill_between(x, [v.min() for v in scan[name]], [v.max() for v in scan[name]],
                        color=COLOURS[name], alpha=0.25, lw=0, zorder=2)
    label(ax, f"({'ab'[k]})  varying {xlabel}")
    ax.set_xlabel(xlabel)
    ax.set_ylim(-1.6, 0)
    if k:
        ax.tick_params(labelleft=False)
mp[0, 0].set_ylabel(r"$\bar{s}$")
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

The lines show the mean of $\bar{s}$, the shaded bands the range from its smallest to its largest value, and the dotted lines the values used above.
$\bar{s}$ becomes more negative with larger $r_m$, as more neighbor shells contribute to the integral.
It becomes much more negative with smaller $\sigma$, because narrower Gaussians make $g_m^i(r)$ more structured, in the liquid as well.
For example, the mean of $\bar{s}$ in the liquid at $\sigma = $ {glue}`sigma_min` $\ell$ is {glue}`liquid_sigma_min`, lower than the mean in the fcc crystal at $\sigma = 0.08\,\ell$, {glue}`fcc_sigma_008`.
Values of $s$ can therefore be compared only when they were computed with the same $r_m$ and $\sigma$.
Over the whole range shown, the crystals and the liquid stay separated.

## Things to watch

```{code-cell} ipython3
:tags: [remove-cell]
liquid = snapshots["liquid"]
ell = (liquid.get_volume() / len(liquid)) ** (1 / 3)
rm, sigma = 2.2 * ell, 0.08 * ell
means = []
for cutoff in (rm, rm + 3 * sigma, rm + 3.0):
    pyscal.find_neighbors(liquid, method="cutoff", cutoff=cutoff)
    means.append(pyscal.entropy(liquid, rm=rm, sigma=sigma).mean())
glue("cut_rm", round(means[0], 3), display=False)
glue("cut_rm3s", round(means[1], 3), display=False)
glue("cut_large", round(means[2], 3), display=False)
```

- **The neighbor cutoff.** Only neighbors contribute to $g_m^i(r)$, and atoms just beyond $r_m$ contribute to it through the tails of their Gaussians. Use a neighbor cutoff of at least $r_m + 3\sigma$, so that these atoms are neighbors. For the liquid above, the mean of $s$ is {glue}`cut_rm` with a cutoff of $r_m$, {glue}`cut_rm3s` with $r_m + 3\sigma$ and {glue}`cut_large` with $r_m + 3$ Å. `entropy` warns only if the cutoff is smaller than $r_m$.
- **Adaptive and other neighbor methods.** The adaptive cutoff, SANN and Voronoi neighbors include about one neighbor shell, so $g_m^i(r)$ is zero beyond it. Use a fixed cutoff.
- **Same $r_m$ and $\sigma$ for all structures.** As shown above, the values depend strongly on $r_m$ and $\sigma$.
- **Units.** The values of pyscal are $s_s^i / (2 \pi k_\mathrm{B})$. Multiply by $2 \pi$ to compare with values in units of $k_\mathrm{B}$, as in [1].

## References

1. P. M. Piaggi and M. Parrinello, Entropy based fingerprint for local crystalline order, *J. Chem. Phys.* **147**, 114112 (2017). [doi:10.1063/1.4998408](https://doi.org/10.1063/1.4998408)
