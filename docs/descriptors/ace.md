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

# Atomic Cluster Expansion descriptors

The Atomic Cluster Expansion (ACE) of Drautz [1] describes the neighborhood of an atom by a vector of many numbers.
The vector does not change when the structure is rotated or translated, or when the neighbors are relabelled.
ACE is the basis of a family of machine learning interatomic potentials, and its descriptors are used as input to models that classify atoms or predict their properties.
Unlike [Steinhardt parameters](steinhardt) or [common neighbor analysis](cna), which reduce the neighborhood to a few numbers chosen by hand, ACE keeps radial and angular information up to a chosen resolution, and leaves it to a model to decide what matters.

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

For an atom $i$, the atomic base $A$ projects the positions of the neighbors $j$ on radial functions $R_n$ and spherical harmonics $Y_{lm}$:

$$
A_{nlm}(i) = \sum_{j} R_n(r_{ij})\, Y_{lm}(\hat{\mathbf{r}}_{ij}), \qquad n = 0, \dots, n_\mathrm{max} - 1, \quad l = 0, \dots, l_\mathrm{max}, \quad m = -l, \dots, l,
$$

where $r_{ij}$ is the distance to neighbor $j$, $\hat{\mathbf{r}}_{ij}$ the unit vector pointing to it, and the sum runs over the neighbors with $r_{ij} < r_c$.
The radial functions are Chebyshev polynomials $T_n$ multiplied by a smooth cutoff function $f_c$:

$$
R_n(r) = T_n(x)\, f_c(r), \qquad x = 2\, \frac{r - r_\mathrm{min}}{r_c - r_\mathrm{min}} - 1, \qquad f_c(r) = \frac{1}{2} \left[ \cos\left( \frac{\pi r}{r_c} \right) + 1 \right],
$$

with $r_\mathrm{min} = 0.5$ Å.
$f_c$ goes smoothly to zero at the cutoff $r_c$, so a neighbor that crosses the cutoff does not change the descriptors abruptly.

The $A_{nlm}$ change when the structure is rotated.
Products of $A_{nlm}$, summed over $m$ with the right coefficients, do not [1, 2].
The number of factors in the product is the body order $\nu$, and pyscal computes $\nu = 1$, 2 and 3:

$$
B^{(1)}_{n}(i) = \mathrm{Re}\, A_{n00}(i),
$$

$$
B^{(2)}_{n_1 n_2 l}(i) = \sum_{m=-l}^{l} \mathrm{Re} \left[ A^*_{n_1 l m}(i)\, A_{n_2 l m}(i) \right], \qquad n_1 \le n_2,
$$

$$
B^{(3)}_{n_1 n_2 n_3 l_1 l_2 l_3}(i) = \sum_{m_1, m_2, m_3} \begin{pmatrix} l_1 & l_2 & l_3 \\ m_1 & m_2 & m_3 \end{pmatrix} \mathrm{Re} \left[ A_{n_1 l_1 m_1}(i)\, A_{n_2 l_2 m_2}(i)\, A_{n_3 l_3 m_3}(i) \right],
$$

where the matrix in the last line is the Wigner 3j symbol, which is zero unless $m_1 + m_2 + m_3 = 0$.
$B^{(1)}$ is a radial density of the neighbors, $B^{(2)}$ is the power spectrum, as in SOAP [3], and $B^{(3)}$ is the bispectrum.
pyscal computes $B^{(3)}$ for $n_1 \le n_2 \le n_3 \le 2$ and $l_1, l_2, l_3 \le 2$, with $l_1 + l_2 + l_3$ even and $|l_1 - l_2| \le l_3 \le l_1 + l_2$, to limit the cost.
The descriptor vector of an atom is the concatenation of the $B^{(\nu)}$ up to $\nu_\mathrm{max}$.
By default, it is divided by its length, so that each atom has a vector of length 1.

The power spectrum contains the [Steinhardt parameters](steinhardt).
If all $N(i)$ neighbors are at the same distance $r$, as in the first shell of a perfect crystal, $R_0 = f_c$ and

$$
B^{(2)}_{00l}(i) = \frac{2l + 1}{4\pi}\, N(i)^2 f_c(r)^2\, q_l(i)^2 .
$$

pyscal implements ACE for a single chemical species, with the radial functions above.
The values cannot be compared with those of other ACE codes, which use other radial functions and normalisations.

## Usage

Find the neighbors first, with a cutoff at least as large as the cutoff of the descriptors:

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=5.5)

result = pyscal.ace(atoms, nmax=4, lmax=4, nu_max=2, cutoff=5.5)
{key: value.shape for key, value in result.items()}
```

```{code-cell} ipython3
:tags: [remove-cell]
# check the relation between the power spectrum and the Steinhardt parameters
from ase.build import bulk
perfect = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4)
pyscal.find_neighbors(perfect, method="cutoff", cutoff=3.0)
q6 = pyscal.steinhardt_parameter(perfect, l=6)[0]
b2 = pyscal.ace(perfect, nmax=2, lmax=6, nu_max=2, normalize=False)["nu2"][:, 6]
r = 3.61 / 2 ** 0.5
fc = 0.5 * (np.cos(np.pi * r / 3.0) + 1)
assert np.allclose(b2, 13 / (4 * np.pi) * 12 ** 2 * fc ** 2 * q6 ** 2)
```

`ace` returns a dictionary with one array per body order and their concatenation.
It stores the concatenated vector and the parameters on `atoms`:

| Key | Shape | Content |
|---|---|---|
| `result["nu1"]` | $(N, n_\mathrm{max})$ | $B^{(1)}$ |
| `result["nu2"]` | $(N, \tfrac{1}{2} n_\mathrm{max}(n_\mathrm{max} + 1)(l_\mathrm{max} + 1))$ | $B^{(2)}$, with `nu_max` ≥ 2 |
| `result["nu3"]` | $(N, n_3)$ | $B^{(3)}$, with `nu_max` ≥ 3 |
| `result["full"]` | $(N, n_\mathrm{features})$ | all blocks, concatenated |
| `atoms.arrays["pyscal_ace"]` | $(N, n_\mathrm{features})$ | the same as `result["full"]` |
| `atoms.info["pyscal_ace_params"]` | dictionary | `nmax`, `lmax`, `nu_max` and `cutoff` |

With `normalize=True` (the default), all blocks are divided by the length of the full vector of the atom.
`normalize=False` returns the unscaled values.

`find_neighbors` must be called first.
Only the neighbors in the neighbor list contribute, so the neighbor list must include all atoms within the cutoff.
Without `cutoff`, `ace` uses the largest cutoff stored by `find_neighbors`, or 5 Å if there is none.

## Telling structures apart

Do the ACE descriptors separate the atoms of an fcc crystal, a bcc crystal and a liquid?
We compute them for the three MD snapshots in the `examples` folder.
The snapshots have different densities, and the descriptors depend on absolute distances.
To compare only the arrangement of the atoms, we scale each snapshot to the volume per atom of the fcc snapshot.

```{code-cell} ipython3
snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
    "liquid": read("conf.lqd.Al.dump", format="lammps-dump-text"),
}
volume = snapshots["fcc"].get_volume() / len(snapshots["fcc"])
for snapshot in snapshots.values():
    scale = (volume / (snapshot.get_volume() / len(snapshot))) ** (1 / 3)
    snapshot.set_cell(snapshot.cell * scale, scale_atoms=True)
    pyscal.find_neighbors(snapshot, method="cutoff", cutoff=5.5)

X = np.vstack([pyscal.ace(s, nmax=4, lmax=4, nu_max=2, cutoff=5.5)["full"]
               for s in snapshots.values()])
y = np.concatenate([[name] * len(s) for name, s in snapshots.items()])
X.shape
```

Each atom has a vector of {glue}`n_features` numbers.
We look at them in two ways with [scikit-learn](https://scikit-learn.org).
Principal component analysis (PCA) finds the two directions in which the vectors vary most, without using the labels.
Linear discriminant analysis (LDA) finds the two directions that separate the labelled groups best.
We fit LDA on a random half of the atoms and test it on the other half.

```{code-cell} ipython3
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

pca = PCA(n_components=2).fit(X)

train = np.random.default_rng(0).random(len(y)) < 0.5
lda = LinearDiscriminantAnalysis(n_components=2).fit(X[train], y[train])
lda.score(X[~train], y[~train])
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
glue("n_features", X.shape[1], display=False)
glue("pca_variance", round(100 * pca.explained_variance_ratio_.sum()), display=False)
def percent(x):
    value = round(100 * x, 1)
    return int(value) if value == int(value) else value
glue("lda_accuracy", percent(lda.score(X[~train], y[~train])), display=False)
glue("n_test", int((~train).sum()), display=False)
predicted = lda.predict(X[~train])
for name in snapshots:
    sel = y[~train] == name
    glue(f"lda_{name}", percent(np.mean(predicted[sel] == name)), display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
projections = [("(a)  PCA", pca.transform(X), ~train | train, "PC"),
               ("(b)  LDA, test atoms", lda.transform(X), ~train, "LD")]
mp = figure(columns=2, ratio=0.5, wspace=0.3)
for k, (title, Z, shown, axis) in enumerate(projections):
    ax = mp[0, k]
    for name in snapshots:
        sel = shown & (y == name)
        ax.scatter(*Z[sel].T, s=6, color=COLOURS[name], alpha=0.5, lw=0, label=name,
                   zorder=2 if name == "liquid" else 3)
    label(ax, title)
    ax.set_xlabel(f"{axis} 1")
    ax.set_ylabel(f"{axis} 2")
legend = mp.fig.legend(*mp[0, 0].get_legend_handles_labels(), frameon=False, ncol=3,
                       loc="upper center", bbox_to_anchor=(0.5, -0.1), markerscale=3)
for handle in legend.legend_handles:
    handle.set_alpha(1)
```

Panel (a) shows the first two principal components, which hold {glue}`pca_variance` % of the variance.
The liquid atoms spread widely, and the crystal atoms form a dense cloud in which fcc and bcc overlap.
The directions of largest variance mostly separate the disordered liquid from the crystals, not the two crystals from each other.

Panel (b) shows the LDA projection of the {glue}`n_test` test atoms, which were not used in the fit.
The three groups are separated.
LDA assigns {glue}`lda_accuracy` % of the test atoms to the right structure: {glue}`lda_fcc` % of the fcc atoms, {glue}`lda_bcc` % of the bcc atoms and {glue}`lda_liquid` % of the liquid atoms.
The descriptors of single atoms contain enough information to tell the three structures apart, even with the thermal disorder of the snapshots.
This information is spread over many components, and PCA alone does not find it.

## Choosing nmax, lmax and nu_max

The number of features grows with the three parameters.
$B^{(1)}$ has $n_\mathrm{max}$ features, and $B^{(2)}$ has $\tfrac{1}{2} n_\mathrm{max} (n_\mathrm{max} + 1)(l_\mathrm{max} + 1)$.
$B^{(3)}$ has a fixed number of features once $n_\mathrm{max} \ge 3$ and $l_\mathrm{max} \ge 2$, because of the limits on $n$ and $l$ given above.
We repeat the LDA test of the previous section for several values of $l_\mathrm{max}$ and $\nu_\mathrm{max}$, with $n_\mathrm{max} = 4$.

```{code-cell} ipython3
lmaxs = [0, 1, 2, 3, 4, 6, 8]
scan = {}
for nu_max in (1, 2, 3):
    for lmax in lmaxs:
        X = np.vstack([pyscal.ace(s, nmax=4, lmax=lmax, nu_max=nu_max, cutoff=5.5)["full"]
                       for s in snapshots.values()])
        lda = LinearDiscriminantAnalysis().fit(X[train], y[train])
        scan[nu_max, lmax] = (X.shape[1], lda.score(X[~train], y[~train]))
```

```{code-cell} ipython3
:tags: [remove-cell]
from _plotstyle import RED, TEAL, PURPLE
NU_STYLE = {1: (RED, "o"), 2: (TEAL, "s"), 3: (PURPLE, "D")}
glue("acc_nu1", percent(scan[1, 4][1]), display=False)
glue("acc_nu2_l2", percent(scan[2, 2][1]), display=False)
glue("acc_nu2_l4", percent(scan[2, 4][1]), display=False)
glue("acc_nu3_l4", percent(scan[3, 4][1]), display=False)
glue("n_nu3", scan[3, 4][0] - scan[2, 4][0], display=False)
glue("n_nu2_l4", scan[2, 4][0], display=False)
glue("n_nu3_l4", scan[3, 4][0], display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.3)
for k, (title, index, ylabel) in enumerate([("(a)", 0, "Number of features"),
                                            ("(b)", 1, "Test atoms assigned correctly  (%)")]):
    ax = mp[0, k]
    for nu_max, (colour, marker) in NU_STYLE.items():
        values = [scan[nu_max, l][index] for l in lmaxs]
        if index:
            values = 100 * np.array(values)
        ax.plot(lmaxs, values, marker=marker, color=colour, mec=DARK, mew=0.7, ms=5,
                lw=1.6, label=rf"$\nu_\mathrm{{max}}$ = {nu_max}")
    label(ax, title)
    ax.set_xlabel(r"$l_\mathrm{max}$")
    ax.set_ylabel(ylabel)
    ax.set_xticks(lmaxs)
mp[0, 0].set_ylim(0, None)
mp[0, 1].set_ylim(80, 101)
mp.fig.legend(*mp[0, 0].get_legend_handles_labels(), frameon=False, ncol=3,
              loc="upper center", bbox_to_anchor=(0.5, -0.1));
```

With $\nu_\mathrm{max} = 1$, the descriptors contain only radial information and do not depend on $l_\mathrm{max}$.
They assign {glue}`acc_nu1` % of the test atoms correctly.
The angular information of $B^{(2)}$ raises this to {glue}`acc_nu2_l2` % with $l_\mathrm{max} = 2$ and {glue}`acc_nu2_l4` % with $l_\mathrm{max} = 4$.
As for the Steinhardt parameters, $l = 4$ is needed to tell the cubic structures apart.
Larger $l_\mathrm{max}$ adds features but hardly improves the result.
$B^{(3)}$ adds {glue}`n_nu3` features, which is more than the {glue}`n_nu2_l4` features of $\nu_\mathrm{max} = 2$ and $l_\mathrm{max} = 4$, and gives {glue}`acc_nu3_l4` % instead of {glue}`acc_nu2_l4` %.
Here, $\nu_\mathrm{max} = 2$ with $l_\mathrm{max} = 4$ is enough.

The computing time grows with the number of neighbors and with the number of spherical harmonics, $(l_\mathrm{max} + 1)^2$.

## Things to watch

- **The neighbor list must reach the cutoff.** Only neighbors found by `find_neighbors` contribute. If the neighbor list is shorter than the cutoff of `ace`, the neighbors in between are missing and $f_c$ does not go to zero at the edge of the neighbor list. Use `method="cutoff"` with the same cutoff in both calls. Without `cutoff`, `ace` takes the largest cutoff stored by `find_neighbors`, which for the other neighbor methods depends on the structure. Pass `cutoff` explicitly.
- **Chemical species are ignored.** All neighbors contribute in the same way, whatever their species.
- **Density.** The descriptors depend on absolute distances. Structures with different densities, or the same structure at different temperatures, differ also through the density. Scale the structures to the same volume per atom when only the arrangement of the atoms matters.
- **The limits of $B^{(3)}$.** The body order 3 block uses only $n \le 2$ and $l \le 2$. Values of `nu_max` larger than 3 give the same result as 3.
- **Normalisation.** With `normalize=True`, the length of the vector, which grows with the number of neighbors, is lost. Use `normalize=False` if it matters, for example to predict a property that depends on the number of neighbors.

## References

1. R. Drautz, Atomic cluster expansion for accurate and transferable interatomic potentials, *Phys. Rev. B* **99**, 014104 (2019). [doi:10.1103/PhysRevB.99.014104](https://doi.org/10.1103/PhysRevB.99.014104)
2. G. Dusson, M. Bachmayr, G. Csányi, R. Drautz, S. Etter, C. van der Oord and C. Ortner, Atomic cluster expansion: Completeness, efficiency and stability, *J. Comput. Phys.* **454**, 110946 (2022). [doi:10.1016/j.jcp.2022.110946](https://doi.org/10.1016/j.jcp.2022.110946)
3. A. P. Bartók, R. Kondor and G. Csányi, On representing chemical environments, *Phys. Rev. B* **87**, 184115 (2013). [doi:10.1103/PhysRevB.87.184115](https://doi.org/10.1103/PhysRevB.87.184115)
