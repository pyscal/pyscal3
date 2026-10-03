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

# Atomic deformation descriptors

The deformation descriptors compare each atom and its neighbors in a deformed structure with the same atoms in a reference structure.
pyscal computes four of them: the atomic strain tensor, the von Mises shear strain, the non-affine displacement $D^2_\mathrm{min}$ and the slip vector.
They are used to find shear bands and plastic events in glasses, to follow slip and stacking faults in crystals, and to measure local strain near defects.
Unlike [common neighbor analysis](cna) or the [centrosymmetry parameter](centrosymmetry.md), which look at one structure, they need a reference structure with the same atoms.
[Wigner–Seitz defect analysis](wigner_seitz) also uses a reference, but counts atoms per lattice site instead of comparing neighbors.

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

Let $\mathbf{R}_{ij}$ be the vector from atom $i$ to its neighbor $j$ in the reference structure, and $\mathbf{r}_{ij}$ the same vector in the current structure.
The descriptors of atom $i$ use the $n_i$ atoms $j$ that are neighbors of $i$ in both structures.

### Atomic strain

The local deformation gradient $\mathbf{F}(i)$ is the $3 \times 3$ matrix that maps the reference vectors onto the current vectors as well as possible [1, 2]:

$$
\mathbf{F}(i) = \arg\min_{\mathbf{F}} \sum_{j=1}^{n_i} \left| \mathbf{r}_{ij} - \mathbf{F}\,\mathbf{R}_{ij} \right|^2
= \left( \sum_{j} \mathbf{r}_{ij} \mathbf{R}_{ij}^\mathsf{T} \right) \left( \sum_{j} \mathbf{R}_{ij} \mathbf{R}_{ij}^\mathsf{T} \right)^{-1}.
$$

The atomic strain is the Green–Lagrange strain tensor of this fit:

$$
\mathbf{E}(i) = \frac{1}{2} \left( \mathbf{F}(i)^\mathsf{T} \mathbf{F}(i) - \mathbf{I} \right),
$$

where $\mathbf{I}$ is the identity matrix.
For a homogeneous deformation $\mathbf{F}$ of the whole structure, $\mathbf{F}(i) = \mathbf{F}$ and $\mathbf{E}(i)$ is the applied strain.

### Von Mises shear strain

`von_mises_strain` reduces $\mathbf{E}$ with components $E_{\alpha\beta}$ to one number per atom:

$$
\eta(i) = \left( \frac{1}{2} \left[ (E_{xx} - E_{yy})^2 + (E_{yy} - E_{zz})^2 + (E_{zz} - E_{xx})^2 \right] + E_{xy}^2 + E_{yz}^2 + E_{xz}^2 \right)^{1/2}.
$$

$\eta$ is zero for a pure change of volume.

### Non-affine displacement

Falk and Langer [1] measure how far the neighbors of an atom have moved away from the best affine fit:

$$
D^2_\mathrm{min}(i) = \frac{1}{n_i} \sum_{j=1}^{n_i} \left| \mathbf{r}_{ij} - \mathbf{F}(i)\,\mathbf{R}_{ij} \right|^2 .
$$

$D^2_\mathrm{min}$ has units of Å².
It is zero for any homogeneous deformation and large where atoms have rearranged, for example in a shear transformation or across a slip plane.
pyscal divides by $n_i$, so $D^2_\mathrm{min}$ is a mean over the neighbors.
Falk and Langer use the sum.

### Slip vector

The slip vector is the mean change of the neighbor vectors, without a fit:

$$
\mathbf{s}(i) = \frac{1}{n_i} \sum_{j=1}^{n_i} \left( \mathbf{r}_{ij} - \mathbf{R}_{ij} \right).
$$

For a homogeneous deformation of a crystal in which every neighbor at $\mathbf{R}$ has a partner at $-\mathbf{R}$, the terms cancel and $\mathbf{s} = 0$.
Next to a plane across which the crystal has slipped by a vector $\mathbf{b}$, the neighbors on the other side contribute $\pm\mathbf{b}$, and $\mathbf{s}$ points along the slip direction.
Zimmerman et al. [3] sum only over the neighbors whose vectors changed by more than a threshold and divide by their number.
pyscal averages over all $n_i$ pairs and uses no threshold, so $|\mathbf{s}|$ is a fraction of $|\mathbf{b}|$.

The strain, $\eta$ and $D^2_\mathrm{min}$ are NaN for atoms with fewer than three pairs.
The slip vector is NaN for atoms with no pairs.

## Usage

We stretch a perfect fcc crystal by 2 % along $x$ and shear it by 4 % in the $xy$ plane.
Both structures need their own neighbor list:

```{code-cell} ipython3
import numpy as np
import pyscal
from ase.build import bulk

reference = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(6)
deformed = reference.copy()
F = np.array([[1.02, 0.04, 0.0],
              [0.0,  1.0,  0.0],
              [0.0,  0.0,  1.0]])
deformed.set_cell(reference.cell @ F.T, scale_atoms=True)

for structure in (reference, deformed):
    pyscal.find_neighbors(structure, method="cutoff", cutoff=3.1)

strain = pyscal.atomic_strain(deformed, reference)
eta = pyscal.von_mises_strain(deformed, reference)
d2 = pyscal.d2min(deformed, reference)
slip = pyscal.slip_vector(deformed, reference)

strain[0].round(4)
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
for key, (r, c) in {"xx": (0, 0), "xy": (0, 1), "yy": (1, 1)}.items():
    glue(f"E_{key}", round(strain[0][r, c], 4), display=False)
```

The strain of every atom is the applied strain $\tfrac{1}{2}(\mathbf{F}^\mathsf{T}\mathbf{F} - \mathbf{I})$, with $E_{xx}$ = {glue}`E_xx`, $E_{xy}$ = {glue}`E_xy` and $E_{yy}$ = {glue}`E_yy`.
$D^2_\mathrm{min}$ and the slip vector are zero to rounding.

```{code-cell} ipython3
:tags: [remove-cell]
assert np.allclose(strain, 0.5 * (F.T @ F - np.eye(3)), atol=1e-10)
assert np.abs(d2).max() < 1e-12 and np.abs(slip).max() < 1e-12
assert "pyscal_strain" not in reference.arrays
```

The first argument is the current structure and the second the reference.
The two must contain the same atoms in the same order.
Each function returns its result and stores it on the current structure:

| Function | Key | Shape | Content |
|---|---|---|---|
| `atomic_strain` | `atoms.arrays["pyscal_strain"]` | $(N, 3, 3)$ | $\mathbf{E}$ |
| `von_mises_strain` | `atoms.arrays["pyscal_von_mises"]` | $(N,)$ | $\eta$ |
| `d2min` | `atoms.arrays["pyscal_d2min"]` | $(N,)$ | $D^2_\mathrm{min}$  (Å²) |
| `slip_vector` | `atoms.arrays["pyscal_slip_vector"]` | $(N, 3)$ | $\mathbf{s}$  (Å) |

`von_mises_strain` also stores `pyscal_strain`.
The reference structure is not changed.

`find_neighbors` must be called on both structures, otherwise the functions raise an error.
The pairs of atom $i$ are the atoms that appear in both of its neighbor lists, so the two lists do not have to be identical.
Any neighbor method works.
Use the same method for both structures, and a cutoff large enough that the reference neighbors are still neighbors in the deformed structure (see [Slip between two blocks](#slip-between-two-blocks)).

## Strain at finite temperature

In a simulation, both the reference and the deformed structure contain thermal vibrations.
We model them by random displacements with a standard deviation $\sigma$ in each direction, drawn independently for the two structures, and apply the deformation of the previous section.
The strain is computed with three neighbor lists: the first shell (12 neighbors), the first two shells (18) and the first three shells (42).

```{code-cell} ipython3
a = 3.61
perfect = bulk("Cu", "fcc", a=a, cubic=True).repeat(6)
exact = 0.5 * (F.T @ F - np.eye(3))

shells = {"1 shell": 3.1, "2 shells": 4.0, "3 shells": 4.7}
sigmas = [0.0, 0.025, 0.05, 0.075, 0.1]

results = {}
for sigma in sigmas:
    reference = perfect.copy()
    deformed = perfect.copy()
    deformed.set_cell(perfect.cell @ F.T, scale_atoms=True)
    reference.rattle(sigma, seed=1)
    deformed.rattle(sigma, seed=2)
    for name, cutoff in shells.items():
        for structure in (reference, deformed):
            pyscal.find_neighbors(structure, method="cutoff", cutoff=cutoff)
        strain = pyscal.atomic_strain(deformed, reference)
        results[sigma, name] = (strain, pyscal.d2min(deformed, reference))
```

```{code-cell} ipython3
:tags: [remove-cell]
from _plotstyle import RED, TEAL, PURPLE
SHELL_STYLE = {"1 shell": (RED, "o"), "2 shells": (TEAL, "s"), "3 shells": (PURPLE, "D")}
nearest = a / np.sqrt(2)
spread = {(s, n): results[s, n][0][:, 0, 1].std() for s in sigmas for n in shells}
glue("sigma_hist", 0.05, display=False)
glue("sigma_hist_pct", round(100 * 0.05 / nearest), display=False)
glue("sigma_max", 0.1, display=False)
glue("spread_1", round(spread[0.05, "1 shell"], 4), display=False)
glue("spread_3", round(spread[0.05, "3 shells"], 4), display=False)
glue("mean_exy_1", round(results[0.05, "1 shell"][0][:, 0, 1].mean(), 4), display=False)
glue("mean_exy_3", round(results[0.05, "3 shells"][0][:, 0, 1].mean(), 4), display=False)
glue("exx_exact", round(exact[0, 0], 4), display=False)
glue("exx_noisy_1", round(results[0.1, "1 shell"][0][:, 0, 0].mean(), 4), display=False)
glue("exx_noisy_3", round(results[0.1, "3 shells"][0][:, 0, 0].mean(), 4), display=False)
glue("d2_thermal_1", round(results[0.05, "1 shell"][1].mean(), 3), display=False)
glue("d2_thermal_3", round(results[0.05, "3 shells"][1].mean(), 3), display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(columns=2, ratio=0.42, wspace=0.3)
ax = mp[0, 0]
bins = np.linspace(-0.02, 0.06, 41)
for name, (colour, marker) in SHELL_STYLE.items():
    ax.hist(results[0.05, name][0][:, 0, 1], bins=bins, histtype="step", lw=1.6,
            color=colour, density=True, label=name)
ax.axvline(exact[0, 1], ls="--", color=DARK, lw=1)
note(ax, f"$\\sigma$ = 0.05 Å", loc="upper right")
label(ax, "(a)")
ax.set_xlabel("$E_{xy}$")
ax.set_ylabel("Probability density")
ax.set_ylim(0, 140)

ax = mp[0, 1]
for name, (colour, marker) in SHELL_STYLE.items():
    ax.plot(sigmas, [spread[s, name] for s in sigmas], marker=marker, color=colour,
            mec=DARK, mew=0.7, ms=5, lw=1.6, label=name)
label(ax, "(b)")
ax.set_xlabel(r"$\sigma$  (Å)")
ax.set_ylabel("Standard deviation of $E_{xy}$")
ax.set_ylim(0, None)
handles, labels = mp[0, 1].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

Panel (a) shows the distribution of $E_{xy}$ over the atoms for $\sigma$ = {glue}`sigma_hist` Å, about {glue}`sigma_hist_pct` % of the nearest neighbor distance.
The dashed line is the applied shear strain.
The mean of $E_{xy}$ over the atoms is {glue}`mean_exy_1` with the first shell and {glue}`mean_exy_3` with three shells, equal to the applied value.
The strain of a single atom scatters around the mean, with a standard deviation of {glue}`spread_1` with the first shell and {glue}`spread_3` with three shells.
Panel (b) shows that the scatter grows linearly with $\sigma$.
A fit to more neighbors averages out more of the thermal noise, at the cost of a coarser spatial resolution.

The thermal noise also shifts the diagonal components.
At $\sigma$ = {glue}`sigma_max` Å, the mean $E_{xx}$ is {glue}`exx_noisy_1` with the first shell and {glue}`exx_noisy_3` with three shells, instead of {glue}`exx_exact`.
The noise in the reference vectors $\mathbf{R}_{ij}$ biases the least squares fit towards smaller $\mathbf{F}$.
The bias is smaller with more shells, whose neighbor vectors are longer.

$D^2_\mathrm{min}$ is not zero in a crystal at finite temperature.
At $\sigma$ = {glue}`sigma_hist` Å, its mean is {glue}`d2_thermal_1` Å² with the first shell and {glue}`d2_thermal_3` Å² with three shells.
This thermal background is the level to compare with when looking for plastic events.

## Slip between two blocks

When a dislocation passes through a crystal, the part above its slip plane is displaced by the Burgers vector $\mathbf{b}$ relative to the part below.
We model this with an fcc crystal whose (111) planes are horizontal, and displace its upper half rigidly by $\mathbf{b}$.
The cell vector along $z$ is tilted by the same $\mathbf{b}$, so that the crystal stays continuous across the periodic boundary and there is only one slip plane.
We compare two slip vectors in the (111) plane:

- a full slip by $\mathbf{b} = \tfrac{a}{2}[1\bar{1}0]$, with $|\mathbf{b}| = a/\sqrt{2}$, after which the crystal is perfect again,
- a slip by the Shockley partial $\mathbf{b} = \tfrac{a}{6}[\bar{1}\bar{1}2]$, with $|\mathbf{b}| = a/\sqrt{6}$, which leaves an intrinsic stacking fault.

Both structures get random displacements with $\sigma$ = 0.05 Å as a rough model of thermal vibrations.

```{code-cell} ipython3
from ase.build import fcc111

a = 3.61
layers = 18
crystal = fcc111("Cu", size=(6, 4, layers), a=a, orthogonal=True, periodic=True)
spacing = crystal.cell[2, 2] / layers
layer = np.round(crystal.positions[:, 2] / spacing).astype(int)
upper = layer >= layers // 2

# in this cell, x is along [1-10] and y along [11-2]
burgers = {
    "full slip": np.array([a / np.sqrt(2), 0.0, 0.0]),
    "stacking fault": np.array([0.0, -a / np.sqrt(6), 0.0]),
}

slipped = {}
for name, b in burgers.items():
    reference = crystal.copy()
    current = crystal.copy()
    current.positions[upper] += b
    current.cell[2] += b
    reference.rattle(0.05, seed=1)
    current.rattle(0.05, seed=2)
    pyscal.find_neighbors(reference, method="cutoff", cutoff=3.1)
    pyscal.find_neighbors(current, method="cutoff", cutoff=4.5)
    s = pyscal.slip_vector(current, reference)
    d2 = pyscal.d2min(current, reference)
    cna = pyscal.common_neighbor_analysis(current)
    slipped[name] = (current, s, d2, cna)
```

The reference neighbors are the first shell.
The current structure uses a larger cutoff, so that the neighbors that moved with the upper block are still paired.

```{code-cell} ipython3
:tags: [remove-cell]
plane = (layer == layers // 2 - 1) | (layer == layers // 2)
bulk_layers = ~plane
for name, b in burgers.items():
    current, s, d2, cna = slipped[name]
    key = name.replace(" ", "_")
    glue(f"s_plane_{key}", round(np.linalg.norm(s[plane], axis=1).mean(), 2), display=False)
    glue(f"s_expected_{key}", round(np.linalg.norm(b) / 4, 2), display=False)
    glue(f"d2_plane_{key}", round(d2[plane].mean(), 2), display=False)
    glue(f"cna_fcc_{key}", cna["fcc"], display=False)
    glue(f"cna_hcp_{key}", cna["hcp"], display=False)
s_all = np.concatenate([np.linalg.norm(slipped[n][1][bulk_layers], axis=1) for n in burgers])
d2_all = np.concatenate([slipped[n][2][bulk_layers] for n in burgers])
glue("s_bulk", round(s_all.mean(), 2), display=False)
glue("d2_bulk", round(d2_all.mean(), 3), display=False)
glue("n_atoms_slip", len(crystal), display=False)
glue("n_per_layer", int(np.sum(layer == 0)), display=False)
```

```{code-cell} ipython3
:tags: [hide-input]
SLIP_STYLE = {"full slip": (RED, "o"), "stacking fault": (TEAL, "s")}
z = crystal.positions[:, 2]
mp = figure(columns=2, ratio=0.42, wspace=0.3)
for k, (title, index, unit) in enumerate([("(a)  slip vector", 1, r"$|\mathbf{s}|$  (Å)"),
                                          ("(b)  $D^2_\mathrm{min}$", 2, r"$D^2_\mathrm{min}$  (Å$^2$)")]):
    ax = mp[0, k]
    for name, (colour, marker) in SLIP_STYLE.items():
        values = slipped[name][index]
        if values.ndim == 2:
            values = np.linalg.norm(values, axis=1)
        ax.scatter(z, values, s=8, color=colour, alpha=0.35, lw=0)
        means = [values[layer == l].mean() for l in range(layers)]
        ax.plot(np.arange(layers) * spacing, means, marker=marker, color=colour,
                mec=DARK, mew=0.7, ms=4.5, lw=1.4, label=name)
    ax.axvline((layers // 2 - 0.5) * spacing, ls="--", color=DARK, lw=1)
    label(ax, title)
    ax.set_xlabel("$z$  along [111]  (Å)")
    ax.set_ylabel(unit)
    ax.set_ylim(0, None)
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=2, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

The figure shows $|\mathbf{s}|$ and $D^2_\mathrm{min}$ for each of the {glue}`n_atoms_slip` atoms (points) and the mean of each (111) layer of {glue}`n_per_layer` atoms (lines).
The dashed line marks the slip plane.
Both descriptors are large only in the two layers next to the slip plane.
Each atom in these layers has 3 of its 12 nearest neighbors on the other side of the plane, so $|\mathbf{s}| \approx |\mathbf{b}|/4$.
The mean is {glue}`s_plane_full_slip` Å for the full slip ($|\mathbf{b}|/4$ = {glue}`s_expected_full_slip` Å) and {glue}`s_plane_stacking_fault` Å for the stacking fault ({glue}`s_expected_stacking_fault` Å).
Away from the plane, the two curves coincide because both structures have the same random displacements.
There, the thermal displacements give a background of {glue}`s_bulk` Å for $|\mathbf{s}|$ and {glue}`d2_bulk` Å² for $D^2_\mathrm{min}$.
The slip vector contains the displacement of the atom itself relative to its neighbors, so its thermal background is large compared to that of $D^2_\mathrm{min}$.

The two descriptors see both slips, but CNA sees only one of them.
After the slip by the partial, CNA labels {glue}`cna_hcp_stacking_fault` atoms hcp, the two layers of the stacking fault.
After the full slip, the crystal is perfect again, and CNA labels {glue}`cna_fcc_full_slip` of the {glue}`n_atoms_slip` atoms fcc.
The slip vector and $D^2_\mathrm{min}$ still show where the crystal has slipped, because they compare with the reference.

The pairs of neighbors depend on the cutoff of the current structure.
We repeat the full slip without random displacements and with the first shell cutoff for both structures:

```{code-cell} ipython3
reference = crystal.copy()
current = crystal.copy()
current.positions[upper] += burgers["full slip"]
current.cell[2] += burgers["full slip"]
for structure in (reference, current):
    pyscal.find_neighbors(structure, method="cutoff", cutoff=3.1)
s = pyscal.slip_vector(current, reference)
np.linalg.norm(s[plane], axis=1).mean()
```

```{code-cell} ipython3
:tags: [remove-cell]
glue("s_plane_same_cutoff", round(np.linalg.norm(s[plane], axis=1).mean(), 2), display=False)
assert np.allclose(np.linalg.norm(s[plane], axis=1), np.linalg.norm(burgers["full slip"]) / 10)
```

Now $|\mathbf{s}|$ on the slip plane is {glue}`s_plane_same_cutoff` Å instead of $|\mathbf{b}|/4$ = {glue}`s_expected_full_slip` Å.
Two of the three neighbors across the plane are no longer within the cutoff after the slip, so they are not paired, and $|\mathbf{s}| = |\mathbf{b}|/10$.

```{code-cell} ipython3
:tags: [remove-cell]
# von Mises strain of the same pure shear, along the axes and rotated by 45 degrees about z
reference = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(5)
shear = np.array([[1.0, 0.02, 0.0], [0.02, 1.0, 0.0], [0.0, 0.0, 1.0]])
c = np.cos(np.pi / 4)
rotation = np.array([[c, -c, 0.0], [c, c, 0.0], [0.0, 0.0, 1.0]])
eta_shear = []
for G in (shear, rotation.T @ shear @ rotation):
    current = reference.copy()
    current.set_cell(reference.cell @ G.T, scale_atoms=True)
    for structure in (reference, current):
        pyscal.find_neighbors(structure, method="cutoff", cutoff=3.1)
    eta_shear.append(pyscal.von_mises_strain(current, reference)[0])
glue("eta_axes", round(eta_shear[0], 3), display=False)
glue("eta_rotated", round(eta_shear[1], 3), display=False)
```

## Things to watch

- **Same atoms, same order.** Atom $i$ of the current structure is compared with atom $i$ of the reference. If the atoms are in a different order, the results are meaningless. ASE sorts the atoms of a LAMMPS dump file by their id when it reads the file.
- **Neighbors that leave the cutoff.** Only neighbors found in both structures are paired. For large deformations or slip, use a larger cutoff for the current structure, as in the example above.
- **Thermal background.** At finite temperature, $D^2_\mathrm{min}$, $|\mathbf{s}|$ and the scatter of the strain are not zero. Compare values with those of a region known to be undeformed, or average over time.
- **Normalisation of $D^2_\mathrm{min}$.** pyscal divides by the number of pairs. Values from codes that use the sum of Falk and Langer are larger by a factor $n_i$.
- **Definition of $\eta$.** The differences of the normal strains enter $\eta$ with a factor 1/2, as in the equation above. Shimizu et al. [2] use a factor 1/6, which makes $\eta$ independent of the orientation of the axes. With the factor 1/2, $\eta$ depends on the orientation: the same pure shear gives $\eta$ = {glue}`eta_axes` along the cube axes and {glue}`eta_rotated` after a rotation by 45° about $z$. Compare values of $\eta$ only for the same orientation, or compute an invariant measure from `pyscal_strain`.

## References

1. M. L. Falk and J. S. Langer, Dynamics of viscoplastic deformation in amorphous solids, *Phys. Rev. E* **57**, 7192 (1998). [doi:10.1103/PhysRevE.57.7192](https://doi.org/10.1103/PhysRevE.57.7192)
2. F. Shimizu, S. Ogata and J. Li, Theory of shear banding in metallic glasses and molecular dynamics calculations, *Mater. Trans.* **48**, 2923 (2007). [doi:10.2320/matertrans.MJ200769](https://doi.org/10.2320/matertrans.MJ200769)
3. J. A. Zimmerman, C. L. Kelchner, P. A. Klein, J. C. Hamilton and S. M. Foiles, Surface step effects on nanoindentation, *Phys. Rev. Lett.* **87**, 165507 (2001). [doi:10.1103/PhysRevLett.87.165507](https://doi.org/10.1103/PhysRevLett.87.165507)
