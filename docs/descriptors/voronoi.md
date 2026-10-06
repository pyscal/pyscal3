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

# Voronoi tessellation

The Voronoi tessellation divides space into one cell per atom, the region that is closer to that atom than to any other.
pyscal uses it to find neighbors without a cutoff, each with a weight (see [Finding neighbors](../guide/neighbors)), to measure the volume available to each atom, and to compute the Voronoi structure vector [1, 2], which identifies the local structure from the shape of the cell.
The face weights are also the basis of the [Minkowski structure metrics](minkowski).

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

The Voronoi cell of atom $i$ at position $\mathbf{r}_i$ is

$$
\mathcal{V}_i = \left\{ \mathbf{x} : |\mathbf{x} - \mathbf{r}_i| \le |\mathbf{x} - \mathbf{r}_j| \ \text{for all atoms}\ j \right\},
$$

where the atoms $j$ include the periodic images.
The cells fill the simulation cell without overlap, so their volumes $V_i$ add up to the volume of the simulation cell.
Each face of $\mathcal{V}_i$ lies on the plane halfway between $i$ and another atom $j$, and these atoms $j$ are the Voronoi neighbors of $i$.
The area $A_{ij}$ of the shared face gives the weight of the neighbor:

$$
w_{ij} = \frac{A_{ij}^{\,p}}{\sum_{k} A_{ik}^{\,p}},
$$

where the sum runs over all neighbors $k$ of $i$ and the exponent $p$ is set by `voroexp` (default 1).

### Voronoi structure vector

The faces of a Voronoi cell are polygons.
The Voronoi structure vector [1, 2] counts them by their number of edges:

$$
\langle n_3\ n_4\ n_5\ n_6 \rangle,
$$

where $n_k$ is the number of faces with $k$ edges.
At finite temperature, thermal motion creates small faces and short edges that are not part of the ideal cell.
`voronoi_vector` ignores a face if its area is not larger than `area_cutoff` (default 0.01) times the surface of the Voronoi cell, and an edge if its length is not larger than `edge_cutoff` (default 0.05) times the perimeter of the face.
Faces that are left with fewer than three or more than six edges are not counted.

## Usage

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="voronoi")
vectors = pyscal.voronoi_vector(atoms)
vectors[:4]
```

`find_neighbors` with `method="voronoi"` computes the tessellation with [Voro++](https://math.lbl.gov/voro++/) [3].
`voronoi_vector` needs this tessellation and raises an error if the neighbors were found with another method.
It returns an integer array with one row $\langle n_3\ n_4\ n_5\ n_6 \rangle$ per atom.
The results are stored on `atoms`:

| Key | Content |
|---|---|
| `atoms.arrays["pyscal_vorovector"]` | $\langle n_3\ n_4\ n_5\ n_6 \rangle$, shape $(N, 4)$, from `voronoi_vector` |
| `atoms.arrays["pyscal_voronoi_volume"]` | volume $V_i$ of the Voronoi cell, shape $(N,)$ |
| `atoms.info["pyscal_bond_weight"]` | weight $w_{ij}$ of each neighbor, in the flat neighbor arrays (see [Finding neighbors](../guide/neighbors)) |
| `pyscal_face_vertices` | number of vertices of each face, in the order of the neighbors |
| `pyscal_face_perimeters` | perimeter of each face |
| `pyscal_vertex_vectors` | coordinates $x, y, z$ of all vertices of the cell, relative to the atom, one after the other |
| `pyscal_vertex_numbers` | for each face, its number of vertices followed by their indices in `pyscal_vertex_vectors` |
| `pyscal_vertex_positions` | positions of the vertices of the cell |

The last five keys hold one row per atom.
When every atom has the same number of faces and vertices, as in a perfect crystal, they are arrays in `atoms.arrays`.
Otherwise they are lists in `atoms.info`, and `pyscal_vertex_positions` is always in `atoms.info`.

The volumes of the cells add up to the volume of the simulation cell:

```{code-cell} ipython3
atoms.arrays["pyscal_voronoi_volume"].sum(), atoms.get_volume()
```

With a positive `cutoff`, `find_neighbors` also merges the vertices of all cells that are closer than `cutoff` and stores one position for each group in `atoms.info["pyscal_unique_vertices"]`.
The vertices of the Voronoi cells are the centres of the interstitial sites of a crystal.
In a perfect fcc crystal there are three per atom, one octahedral and two tetrahedral sites:

```{code-cell} ipython3
from pyscal.structures import make_crystal

fcc = make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4))
pyscal.find_neighbors(fcc, method="voronoi", cutoff=0.1)
len(fcc.info["pyscal_unique_vertices"]) / len(fcc)
```

## Values for perfect crystals

```{code-cell} ipython3
import numpy as np
import pandas as pd
from ase.cluster import Icosahedron

crystals = {
    "fcc": make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4)),
    "hcp": make_crystal("hcp", lattice_constant=2.51, repetitions=(4, 4, 4)),
    "bcc": make_crystal("bcc", lattice_constant=2.87, repetitions=(4, 4, 4)),
}
perfect = {}
for name, crystal in crystals.items():
    pyscal.find_neighbors(crystal, method="voronoi")
    perfect[name] = pyscal.voronoi_vector(crystal)[0]

cluster = Icosahedron("Cu", noshells=2)
cluster.center(vacuum=8)
cluster.pbc = True
pyscal.find_neighbors(cluster, method="voronoi")
centre = np.argmin(np.linalg.norm(cluster.positions - cluster.get_center_of_mass(), axis=1))
perfect["ico"] = pyscal.voronoi_vector(cluster)[centre]

pd.DataFrame(perfect, index=["n3", "n4", "n5", "n6"]).T
```

The table shows the vector of one atom of each perfect crystal and of the central atom of a 13 atom icosahedral cluster.
The Voronoi cell of fcc is a rhombic dodecahedron with 12 faces of four edges, and that of bcc a truncated octahedron with 6 square and 8 hexagonal faces.
The cell of an atom at the centre of an icosahedron is a pentagonal dodecahedron, $\langle 0\ 0\ 12\ 0 \rangle$.
The cell of hcp also has 12 faces of four edges, so **the Voronoi vector does not distinguish fcc from hcp**.

## Crystals at finite temperature

How well does the vector identify a crystal when the atoms vibrate?
We first add random displacements to perfect fcc, hcp and bcc crystals, a rough model of thermal vibrations, and count the atoms that keep the vector of the perfect crystal.
The standard deviation of the displacements in each direction is given in percent of the nearest neighbor distance $d_1$.

```{code-cell} ipython3
from ase.build import bulk

nearest = {"fcc": 3.61 / np.sqrt(2), "hcp": 2.55, "bcc": 2.87 * np.sqrt(3) / 2}
noises = np.linspace(0, 0.08, 17)

def build(name):
    if name == "hcp":
        return bulk("Cu", "hcp", a=2.55).repeat((8, 8, 5))
    a = 3.61 if name == "fcc" else 2.87
    return bulk("Cu", name, a=a, cubic=True).repeat(6)

def fraction_ideal(atoms, name, **thresholds):
    vectors = pyscal.voronoi_vector(atoms, **thresholds)
    return np.mean(np.all(vectors == perfect[name], axis=1))

noise_scan = {}
for name in ("fcc", "hcp", "bcc"):
    noise_scan[name] = []
    for noise in noises:
        crystal = build(name)
        crystal.rattle(noise * nearest[name], seed=1)
        pyscal.find_neighbors(crystal, method="voronoi")
        noise_scan[name].append(fraction_ideal(crystal, name))
```

We then take the fcc and bcc MD snapshots from the `examples` folder and vary `edge_cutoff`, which sets how short an edge has to be to be ignored.

```{code-cell} ipython3
snapshots = {
    "fcc": read("conf.fcc.dump", format="lammps-dump-text"),
    "bcc": read("conf.bcc.dump", format="lammps-dump-text"),
}
edge_cutoffs = np.linspace(0, 0.15, 31)
edge_scan = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="voronoi")
    edge_scan[name] = [fraction_ideal(snapshot, name, edge_cutoff=e) for e in edge_cutoffs]
```

```{code-cell} ipython3
:tags: [hide-input]
markers = {"fcc": "o", "hcp": "s", "bcc": "D"}
mp = figure(columns=2, ratio=0.42, wspace=0.15)
ax = mp[0, 0]
for name, values in noise_scan.items():
    ax.plot(100 * noises, values, marker=markers[name], color=COLOURS[name], mec=DARK,
            mew=0.7, ms=5, lw=1.6, label=name)
label(ax, "(a)  displaced crystals")
ax.set_xlabel("Standard deviation  (% of $d_1$)")
ax.set_ylabel("Fraction with the ideal vector")
ax.set_ylim(-0.03, 1.05)

ax = mp[0, 1]
for name, values in edge_scan.items():
    ax.plot(edge_cutoffs, values, marker=markers[name], markevery=2, color=COLOURS[name],
            mec=DARK, mew=0.7, ms=5, lw=1.6)
ax.axvline(0.05, ls=":", color=DARK, lw=0.8)
ax.text(0.053, 0.97, "default", transform=ax.get_xaxis_transform(), ha="left", va="top",
        fontsize=9, color=DARK)
label(ax, "(b)  MD snapshots")
ax.set_xlabel("edge_cutoff")
ax.set_ylim(-0.03, 1.05)
ax.tick_params(labelleft=False)
handles, labels = mp[0, 0].get_legend_handles_labels()
mp.fig.legend(handles, labels, frameon=False, ncol=3, loc="upper center",
              bbox_to_anchor=(0.5, -0.1));
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue

def at_noise(name, noise):
    return round(100 * noise_scan[name][int(np.argmin(np.abs(noises - noise)))])

glue("fcc_at_4", at_noise("fcc", 0.04), display=False)
glue("hcp_at_4", at_noise("hcp", 0.04), display=False)
glue("bcc_at_6", at_noise("bcc", 0.06), display=False)
glue("bcc_at_8", at_noise("bcc", 0.08), display=False)

default = int(np.argmin(np.abs(edge_cutoffs - 0.05)))
glue("fcc_md_default", round(100 * edge_scan["fcc"][default], 1), display=False)
glue("bcc_md_default", round(100 * edge_scan["bcc"][default]), display=False)
glue("fcc_md_best", round(100 * max(edge_scan["fcc"])), display=False)
glue("fcc_md_best_cutoff", round(float(edge_cutoffs[np.argmax(edge_scan["fcc"])]), 3), display=False)
glue("bcc_md_zero", round(100 * edge_scan["bcc"][0]), display=False)

# best common setting of both thresholds
grid = [(e, a) for e in np.linspace(0, 0.15, 16) for a in np.linspace(0, 0.05, 11)]
common = max(min(fraction_ideal(snapshots[n], n, edge_cutoff=e, area_cutoff=a)
                 for n in snapshots) for e, a in grid)
glue("best_common", round(100 * common), display=False)
```

Panel (a) shows that the vector of bcc survives much larger displacements than those of fcc and hcp.
At a standard deviation of 4 %, only {glue}`fcc_at_4` % of the fcc atoms and {glue}`hcp_at_4` % of the hcp atoms keep the ideal vector, while bcc still has {glue}`bcc_at_6` % at 6 % and {glue}`bcc_at_8` % at 8 %.
The reason is the shape of the cells.
In the rhombic dodecahedron of fcc and in the cell of hcp, four faces meet at six of the vertices.
These vertices are shared by six cells, and almost any displacement splits each of them into several vertices joined by short new edges, and creates small new faces.
In the truncated octahedron of bcc, exactly three faces meet at every vertex, and small displacements only move the vertices.

In the MD snapshots, panel (b), the default `edge_cutoff` gives the ideal vector to {glue}`bcc_md_default` % of the bcc atoms but only {glue}`fcc_md_default` % of the fcc atoms.
Ignoring longer edges helps fcc, up to {glue}`fcc_md_best` % at `edge_cutoff` = {glue}`fcc_md_best_cutoff`, but at that value the edges of the hexagons of bcc are ignored as well.
With the default `area_cutoff`, bcc does best without an edge threshold ({glue}`bcc_md_zero` %).
We also scanned `area_cutoff` from 0 to 0.05 together with `edge_cutoff`.
No pair of values gives the ideal vector to more than {glue}`best_common` % of the atoms in both crystals.
At finite temperature, the Voronoi vector identifies most bcc atoms, but few fcc or hcp atoms.
For close packed structures, use [common neighbor analysis](cna) or [averaged Steinhardt parameters](steinhardt).

## Voronoi volume

The volume of the Voronoi cell measures the space available to each atom.
It is used, for example, to measure the local density or the free volume in a glass.
We compare the volumes in the fcc, bcc and liquid snapshots.

```{code-cell} ipython3
snapshots["liquid"] = read("conf.lqd.Al.dump", format="lammps-dump-text")
volumes = {}
for name, snapshot in snapshots.items():
    pyscal.find_neighbors(snapshot, method="voronoi")
    volumes[name] = snapshot.arrays["pyscal_voronoi_volume"]
```

```{code-cell} ipython3
:tags: [hide-input]
mp = figure(ratio=0.42)
ax = mp[0, 0]
bins = np.linspace(13, 31, 73)
for name, values in volumes.items():
    ax.hist(values, bins=bins, density=True, color=COLOURS[name], alpha=0.7, ec=DARK,
            lw=0.5, label=name)
ax.set_xticks(range(14, 32, 2))
ax.set_xlabel("Voronoi volume  $V_i$  (Å$^3$)")
ax.set_ylabel("Probability density")
ax.legend(frameon=False);
```

```{code-cell} ipython3
:tags: [remove-cell]
for name, values in volumes.items():
    glue(f"{name}_volume_mean", round(float(values.mean()), 1), display=False)
    glue(f"{name}_volume_spread", round(float(100 * values.std() / values.mean()), 1), display=False)
```

In the crystals, the volumes scatter around their mean, {glue}`fcc_volume_mean` Å$^3$ in fcc and {glue}`bcc_volume_mean` Å$^3$ in bcc, with a standard deviation of {glue}`fcc_volume_spread` % and {glue}`bcc_volume_spread` % of the mean.
In the liquid, the mean is {glue}`liquid_volume_mean` Å$^3$ and the standard deviation is {glue}`liquid_volume_spread` % of the mean.
The broad distribution reflects the disorder of the liquid, where some atoms have much more space than others.

## Things to watch

```{code-cell} ipython3
:tags: [remove-cell]
liquid_faces = np.concatenate([np.ravel(f) for f in snapshots["liquid"].info["pyscal_face_vertices"]])
glue("liquid_large_faces", round(100 * np.mean(liquid_faces > 6)), display=False)

```

- **fcc and hcp have the same vector.** Both are $\langle 0\ 12\ 0\ 0 \rangle$. Use another descriptor to tell them apart.
- **Thresholds.** The vectors depend on `edge_cutoff` and `area_cutoff`, and no single setting works for all crystal structures at finite temperature. Report the values used.
- **Faces with more than six edges.** These are not counted in the vector. In the liquid snapshot, {glue}`liquid_large_faces` % of the faces have more than six vertices.
- **More neighbors than expected.** The small faces created by thermal motion make extra Voronoi neighbors, about 14 per atom in an fcc crystal at finite temperature (see [Finding neighbors](../guide/neighbors)). Descriptors that weight neighbors by the face area, such as the [Minkowski structure metrics](minkowski), are less affected than those that count them.

## References

1. J. L. Finney, Random packings and the structure of simple liquids. I. The geometry of random close packing, *Proc. R. Soc. Lond. A* **319**, 479 (1970). [doi:10.1098/rspa.1970.0189](https://doi.org/10.1098/rspa.1970.0189)
2. M. Tanemura, Y. Hiwatari, H. Matsuda, T. Ogawa, N. Ogita and A. Ueda, Geometrical analysis of crystallization of the soft-core model, *Prog. Theor. Phys.* **58**, 1079 (1977). [doi:10.1143/PTP.58.1079](https://doi.org/10.1143/PTP.58.1079)
3. C. H. Rycroft, VORO++: A three-dimensional Voronoi cell library in C++, *Chaos* **19**, 041111 (2009). [doi:10.1063/1.3215722](https://doi.org/10.1063/1.3215722)
