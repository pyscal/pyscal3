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

# Building structures

pyscal analyses ASE `Atoms` objects, so structures can come from three places:
a file written by a simulation code (see [Reading files and trajectories](files)),
the builders of ASE,
or the builders of pyscal in `pyscal.structures`.
This page shows the last two.

```{code-cell} ipython3
:tags: [remove-cell]
import os, sys, warnings
sys.path.insert(0, os.path.abspath(".."))
from _plotstyle import EXAMPLES, COLOURS, DARK, figure, label, note
os.chdir(EXAMPLES)
warnings.simplefilter("ignore")
%config InlineBackend.figure_format = "retina"
```

## ASE builders

[`ase.build`](https://wiki.fysik.dtu.dk/ase/ase/build/build.html) creates crystals, surfaces and molecules.
`bulk` gives the primitive cell by default and the conventional cubic cell with `cubic=True`.
`repeat` makes a supercell.

```{code-cell} ipython3
from ase.build import bulk

atoms = bulk("Cu", "fcc", a=3.61, cubic=True).repeat((4, 4, 4))
len(atoms), atoms.cell.lengths()
```

## pyscal builders

`pyscal.structures` creates crystals from a structure name or an element, custom lattices, and grain boundaries.

```{code-cell} ipython3
from pyscal.structures import make_crystal, available_structures

available_structures()
```

`make_crystal` takes the structure, the lattice constant and the number of unit cells in each direction:

```{code-cell} ipython3
fcc = make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4), element="Cu")
hcp = make_crystal("hcp", lattice_constant=2.51, repetitions=(4, 4, 4), ca_ratio=1.633)
l12 = make_crystal("l12", lattice_constant=3.57, repetitions=(3, 3, 3), element=["Al", "Ni"])
fcc.get_chemical_formula(), hcp.get_chemical_formula(), l12.get_chemical_formula()
```

- `element` gives the chemical element, or a list of elements for structures with several sublattices such as `l12` and `b2`. Without it, the atoms of sublattice $k$ get the atomic number $k$, so that sublattices can be told apart.
- `ca_ratio` sets $c/a$ for `hcp` and `dhcp`.
- `primitive=True` gives the primitive cell.

For random displacements, use `atoms.rattle(stdev, seed=...)` from ASE on the finished structure, as below.

`make_element` looks up a tabulated structure and lattice constant for an element:

```{code-cell} ipython3
from pyscal.structures import make_element

iron = make_element("Fe", repetitions=(4, 4, 4))
iron.get_chemical_formula(), iron.cell.lengths()
```

`make_general_lattice` builds a lattice from fractional positions, types and a box, for structures that are not in the list.

## Grain boundaries

`make_grain_boundary` builds a bicrystal with a symmetric tilt grain boundary, from the tilt axis, the $\Sigma$ value of the coincidence site lattice and the boundary plane.
The cell is periodic, so it contains two grain boundaries.
Common neighbor analysis finds the atoms at the boundaries, which are not labelled fcc:

```{code-cell} ipython3
import pyscal
from pyscal.structures import make_grain_boundary

bicrystal = make_grain_boundary(axis=[1, 0, 0], sigma=5, gb_plane=[0, 1, 3],
                                structure="fcc", lattice_constant=4.05,
                                repetitions=(2, 1, 3))
pyscal.common_neighbor_analysis(bicrystal)
```

```{code-cell} ipython3
:tags: [hide-input]
import numpy as np

names = {0: "others", 1: "fcc", 2: "hcp", 3: "bcc", 4: "ico"}
labels = bicrystal.arrays["pyscal_structure"]
mp = figure(ratio=0.45)
ax = mp[0, 0]
for code in np.unique(labels):
    sel = labels == code
    ax.scatter(bicrystal.positions[sel, 0], bicrystal.positions[sel, 2], s=14,
               color=COLOURS[names[code]], ec=DARK, lw=0.3, label=names[code])
ax.set_xlabel("$x$  (Å)")
ax.set_ylabel("$z$  (Å)")
ax.set_aspect("equal")
mp.fig.legend(*ax.get_legend_handles_labels(), frameon=False, ncol=2,
              loc="upper center", bbox_to_anchor=(0.5, -0.02));
```

The figure shows the atoms projected along the tilt axis.
The grey atoms are the two grain boundaries, one in the middle of the cell and one at its periodic edge.
The blue atoms are the two grains.

## Changing a structure

ASE `Atoms` objects can be changed with plain NumPy operations, for example to create defects:

```{code-cell} ipython3
crystal = make_crystal("fcc", lattice_constant=3.61, repetitions=(4, 4, 4), element="Cu")

del crystal[0]                                  # a vacancy
crystal.rattle(stdev=0.05, seed=1)              # random displacements, 0.05 Å
crystal.set_cell(crystal.cell * 1.01, scale_atoms=True)   # 1 % expansion
len(crystal)
```

[Common neighbor analysis](../descriptors/cna) builds a crystal with stacking faults by removing close-packed layers, and [Wigner–Seitz analysis](../descriptors/wigner_seitz) creates vacancies and interstitials.
