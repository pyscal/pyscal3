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

# Reading files and trajectories

pyscal reads structures through ASE, so every [file format that ASE supports](https://wiki.fysik.dtu.dk/ase/ase/io/io.html) can be analysed.
For long LAMMPS trajectories, `pyscal.Trajectory` reads one frame at a time.
This page also shows how to write the results back to a file.

```{code-cell} ipython3
:tags: [remove-cell]
import os, sys, warnings
sys.path.insert(0, os.path.abspath(".."))
from _plotstyle import EXAMPLES, COLOURS, DARK, figure, label, note
os.chdir(EXAMPLES)
warnings.simplefilter("ignore")
%config InlineBackend.figure_format = "retina"
```

## Single structures

```python
from ase.io import read

atoms = read("dump.lammpstrj", format="lammps-dump-text")   # LAMMPS dump, last frame
atoms = read("POSCAR")                                        # VASP
atoms = read("structure.cif")                                 # CIF
atoms = read("structure.extxyz")                              # extended XYZ
```

For LAMMPS dump files, ASE sets the chemical elements from the `element` column if there is one, otherwise from the `mass` column, and otherwise from the atom types.
`read` returns the last frame of a file with several frames.
`read(..., index=0)` returns the first, and `read(..., index=":")` returns a list of all frames.

### Cell and periodic boundaries

pyscal uses `atoms.cell` and `atoms.pbc`:

- A plain XYZ file has no cell, and ASE sets `pbc=False`. pyscal then treats the structure as isolated, without periodic images. If the file is a periodic crystal, set the cell before the analysis:

  ```python
  atoms.set_cell([[a, 0, 0], [0, b, 0], [0, 0, c]])
  atoms.set_pbc(True)
  ```

- Mixed periodicity, for example a slab with `pbc=[True, True, False]`, is supported. There are no periodic images along the directions with `pbc=False`.
- Triclinic cells are supported.
- `radial_distribution_function` and `entropy` need the density, and therefore a cell that is periodic in all three directions.

## Trajectories

A trajectory with many frames can be read with `read(..., index=":")`, but this keeps every frame in memory.
`pyscal.Trajectory` reads a LAMMPS dump file lazily.
It finds where each frame starts and reads a frame only when it is used.

```{code-cell} ipython3
import pyscal

trajectory = pyscal.Trajectory("traj.light")
trajectory
```

Indexing a `Trajectory` gives a `Timeslice`, and `to_atoms()` converts it to a list of ASE `Atoms` objects, one per frame:

```{code-cell} ipython3
first = trajectory[0].to_atoms()[0]
last_three = trajectory[-3:].to_atoms()
len(first), len(last_three)
```

| Operation | Result |
|---|---|
| `trajectory[i]`, `trajectory[i:j]` | a `Timeslice` with one or several frames |
| `timeslice.to_atoms(species=["Cu", "Zr"])` | a list of `Atoms`, with LAMMPS types 1, 2 mapped to Cu, Zr |
| `timeslice.to_atoms(customkeys=["c_pe"])` | also reads the per-atom column `c_pe` into `atoms.arrays` |
| `timeslice_a + timeslice_b` | a `Timeslice` with the frames of both |
| `timeslice.to_file("frames.dump")` | writes the frames to a LAMMPS dump file |

A typical analysis loops over the frames and keeps only the results:

```{code-cell} ipython3
import numpy as np

q6_mean = []
for i in range(trajectory.nblocks):
    atoms = trajectory[i].to_atoms()[0]
    pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)
    q6 = pyscal.steinhardt_parameter(atoms, l=6, averaged=True)[0]
    q6_mean.append(q6.mean())

np.round(q6_mean, 3)
```

The file holds 10 frames of a liquid, and the mean $\bar{q}_6$ stays at the value of a liquid in every frame.
In a solidification run, the same loop would show $\bar{q}_6$ rising as crystals grow.

## Writing results

The per-atom results are in `atoms.arrays`, so they can be written with the positions.
The extended XYZ format keeps them as named columns that [OVITO](https://www.ovito.org) and other programs read:

```{code-cell} ipython3
from ase.io import write

write("analysed.extxyz", atoms,
      columns=["symbols", "positions", "pyscal_avg_q6"], write_info=False)
```

```{code-cell} ipython3
:tags: [remove-cell]
os.remove("analysed.extxyz")
```

`write_info=False` leaves out the neighbor data in `atoms.info`.
With [`store_rows=False`](neighbors.md#how-the-neighbors-are-stored) in `find_neighbors`, `atoms.write("analysed.extxyz")` also works without selecting columns.

To keep the results of every frame of a trajectory, collect the arrays in a list, or write one file per frame.
