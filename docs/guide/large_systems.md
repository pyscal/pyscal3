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

# Large systems

pyscal analyses systems of millions of atoms on a laptop.
For one million fcc atoms with 12 neighbors each, the neighbor search, common neighbor analysis and $q_6$ each take a fraction of a second on a 14 core laptop.
This page describes what determines the time and memory, and the options that matter for large systems.

```{code-cell} ipython3
:tags: [remove-cell]
import os, sys, warnings
sys.path.insert(0, os.path.abspath(".."))
from _plotstyle import EXAMPLES, COLOURS, DARK, figure, label, note
os.chdir(EXAMPLES)
warnings.simplefilter("ignore")
%config InlineBackend.figure_format = "retina"
```

## Threads

The neighbor search, common neighbor analysis and the per-atom loops of the descriptors run on all CPU cores that are available to the process.
The results are the same for any number of threads, to the last bit.

```{code-cell} ipython3
import pyscal

pyscal.get_num_threads()
```

`pyscal.set_num_threads(n)` changes the number of threads, and `pyscal.set_num_threads()` restores the default.
The environment variable `PYSCAL_NUM_THREADS`, or else `OMP_NUM_THREADS`, sets the default when pyscal is imported.

When several analyses run in parallel processes, for example with `multiprocessing` or on a cluster with one job per core, set the number of threads to 1 in each process.
Otherwise every process starts one thread per core, and the threads compete for the cores.

## Time and memory

The cost of most calculations is proportional to the number of neighbor pairs, which is the number of atoms times the number of neighbors per atom.
The number of neighbors grows with the cube of the cutoff: in fcc, a cutoff of 3 Å around a copper atom contains 12 neighbors and a cutoff of 5 Å contains 42.
Use the smallest cutoff that contains the shells the descriptor needs.

`find_neighbors` stores, for every neighbor pair, its index, distance, weight, angles and vector: 64 bytes per pair in the flat arrays.
For one million atoms with 12 neighbors this is about 0.8 GB.

## Skipping the per-atom rows

`find_neighbors` also stores the neighbors as rows, one per atom (see [How the neighbors are stored](neighbors.md#how-the-neighbors-are-stored)).
When all atoms have the same number of neighbors, the rows are views of the flat arrays and cost almost nothing.
When the numbers differ, the rows are Python lists, which take much more time and memory than the flat arrays.
All descriptors read the flat arrays, so the rows can be skipped with `store_rows=False`.

The cell below times both options for a crystal with random displacements, where atoms have different numbers of neighbors within the cutoff.

```{code-cell} ipython3
import time
import pandas as pd
from ase.build import bulk

crystal = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(20)
crystal.rattle(0.1, seed=1)

timings = {}
for store_rows in (True, False):
    atoms = crystal.copy()
    start = time.perf_counter()
    pyscal.find_neighbors(atoms, method="cutoff", cutoff=5.0, store_rows=store_rows)
    timings[f"store_rows={store_rows}"] = time.perf_counter() - start

pd.Series(timings, name="time (s)").round(3)
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
glue("rows_ratio", round(timings["store_rows=True"] / timings["store_rows=False"]), display=False)
glue("natoms_rows", len(crystal), display=False)
```

For these {glue}`natoms_rows` atoms, skipping the rows makes the search about {glue}`rows_ratio` times faster.
The times are for the computer that built this page.
`store_rows=False` also makes it possible to write `atoms` to an extended XYZ file directly.

## Recommendations for large systems

- Use `store_rows=False`, unless a script reads `atoms.info["pyscal_neighbors"]` and the other row keys.
- Choose the smallest cutoff that the descriptor needs.
- Read long trajectories frame by frame with [`pyscal.Trajectory`](files.md#trajectories).
- Keep only the results you need from each frame, for example `atoms.arrays["pyscal_structure"].copy()`, and let the `Atoms` object with its neighbor data be freed.
