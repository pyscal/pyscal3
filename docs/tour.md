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

# A first analysis

This page goes through a complete analysis in five steps:
read a structure, find the neighbors of each atom, compute two descriptors, look at the results, and save them.
Every code cell on this page runs as shown.

```{code-cell} ipython3
:tags: [remove-cell]
import os, sys, warnings
sys.path.insert(0, os.path.abspath("."))
from _plotstyle import EXAMPLES, COLOURS, DARK, figure, label, note
os.chdir(EXAMPLES)
warnings.simplefilter("ignore")
%config InlineBackend.figure_format = "retina"
```

## 1. Read a structure

pyscal works on ASE [`Atoms`](https://wiki.fysik.dtu.dk/ase/ase/atoms.html) objects.
Anything that ASE can read or build can be analysed, without conversion.
Here we read a snapshot from a molecular dynamics (MD) simulation of an fcc crystal at finite temperature.
The file is in the `examples` folder of the pyscal repository.

```{code-cell} ipython3
import pyscal
from ase.io import read

atoms = read("conf.fcc.dump", format="lammps-dump-text")
atoms
```

## 2. Find the neighbors

Most descriptors describe an atom through its neighbors, so the first step is almost always `find_neighbors`.

```{code-cell} ipython3
pyscal.find_neighbors(atoms, method="cutoff", cutoff=0)
```

`cutoff=0` selects the adaptive cutoff.
Each atom gets its own cutoff, computed from the distances to its six nearest atoms, so no distance has to be chosen by hand.
The neighbor list is stored on `atoms` and every descriptor below reuses it.
[Finding neighbors](guide/neighbors) compares the available methods.

## 3. Compute descriptors

The Steinhardt parameter $q_6$ measures the orientational order of the neighbors of an atom.
It is close to 0.57 for an atom in a perfect fcc crystal and lower in a liquid.

```{code-cell} ipython3
q6 = pyscal.steinhardt_parameter(atoms, l=6)[0]
q6.mean()
```

Common neighbor analysis (CNA) assigns a crystal structure to each atom.
It finds its own neighbors and does not change the neighbor list stored in step 2.

```{code-cell} ipython3
counts = pyscal.common_neighbor_analysis(atoms)
counts
```

```{code-cell} ipython3
:tags: [remove-cell]
from myst_nb import glue
glue("fcc_percent", round(100 * counts["fcc"] / len(atoms)), display=False)
```

Only {glue}`fcc_percent` % of the atoms are labelled fcc, and most of the others are labelled `others`.
The thermal displacements are large, and many atoms do not have the exact CNA signature of fcc.
[Common neighbor analysis](descriptors/cna) shows how the labels depend on temperature.

## 4. Look at the results

Every function returns its result and also stores it on `atoms`, with keys that start with `pyscal_`.
Per-atom values are in `atoms.arrays`, which ASE keeps aligned with the atoms when it reorders, deletes or repeats them.
Values that are not one number per atom are in `atoms.info`.

```{code-cell} ipython3
sorted(key for key in atoms.arrays if key.startswith("pyscal_"))
```

```{code-cell} ipython3
atoms.arrays["pyscal_structure"][:10]
```

`pyscal_structure` holds the CNA label of each atom: 0 others, 1 fcc, 2 hcp, 3 bcc, 4 icosahedral.

## 5. Compare with a liquid

A descriptor becomes useful when it separates the structures of interest.
We repeat the steps for a liquid snapshot and compare the distributions of $q_6$.
We also compute the averaged parameter $\bar{q}_6$, which averages $q_6$ over each atom and its neighbors.

```{code-cell} ipython3
liquid = read("conf.lqd.Al.dump", format="lammps-dump-text")

results = {}
for name, structure in [("fcc", atoms), ("liquid", liquid)]:
    pyscal.find_neighbors(structure, method="cutoff", cutoff=0)
    q6 = pyscal.steinhardt_parameter(structure, l=6)[0]
    q6_averaged = pyscal.steinhardt_parameter(structure, l=6, averaged=True)[0]
    results[name] = (q6, q6_averaged)
```

```{code-cell} ipython3
:tags: [hide-input]
import numpy as np

mp = figure(columns=2, ratio=0.42, wspace=0.25)
bins = np.linspace(0, 0.7, 50)
for k, title in enumerate(["(a)  $q_6$", r"(b)  $\bar{q}_6$"]):
    ax = mp[0, k]
    for name, values in results.items():
        ax.hist(values[k], bins=bins, density=True, color=COLOURS[name], alpha=0.75,
                ec=DARK, lw=0.5, label=name)
    label(ax, title)
    ax.set_xlabel(title.split()[-1])
    ax.set_xlim(0, 0.7)
mp[0, 0].set_ylabel("Probability density")
mp[0, 0].legend(frameon=False);
```

The distributions of $q_6$ (a) overlap, so a threshold on $q_6$ would misclassify many atoms.
The distributions of $\bar{q}_6$ (b) are well separated, because averaging over the neighbors removes much of the thermal noise.
[Steinhardt parameters](descriptors/steinhardt) explains the averaging.

## 6. Save the results

To look at the results in a visualisation program such as [OVITO](https://www.ovito.org), write the positions and the per-atom results to an extended XYZ file.

```{code-cell} ipython3
from ase.io import write

write("fcc_analysed.extxyz", atoms,
      columns=["symbols", "positions", "pyscal_q6", "pyscal_structure"],
      write_info=False)
```

`columns` selects the per-atom arrays to write, and `write_info=False` leaves out the neighbor data in `atoms.info`, which the extended XYZ format cannot store when atoms have different numbers of neighbors.
In OVITO, the two columns appear as particle properties that can be used for colouring and selection.

```{code-cell} ipython3
:tags: [remove-cell]
os.remove("fcc_analysed.extxyz")
```

## Next steps

- [Finding neighbors](guide/neighbors): choose the neighbor method that fits your system.
- [Descriptors](descriptors/index): what each descriptor measures and when to use it.
- [API reference](api): all functions and their parameters.
