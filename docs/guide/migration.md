# Migrating from pyscal 3

pyscal 4 has a new interface.
Scripts written for pyscal 3 need small changes to run.
This page lists the differences and shows how to translate a typical script.

## What changed

- **ASE `Atoms` replaces `System`.** pyscal no longer has its own `System` and `Atoms` classes. Every function takes an ASE `Atoms` object as its first argument. Structures are read with `ase.io.read` or built with `ase.build` or `pyscal.structures`.
- **Functions replace methods.** `sys.find.neighbors(...)` becomes `pyscal.find_neighbors(atoms, ...)`, and `sys.calculate.steinhardt_parameter(...)` becomes `pyscal.steinhardt_parameter(atoms, ...)`.
- **Results are stored on the `Atoms` object.** Per-atom results are in `atoms.arrays["pyscal_<name>"]`, other results in `atoms.info["pyscal_<name>"]`. The functions also return them.
- **New descriptors.** Wigner $W_l$, Minkowski structure metrics, Ackland–Jones classification, coordination numbers, angular and bond length distributions, deformation descriptors, Wigner–Seitz defect analysis and ACE descriptors are new in pyscal 4.

## A script in both versions

pyscal 3:

```python
from pyscal3 import System

sys = System("conf.dump")
sys.find.neighbors(method="cutoff", cutoff=3)
sys.calculate.steinhardt_parameter([4, 6])
sys.find.solids()
solid = sys.atoms.solid
```

pyscal 4:

```python
import pyscal
from ase.io import read

atoms = read("conf.dump", format="lammps-dump-text")
pyscal.find_neighbors(atoms, method="cutoff", cutoff=3)
pyscal.steinhardt_parameter(atoms, l=[4, 6])
pyscal.find_solids(atoms)
solid = atoms.arrays["pyscal_solid"]
```

## Results that changed

Some results differ from pyscal 3 because errors were corrected.
The [changelog](https://github.com/pyscal/pyscal3/blob/main/CHANGELOG.md) lists them.
The most visible ones:

- `radial_distribution_function` returns a normalised $g(r)$, which tends to 1 at large distances.
- `short_range_order` computes the Warren–Cowley parameter $\alpha_{AB} = 1 - p_{AB} / c_B$.
- SANN neighbor lists include the correct number of atoms.
- `W_l` values do not depend on the orientation of the crystal.

To keep using pyscal 3, install `pyscal3<4`.
