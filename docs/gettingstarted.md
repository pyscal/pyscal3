 # Installation

`pyscal3` can be installed on Linux, macOS, and Windows. The following instructions will help install `pyscal3`:

````{tab-set}
```{tab-item} pip
`pip install pyscal3`
```

```{tab-item} conda
`conda install -c conda-forge pyscal3`
```

```{tab-item} from source
We strongly recommend creating a conda environment for the installation. To see how you can install conda see [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

Once a conda distribution is available, the following steps will help set up an environment to use `pyscal3`. First step is to clone the repository.

`git clone https://github.com/pyscal/pyscal3.git`

`cd pyscal3`  
`pip install .`
```
````

# Reading structures

pyscal works on [ASE](https://wiki.fysik.dtu.dk/ase/) `Atoms` objects, so any file format that ASE can read is available without conversion:

``` python
import pyscal
from ase.io import read

atoms = read('dump.lammpstrj', format='lammps-dump-text')   # LAMMPS dump
atoms = read('POSCAR')                                        # VASP
atoms = read('structure.cif')                                 # CIF
atoms = read('cluster.xyz')                                   # plain xyz

pyscal.find_neighbors(atoms, method='cutoff', cutoff=0)
q6 = pyscal.steinhardt_parameter(atoms, l=6)[0]
```

Plain `.xyz` files carry no cell, so ASE returns a zero cell with `pbc=False`. pyscal treats such structures as isolated (non-periodic) systems: neighbors are found from the coordinates alone and no periodic images are created. If the file describes a periodic crystal, set the cell and periodicity before the analysis, or use the extended xyz format whose `Lattice="..."` header ASE reads automatically:

``` python
atoms = read('crystal.xyz')
atoms.set_cell([[a, 0, 0], [0, b, 0], [0, 0, c]])
atoms.set_pbc(True)
```

Mixed periodicity (for example a slab with `pbc=[True, True, False]`) is respected as well. Descriptors that need the global density, such as the radial distribution function and the entropy parameter, require a fully periodic cell.

Multi-frame LAMMPS dump files can be read frame by frame with ASE (`read(..., index=':')`) or, for large files, with the lazy `pyscal.Trajectory` reader described in the [trajectory example](../examples/13_trajectory_module).
