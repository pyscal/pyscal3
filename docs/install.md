# Installation

pyscal runs on Linux, macOS and Windows with Python 3.10 or later.
It is distributed as **`pyscal3`**, because the name `pyscal` on PyPI belongs to an unrelated package.
After installation, `import pyscal` and `import pyscal3` both work and give the same module.

````{tab-set}
```{tab-item} pip
    pip install pyscal3

Wheels with the compiled C++ code are available for the common platforms, so no compiler is needed.
```

```{tab-item} conda
    conda install -c conda-forge pyscal3
```

```{tab-item} from source
Building from source needs a C++17 compiler.

    git clone https://github.com/pyscal/pyscal3.git
    cd pyscal3
    pip install .

To run the tests, install `pytest` and run `pytest tests` in the repository.
```
````

pyscal depends on `numpy`, `scipy`, `ase` and `pyyaml`, which are installed with it.

## Checking the installation

```python
import pyscal
from ase.build import bulk

atoms = bulk("Cu", "fcc", a=3.61, cubic=True).repeat(4)
pyscal.common_neighbor_analysis(atoms)
```

This prints `{'others': 0, 'fcc': 256, 'hcp': 0, 'bcc': 0, 'ico': 0}`.

## Number of threads

pyscal uses all CPU cores that are available to the process.
The results do not depend on the number of threads.
To use fewer threads, for example when several analyses run side by side, call

```python
pyscal.set_num_threads(4)
```

or set the environment variable `PYSCAL_NUM_THREADS` (or `OMP_NUM_THREADS`) before Python starts.
`pyscal.get_num_threads()` returns the current number.

## Coming from pyscal 3

pyscal 4 has a new interface, built around ASE `Atoms` and functions instead of the `System` class.
[Migrating from pyscal 3](guide/migration) shows how to translate existing scripts.
To keep using the old interface, install `pyscal3<4`.
