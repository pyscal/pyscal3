# Contributing

Bug reports, questions, fixes and new descriptors are welcome.

## Reporting a bug or asking for a feature

Open an issue on the [issue page](https://github.com/pyscal/pyscal3/issues) of the repository.
For a bug, include a short script that shows it, the output you expected, and the versions of pyscal, ASE and Python (`pyscal.__version__`).
For other questions, write to [sarath.menon@pyscal.org](mailto:sarath.menon@pyscal.org).

## Development setup

1. Fork the repository on GitHub and clone your fork.
2. Create a separate environment, for example with conda, and install pyscal from the clone. The build needs a C++17 compiler.

   ```bash
   cd pyscal3
   pip install -e .
   ```

3. Create a branch for your change: `git checkout -b my-change`.
4. Run the tests with `pytest tests`. After changes to the C++ code, rebuild with `pip install -e .` before testing.

## Adding a descriptor

- Write the descriptor as a function in `src/pyscal3/descriptors.py` that takes an ASE `Atoms` object as its first argument, returns its result and stores it on the `Atoms` object under a key that starts with `pyscal_`. Add it to `__all__` in `src/pyscal3/__init__.py`.
- Descriptors that use neighbors read them with `neighbor_arrays` from `src/pyscal3/_bridge.py`, which gives the flat neighbor arrays described in [Finding neighbors](guide/neighbors.md#how-the-neighbors-are-stored).
- Per-atom loops that need speed go into C++ in `src/pyscal3/`, with pybind11 bindings in `system_binding.cpp`. Use `pyscal::parallel_for` from `parallel.h` for loops over atoms, and write each atom's result only from the thread that handles it, so that the result does not depend on the number of threads.
- Write a docstring in the [numpydoc format](https://numpydoc.readthedocs.io/en/latest/format.html), with the reference to the original publication.
- Add tests in `tests/`. Compare with published values or an independent implementation where possible.

## Documentation

The documentation is a [Jupyter Book](https://jupyterbook.org).
The pages in `docs/` are Markdown files, and most of them are notebooks in the [MyST format](https://myst-nb.readthedocs.io): their code cells run when the book is built, so the results and figures always match the code.

```bash
pip install -r requirements.txt
jupyter-book build .
```

The book is written to `_build/html`.
A single page can be built with `jupyter-book build docs/descriptors/cna.md`.

A new descriptor page goes into `docs/descriptors/`, and is added to `_toc.yml` and to the overview in `docs/descriptors/index.md`.
Use the page of a similar descriptor as a template.
Each page has a definition, a short usage example with the stored keys, one or two examples on realistic structures with figures, the pitfalls, and the references.
Figures use the shared style in `docs/_plotstyle.py`, and their plotting code goes into cells tagged `hide-input`.

## Pull requests

Open a pull request against the `main` branch.
Before you open it, make sure the tests pass and the new feature has tests, a docstring and documentation.
The tests run automatically on Linux, macOS and Windows.
