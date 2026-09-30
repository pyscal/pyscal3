# Vendored matscipy-neighbours

This directory is a copy of the CPU part of the C++ core of
[matscipy-neighbours](https://github.com/libAtoms/matscipy-neighbours),
which pyscal3 uses for its neighbour search. It is compiled into
`pyscal3.csystem` by `setup.py`.

- Upstream: https://github.com/libAtoms/matscipy-neighbours
- Commit: `1d38a67` (2026-09-28)
- Licence: MIT, see `LICENSE.md`. Authors are listed in `AUTHORS`.
- Local changes: none.

## Files

The sources of the upstream `neighbours` library target without the GPU backend,
from `src/libneighbours/`, and the headers they include:

    cell_list.cc  cell_list.hh
    error.cc  error.hh
    first_neighbours.cc  first_neighbours.hh
    memory_space.cc  memory_space.hh
    neighbour_list.cc  neighbour_list.hh  neighbour_visit.hh
    tools.cc  tools.hh
    triplet_list.cc  triplet_list.hh
    types.hh

The GPU sources (`memory_space_gpu.cc`, `device_primitives.cc`,
`neighbour_list_gpu.cc`, `device.hh`) are not copied. The GPU code in the
copied headers is behind `MATSCIPY_ENABLE_CUDA` and `MATSCIPY_ENABLE_HIP`,
which pyscal3 does not define.

pyscal3 builds these files without OpenMP, so the `#pragma omp` lines are
ignored and the search runs on one thread.

## Updating

1. Clone upstream and check out the new commit.
2. Copy the files listed above from `src/libneighbours/`, and `LICENSE.md`
   and `AUTHORS` from the repository root. If upstream adds a source file to
   the CPU `neighbours` target in `src/libneighbours/CMakeLists.txt`, add it
   here and to `setup.py`.
3. Update the commit and date above.
4. Rebuild and run the test suite, in particular `tests/test_neighbor_backend.py`,
   which compares the search with stored reference results.
