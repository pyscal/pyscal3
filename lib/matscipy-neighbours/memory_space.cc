/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 *
 * Host backend for the memory-space abstraction. Device backends (CUDA/HIP)
 * live in their own translation units, compiled only when enabled.
 */

#include <cstdlib>
#include <new>

#include "memory_space.hh"

namespace matscipy {
namespace detail {

/* Out of memory throws std::bad_alloc, like the device hooks. */
void *alloc_host(std::size_t bytes) {
    void *ptr = std::malloc(bytes);
    if (!ptr) throw std::bad_alloc();
    return ptr;
}

void free_host(void *ptr) { std::free(ptr); }

}  // namespace detail
}  // namespace matscipy
