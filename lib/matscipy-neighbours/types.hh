/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_TYPES_HH
#define MATSCIPY_TYPES_HH

#include <cstdint>

/* Mark a small leaf function callable from both host and device so the CPU and
   GPU paths share one definition. Expands to nothing for the host compiler and
   to `__host__ __device__` when nvcc/hipcc compiles the translation unit. */
#if defined(__CUDACC__) || defined(__HIPCC__)
#define MATSCIPY_HD __host__ __device__
#else
#define MATSCIPY_HD
#endif

namespace matscipy {

/* Integer and floating-point types used throughout the core. Plain aliases so
   the core never needs to include the NumPy headers. `index_t` is 64-bit so
   atom, cell and pair counts cannot overflow; the Python layer exposes it as
   NumPy int64 / DLPack int64. */
using index_t = std::int64_t;
using real_t = double;

/* Error code returned by core routines. */
using error_t = int;
constexpr error_t NL_SUCCESS = 0;
constexpr error_t NL_ERROR = -1;            /* internal / runtime failure */
constexpr error_t NL_INVALID_ARGUMENT = -2; /* unusable caller input */
constexpr error_t NL_OUT_OF_MEMORY = -3;    /* an allocation failed */

/* Bit flags selecting which per-pair quantities a neighbour-list call computes.
   The Python layer maps the "ijdDS" quantity string onto these. */
enum Quantity : int {
    QUANTITY_FIRST = 1 << 0,    /* i: first atom index */
    QUANTITY_SECOND = 1 << 1,   /* j: second atom index */
    QUANTITY_DISTVEC = 1 << 2,  /* D: distance vector */
    QUANTITY_ABSDIST = 1 << 3,  /* d: absolute distance */
    QUANTITY_SHIFT = 1 << 4,    /* S: cell shift */
};

}  // namespace matscipy

#endif
