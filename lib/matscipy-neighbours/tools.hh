/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_TOOLS_HH
#define MATSCIPY_TOOLS_HH

#include <cmath>

#include "types.hh"

namespace matscipy {

/* c = a x b for 3-vectors. */
void cross_product(const real_t *a, const real_t *b, real_t *c);

/* Euclidean norm (length) of a 3-vector. */
real_t norm(const real_t *a);

/* The cell-indexing leaf helpers below are host+device header-inline functions,
   shared by the CPU loops and the GPU kernels from a single definition each. */

/* vout = mat . vin, with mat a row-major 3x3 matrix. */
MATSCIPY_HD inline void mat_mul_vec(const real_t *mat, const real_t *vin,
                                    real_t *vout) {
    for (int i = 0; i < 3; i++) {
        vout[i] = 0.0;
        for (int j = 0; j < 3; j++) vout[i] += mat[3 * i + j] * vin[j];
    }
}

/* Map cell index i back into [0, n) by shifting by multiples of n (periodic).
   O(1) for any distance from the cell. */
MATSCIPY_HD inline index_t bin_wrap(index_t i, index_t n) {
    index_t r = i % n;
    return r < 0 ? r + n : r;
}

/* Clamp cell index i into [0, n) (non-periodic). */
MATSCIPY_HD inline index_t bin_trunc(index_t i, index_t n) {
    if (i < 0)
        i = 0;
    else if (i >= n)
        i = n - 1;
    return i;
}

/* True for a finite value. Written without <cmath> so it is usable in device
   code under both nvcc and hipcc: x - x is 0 for finite x, NaN for +-inf and
   NaN, and IEEE semantics (no fast-math) forbid folding it to 0. */
MATSCIPY_HD inline bool is_finite(real_t x) { return x - x == 0.0; }

/* Largest magnitude of a raw (unwrapped) cell coordinate. Scaled fractional
   coordinates beyond this (atoms astronomically far outside the cell, or
   non-finite positions the caller did not reject) are clamped before the
   float->integer cast, which would otherwise be undefined behaviour. */
constexpr real_t MAX_RAW_CELL_COORD = 1099511627776.0; /* 2^40 */

MATSCIPY_HD inline index_t scaled_coord_to_cell(real_t x) {
    /* NaN compares false to everything, so it falls through to the cast; map
       it to the upper clamp explicitly. */
    if (!(x > -MAX_RAW_CELL_COORD)) {
        return x != x ? static_cast<index_t>(MAX_RAW_CELL_COORD)
                      : -static_cast<index_t>(MAX_RAW_CELL_COORD);
    }
    if (!(x < MAX_RAW_CELL_COORD)) return static_cast<index_t>(MAX_RAW_CELL_COORD);
    return static_cast<index_t>(std::floor(x));
}

/* Map a Cartesian position to (unwrapped) integer cell indices. */
MATSCIPY_HD inline void position_to_cell_index(const real_t *cell_origin,
                                               const real_t *inv_cell,
                                               const real_t *ri, index_t n1,
                                               index_t n2, index_t n3,
                                               index_t *c1, index_t *c2,
                                               index_t *c3) {
    real_t dri[3], si[3];
    for (int i = 0; i < 3; i++) dri[i] = ri[i] - cell_origin[i];
    mat_mul_vec(inv_cell, dri, si);
    *c1 = scaled_coord_to_cell(si[0] * n1);
    *c2 = scaled_coord_to_cell(si[1] * n2);
    *c3 = scaled_coord_to_cell(si[2] * n3);
}

}  // namespace matscipy

#endif
