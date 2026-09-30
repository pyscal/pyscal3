/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#include "tools.hh"

#include <cmath>

namespace matscipy {

void cross_product(const real_t *a, const real_t *b, real_t *c) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

real_t norm(const real_t *a) {
    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

/* mat_mul_vec, bin_wrap, bin_trunc and position_to_cell_index are host+device
   header-inline functions in tools.hh, shared by the CPU and GPU paths. */

}  // namespace matscipy
