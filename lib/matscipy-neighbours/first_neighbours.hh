/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_FIRST_NEIGHBOURS_HH
#define MATSCIPY_FIRST_NEIGHBOURS_HH

#include <vector>

#include "types.hh"

namespace matscipy {

/*
 * Build the row-start ("seed") array of a neighbour list whose first index
 * array i_n is sorted and non-decreasing.
 *
 * n      number of rows (atoms); seed must have room for n + 1 entries
 * nn     length of i_n
 * i_n    [nn] sorted first-atom indices, each in [0, n)
 * seed   [n + 1] output: seed[k] is the position in i_n where atom k starts,
 *        seed[n] == nn. Atoms before the first entry keep seed value -1.
 *
 * Handles nn == 0 (empty list) without dereferencing i_n. Returns NL_SUCCESS,
 * or NL_INVALID_ARGUMENT (with a message via set_error) for n < 0, an index
 * outside [0, n) or an unsorted i_n; seed is untouched in that case.
 */
error_t first_neighbours(index_t n, index_t nn, const index_t *i_n,
                         index_t *seed);

/*
 * Given a sorted, contiguous array starting at 0, fill `seed` with the array
 * pointing to the index jumps (the row-start array for the implied number of
 * distinct values). `seed` gets length (number of distinct values) + 1.
 * Returns NL_INVALID_ARGUMENT if the array does not start at 0 or has a gap
 * or a decrease between consecutive entries.
 */
error_t get_jump_indicies(index_t nn, const index_t *sorted,
                          std::vector<index_t> &seed);

}  // namespace matscipy

#endif
