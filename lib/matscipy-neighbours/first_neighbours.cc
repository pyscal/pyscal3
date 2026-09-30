/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#include "first_neighbours.hh"

#include "error.hh"

namespace matscipy {

error_t first_neighbours(index_t n, index_t nn, const index_t *i_n,
                         index_t *seed) {
    clear_error();

    if (n < 0) {
        return set_invalid_argument(
            "first_neighbours: number of atoms must be non-negative.");
    }
    if (nn < 0 || (nn > 0 && !i_n)) {
        return set_invalid_argument(
            "first_neighbours: invalid neighbour index array.");
    }
    /* Validate before writing anything: every index must address a row and
       the array must be sorted (the seeds are computed from the jumps). */
    for (index_t k = 0; k < nn; k++) {
        if (i_n[k] < 0 || i_n[k] >= n) {
            return set_invalid_argumentf(
                "first_neighbours: index %lld at position %lld is outside "
                "[0, %lld).",
                static_cast<long long>(i_n[k]), static_cast<long long>(k),
                static_cast<long long>(n));
        }
        if (k > 0 && i_n[k] < i_n[k - 1]) {
            return set_invalid_argumentf(
                "first_neighbours: index array must be sorted (decrease at "
                "position %lld).",
                static_cast<long long>(k));
        }
    }

    for (index_t k = 0; k <= n; k++) {
        seed[k] = -1;
    }

    /* Empty neighbour list: every row starts (and ends) at 0. */
    if (nn == 0) {
        for (index_t k = 0; k <= n; k++) seed[k] = 0;
        return NL_SUCCESS;
    }

    seed[i_n[0]] = 0;

    for (index_t k = 1; k < nn; k++) {
        if (i_n[k] != i_n[k - 1]) {
            for (index_t l = i_n[k - 1] + 1; l <= i_n[k]; l++) {
                seed[l] = k;
            }
        }
    }

    for (index_t k = i_n[nn - 1] + 1; k <= n; k++) {
        seed[k] = nn;
    }
    return NL_SUCCESS;
}

error_t get_jump_indicies(index_t nn, const index_t *sorted,
                          std::vector<index_t> &seed) {
    clear_error();
    seed.clear();

    if (nn < 0 || (nn > 0 && !sorted)) {
        return set_invalid_argument("get_jump_indicies: invalid array.");
    }
    if (nn > 0 && sorted[0] != 0) {
        return set_invalid_argument(
            "get_jump_indicies: array must start at 0.");
    }
    /* Number of distinct values = number of jumps + 1; each jump must be
       exactly +1 so the values are contiguous. */
    index_t n = 0;
    for (index_t i = 0; i < nn - 1; i++) {
        const index_t step = sorted[i + 1] - sorted[i];
        if (step < 0 || step > 1) {
            return set_invalid_argumentf(
                "get_jump_indicies: array must be sorted and contiguous "
                "(step %lld at position %lld).",
                static_cast<long long>(step), static_cast<long long>(i));
        }
        if (step == 1) n++;
    }
    n++;

    seed.assign(n + 1, 0);
    return first_neighbours(n, nn, sorted, seed.data());
}

}  // namespace matscipy
