/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#include "triplet_list.hh"

#include <cmath>

#include "error.hh"

namespace matscipy {

error_t triplet_list(index_t n_first, const index_t *first_i,
                     index_t n_absdist, const real_t *absdist, real_t cutoff,
                     std::vector<index_t> &ij_t, std::vector<index_t> &ik_t) {
    clear_error();
    ij_t.clear();
    ik_t.clear();

    if (n_first < 0 || (n_first > 0 && !first_i)) {
        return set_invalid_argument("triplet_list: invalid row-start array.");
    }
    if (absdist && (n_absdist < 0 || std::isnan(cutoff))) {
        return set_invalid_argument(
            "triplet_list: invalid distance array or cutoff.");
    }
    /* first_neighbours() emits -1 for atoms before the first pair and then
       0 for the first atom with one; -1 is only meaningful in that form. */
    for (index_t r = 0; r < n_first; r++) {
        if (r > 0 && first_i[r - 1] == -1 && first_i[r] > 0) {
            return set_invalid_argumentf(
                "triplet_list: leading -1 row starts must be followed by 0, "
                "got %lld at position %lld.",
                static_cast<long long>(first_i[r]), static_cast<long long>(r));
        }
        if (first_i[r] < -1) {
            return set_invalid_argumentf(
                "triplet_list: row start %lld at position %lld is negative.",
                static_cast<long long>(first_i[r]), static_cast<long long>(r));
        }
        if (r > 0 && first_i[r] < first_i[r - 1]) {
            return set_invalid_argumentf(
                "triplet_list: row starts must be non-decreasing (decrease at "
                "position %lld).",
                static_cast<long long>(r));
        }
        if (absdist && first_i[r] > n_absdist) {
            return set_invalid_argumentf(
                "triplet_list: row start %lld at position %lld exceeds the "
                "number of distances (%lld).",
                static_cast<long long>(first_i[r]), static_cast<long long>(r),
                static_cast<long long>(n_absdist));
        }
    }

    for (index_t r = 0; r < n_first - 1; r++) {
        /* A -1 marks "no entries"; never let it index absdist. */
        const index_t begin = first_i[r] < 0 ? 0 : first_i[r];
        const index_t end = first_i[r + 1] < 0 ? 0 : first_i[r + 1];
        for (index_t ij = begin; ij < end; ij++) {
            for (index_t ik = begin; ik < end; ik++) {
                if (ij == ik) continue;
                if (absdist &&
                    (absdist[ij] >= cutoff || absdist[ik] >= cutoff)) {
                    continue;
                }
                ij_t.push_back(ij);
                ik_t.push_back(ik);
            }
        }
    }
    return NL_SUCCESS;
}

}  // namespace matscipy
