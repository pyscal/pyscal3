/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_TRIPLET_LIST_HH
#define MATSCIPY_TRIPLET_LIST_HH

#include <vector>

#include "types.hh"

namespace matscipy {

/*
 * Build a triplet list from a first-neighbour (row-start) array. For each row
 * r, every ordered pair (ij, ik) of distinct neighbour-list entries in
 * [first_i[r], first_i[r + 1]) becomes a triplet.
 *
 * n_first    length of first_i
 * first_i    [n_first] row-start array, non-decreasing. A leading -1 (what
 *            first_neighbours() emits for atoms before the first entry) is
 *            accepted and treated as "no entries".
 * n_absdist  length of absdist (ignored when absdist is nullptr)
 * absdist    [n_absdist] absolute distances indexed by neighbour-list entry,
 *            or nullptr
 * cutoff     if absdist is given, skip triplets where either leg >= cutoff
 * ij_t       output: first leg of each triplet
 * ik_t       output: second leg of each triplet
 *
 * Returns NL_SUCCESS, or NL_INVALID_ARGUMENT (message via set_error) if the
 * row starts decrease, or exceed n_absdist when distances are given.
 */
error_t triplet_list(index_t n_first, const index_t *first_i,
                     index_t n_absdist, const real_t *absdist, real_t cutoff,
                     std::vector<index_t> &ij_t, std::vector<index_t> &ik_t);

}  // namespace matscipy

#endif
