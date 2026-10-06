/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_NEIGHBOUR_LIST_HH
#define MATSCIPY_NEIGHBOUR_LIST_HH

#include <cstddef>
#include <vector>

#include "cell_list.hh"
#include "types.hh"

namespace matscipy {

/* Minimal non-owning view over a contiguous buffer (a pointer + length). Works
   for host buffers and device pointers alike; on a device pointer only size()
   / data() are meaningful on the host (indexing happens in kernels). */
template <typename T>
struct Span {
    T *ptr = nullptr;
    std::size_t len = 0;
    T *data() const { return ptr; }
    std::size_t size() const { return len; }
    bool empty() const { return len == 0; }
    T &operator[](std::size_t i) const { return ptr[i]; }
    T *begin() const { return ptr; }
    T *end() const { return ptr + len; }
};

/* A 3-vector value (distance vector or shift) read out of a flattened buffer. */
template <typename T>
struct Vec3 {
    T x, y, z;
};

/* Result of a neighbour-list computation. Only the buffers whose quantity flag
   was requested are filled; the rest stay empty. `npairs` is the number of
   neighbour pairs found, sorted by the first index i. Distance vectors and
   shifts are stored row-major with three entries per pair.

   The *_view() accessors give typed Spans aliasing the buffers instead of raw
   .data(); distvec_at()/shift_at() read one pair's 3-vector. For each pair the
   distance vector satisfies D == r[j] - r[i] + S @ cell, where S is the shift. */
struct NeighbourList {
    std::vector<index_t> first;    /* [npairs]      i indices */
    std::vector<index_t> secnd;    /* [npairs]      j indices */
    std::vector<real_t> distvec;   /* [3 * npairs]  distance vectors */
    std::vector<real_t> absdist;   /* [npairs]      absolute distances */
    std::vector<index_t> shift;    /* [3 * npairs]  cell shifts */
    index_t npairs = 0;

    Span<const index_t> first_view() const { return {first.data(), first.size()}; }
    Span<const index_t> secnd_view() const { return {secnd.data(), secnd.size()}; }
    Span<const real_t> distvec_view() const {
        return {distvec.data(), distvec.size()};
    }
    Span<const real_t> absdist_view() const {
        return {absdist.data(), absdist.size()};
    }
    Span<const index_t> shift_view() const { return {shift.data(), shift.size()}; }

    Vec3<real_t> distvec_at(index_t p) const {
        return {distvec[3 * p], distvec[3 * p + 1], distvec[3 * p + 2]};
    }
    Vec3<index_t> shift_at(index_t p) const {
        return {shift[3 * p], shift[3 * p + 1], shift[3 * p + 2]};
    }
};

/*
 * Validate the arguments shared by neighbour_list(), neighbour_matrix() and the
 * GPU builds: nat >= 0; cutoff finite and > 0; per-atom radii finite and >= 0;
 * per-type cutoffs finite with ncutoffs > 0 and a types array whose entries lie
 * in [0, ncutoffs). Returns NL_SUCCESS or NL_INVALID_ARGUMENT with a message.
 */
error_t validate_neighbour_args(index_t nat, real_t cutoff,
                                const real_t *per_atom_cutoff,
                                const real_t *per_type_cutoff_sq,
                                index_t ncutoffs, const index_t *types);

/*
 * Build a neighbour list by spatial binning. Output pairs are sorted by the
 * first index i; for each pair the distance vector D == r[j] - r[i] + S @ cell.
 *
 * quantities          bitwise-or of Quantity flags selecting outputs
 * cell_origin[3]      lower-left corner of the simulation cell
 * cell[9]             row-major 3x3 cell matrix (rows are lattice vectors)
 * inv_cell[9]         row-major inverse such that inv_cell . r = fractional
 * pbc[3]              periodicity per lattice direction
 * nat                 number of atoms
 * positions[3 * nat]  Cartesian atom positions
 * cutoff              binning / global cutoff (max pair cutoff for arrays)
 * per_atom_cutoff     [nat] per-atom radii (pair cutoff = sum), or nullptr
 * per_type_cutoff_sq  [ncutoffs * ncutoffs] squared per-type cutoffs, or nullptr
 * ncutoffs            dimension of the per-type matrix
 * types               [nat] per-atom type indices, or nullptr
 * out                 result (cleared and filled)
 *
 * Returns NL_SUCCESS; NL_INVALID_ARGUMENT (see validate_neighbour_args(), plus
 * non-finite positions); or NL_ERROR for a degenerate (zero-volume) cell. A
 * message is recorded via set_error() in both failure cases.
 */
error_t neighbour_list(int quantities, const real_t cell_origin[3],
                       const real_t cell[9], const real_t inv_cell[9],
                       const bool pbc[3], index_t nat, const real_t *positions,
                       real_t cutoff, const real_t *per_atom_cutoff,
                       const real_t *per_type_cutoff_sq, index_t ncutoffs,
                       const index_t *types, NeighbourList &out,
                       CellOrder order = CellOrder::Linear);

/*
 * Fixed-capacity ("dense") neighbour list: each atom's neighbours occupy a row
 * of an n x max_neighbours matrix, so the output shape is static (it depends
 * only on n and max_neighbours, not on the number of pairs). This suits array
 * frameworks that compile for fixed shapes (e.g. JAX): forces are a masked sum
 * over the neighbour axis, with no scatter.
 *
 *   idx   [n * max_neighbours]      neighbour atom indices, row-major
 *   dist  [n * max_neighbours * 3]  distance vectors D == r[j] - r[i] + S @ cell
 *                                   (if QUANTITY_DISTVEC is requested)
 *   shift [n * max_neighbours * 3]  cell shifts S (if QUANTITY_SHIFT is
 *                                   requested)
 *   count [n]                       true neighbour count per atom (may exceed
 *                                   max_neighbours)
 *   overflow                        true if any count > max_neighbours, i.e. the
 *                                   capacity was too small and rows are clipped;
 *                                   the caller should retry with a larger capacity
 *
 * Only the first min(count, max_neighbours) slots of a row are defined; the
 * contents of the unused slots are unspecified, so mask them with `count`.
 */
struct NeighbourMatrix {
    index_t n = 0;
    index_t max_neighbours = 0;
    std::vector<index_t> idx;
    std::vector<real_t> dist;
    std::vector<index_t> shift;
    std::vector<index_t> count;
    bool overflow = false;
};

/* Build the dense neighbour list (parameters as neighbour_list()). The indices
   and counts are always returned; `quantities` selects the per-slot extras,
   any combination of QUANTITY_DISTVEC and QUANTITY_SHIFT (others are
   ignored). */
error_t neighbour_matrix(const real_t cell_origin[3], const real_t cell[9],
                         const real_t inv_cell[9], const bool pbc[3], index_t nat,
                         const real_t *positions, real_t cutoff,
                         const real_t *per_atom_cutoff,
                         const real_t *per_type_cutoff_sq, index_t ncutoffs,
                         const index_t *types, index_t max_neighbours,
                         NeighbourMatrix &out, CellOrder order = CellOrder::Linear,
                         int quantities = QUANTITY_DISTVEC);

}  // namespace matscipy

#endif
