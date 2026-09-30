/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 *
 * Binning neighbour-list construction as data-parallel primitives:
 *   - the per-atom cell index is computed once (cell_origin/inv_cell);
 *   - atoms are sorted into a CSR cell list (cell_list.{hh,cc}), dense or, for
 *     huge/sparse grids, a hashed compact table;
 *   - per-atom data is stored in cell-sorted order so the inner pair loop reads
 *     contiguous memory and does no matrix-vector products;
 *   - pairs are produced in two passes (count, then fill into exactly-sized
 *     buffers); each atom's work is independent and parallelised with OpenMP.
 */

#include "neighbour_list.hh"

/* pyscal3: the two pair loops of build_pairs run on the threads of
   pyscal3's parallel.h instead of OpenMP (see VENDORED.md). */
#include "parallel.h"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cmath>

#include "cell_list.hh"
#include "error.hh"
#include "neighbour_visit.hh"
#include "tools.hh"

namespace matscipy {

namespace {

/* The two passes (count, then fill) for a chosen cell backend. The traversal
   (visit_neighbours) and NeighbourContext in neighbour_visit.hh are shared with
   the GPU path. The first pass counts each atom's neighbours; the counts are
   prefix-summed into write offsets, the output buffers are sized exactly, and
   the second pass fills them. Pairs come out sorted by the first index i. */
template <typename Query>
void build_pairs(const NeighbourContext &ctx, const Query &q, index_t nat,
                 int quantities, NeighbourList &out) {
    std::vector<index_t> offset(nat + 1, 0);
    pyscal::parallel_for(nat, [&](index_t begin, index_t end) {
    for (index_t si = begin; si < end; si++) {
        index_t count = 0;
        visit_neighbours(ctx, q, si, [&count](index_t, const real_t *, real_t,
                                              const index_t *) { count++; });
        offset[ctx.sorted_atom[si] + 1] = count;
    }
    });

    for (index_t i = 0; i < nat; i++) offset[i + 1] += offset[i];
    const index_t npairs = offset[nat];
    out.npairs = npairs;

    const bool want_first = quantities & QUANTITY_FIRST;
    const bool want_secnd = quantities & QUANTITY_SECOND;
    const bool want_distvec = quantities & QUANTITY_DISTVEC;
    const bool want_absdist = quantities & QUANTITY_ABSDIST;
    const bool want_shift = quantities & QUANTITY_SHIFT;
    if (want_first) out.first.resize(npairs);
    if (want_secnd) out.secnd.resize(npairs);
    if (want_distvec) out.distvec.resize(3 * npairs);
    if (want_absdist) out.absdist.resize(npairs);
    if (want_shift) out.shift.resize(3 * npairs);

    pyscal::parallel_for(nat, [&](index_t begin, index_t end) {
    for (index_t si = begin; si < end; si++) {
        const index_t i = ctx.sorted_atom[si];
        index_t w = offset[i];
        visit_neighbours(
            ctx, q, si,
            [&](index_t sj, const real_t *dr, real_t abs_dr_sq,
                const index_t *shift) {
                if (want_first) out.first[w] = i;
                if (want_secnd) out.secnd[w] = ctx.sorted_atom[sj];
                if (want_distvec) {
                    out.distvec[3 * w + 0] = dr[0];
                    out.distvec[3 * w + 1] = dr[1];
                    out.distvec[3 * w + 2] = dr[2];
                }
                if (want_absdist) out.absdist[w] = std::sqrt(abs_dr_sq);
                if (want_shift) {
                    out.shift[3 * w + 0] = shift[0];
                    out.shift[3 * w + 1] = shift[1];
                    out.shift[3 * w + 2] = shift[2];
                }
                w++;
            });
    }
    });
}

}  // namespace

error_t validate_neighbour_args(index_t nat, real_t cutoff,
                                const real_t *per_atom_cutoff,
                                const real_t *per_type_cutoff_sq,
                                index_t ncutoffs, const index_t *types) {
    if (nat < 0) {
        return set_invalid_argument("Number of atoms must be non-negative.");
    }
    if (!(cutoff > 0) || !std::isfinite(cutoff)) {
        return set_invalid_argument(
            "Cutoff must be finite and positive (are all per-atom radii or "
            "per-type cutoffs zero?).");
    }
    if (per_atom_cutoff) {
        for (index_t a = 0; a < nat; a++) {
            if (!(per_atom_cutoff[a] >= 0) || !std::isfinite(per_atom_cutoff[a])) {
                return set_invalid_argumentf(
                    "Per-atom cutoff radius of atom %lld is not finite and "
                    "non-negative.",
                    static_cast<long long>(a));
            }
        }
    }
    if (per_type_cutoff_sq) {
        if (ncutoffs <= 0) {
            return set_invalid_argument(
                "Per-type cutoff matrix must have a positive dimension.");
        }
        for (index_t k = 0; k < ncutoffs * ncutoffs; k++) {
            if (!std::isfinite(per_type_cutoff_sq[k])) {
                return set_invalid_argument(
                    "Per-type cutoff matrix contains a non-finite entry.");
            }
        }
        if (!types) {
            return set_invalid_argument(
                "Per-type cutoffs require a per-atom type array.");
        }
        for (index_t a = 0; a < nat; a++) {
            if (types[a] < 0 || types[a] >= ncutoffs) {
                return set_invalid_argumentf(
                    "Type %lld of atom %lld is outside the %lldx%lld per-type "
                    "cutoff matrix.",
                    static_cast<long long>(types[a]), static_cast<long long>(a),
                    static_cast<long long>(ncutoffs),
                    static_cast<long long>(ncutoffs));
            }
        }
    }
    return NL_SUCCESS;
}

static error_t neighbour_list_body(int quantities, const real_t cell_origin[3],
                                   const real_t cell[9], const real_t inv_cell[9],
                                   const bool pbc[3], index_t nat, const real_t *r,
                                   real_t cutoff, const real_t *per_atom_cutoff,
                                   const real_t *per_type_cutoff_sq,
                                   index_t ncutoffs, const index_t *types,
                                   NeighbourList &out, CellOrder order) {

    out.first.clear();
    out.secnd.clear();
    out.distvec.clear();
    out.absdist.clear();
    out.shift.clear();
    out.npairs = 0;

    error_t status = validate_neighbour_args(nat, cutoff, per_atom_cutoff,
                                             per_type_cutoff_sq, ncutoffs, types);
    if (status != NL_SUCCESS) return status;

    /* The grid (definition + binning). Resolution/bins/len are derived from the
       cutoff; the cutoff itself stays a query parameter (below), not in the
       grid. No bin-count reduction here: a huge/sparse grid is handled by the
       compact (hashed) backend. */
    CellGrid cg;
    if (!cell_grid_geometry(cell_origin, cell, inv_cell, pbc, cutoff, cg)) {
        return set_error("Zero cell volume.");
    }
    const index_t n1 = cg.n1, n2 = cg.n2, n3 = cg.n3;

    /* Raw (un-wrapped) cell coordinate of every atom, computed once. A
       non-finite position would bin to garbage and (as NaN) defeat the
       distance test, so reject it here. */
    std::vector<index_t> raw(3 * nat);
    std::atomic<bool> all_finite_atoms(true);
    pyscal::parallel_for(nat, [&](index_t begin, index_t end) {
    bool all_finite = true;
    for (index_t a = begin; a < end; a++) {
        all_finite = all_finite && std::isfinite(r[3 * a]) &&
                     std::isfinite(r[3 * a + 1]) && std::isfinite(r[3 * a + 2]);
        position_to_cell_index(cell_origin, inv_cell, &r[3 * a], n1, n2, n3,
                               &raw[3 * a], &raw[3 * a + 1], &raw[3 * a + 2]);
    }
    if (!all_finite) all_finite_atoms.store(false);
    }, 4096);
    if (!all_finite_atoms.load()) {
        return set_invalid_argument("Positions must be finite.");
    }

    /* Bin the atoms (dense CSR, or hashed compact for huge/sparse grids). */
    build_cell_grid(raw.data(), nat, order, cg);

    /* Per-atom data in cell-sorted order (contiguous inner loop). */
    std::vector<index_t> raw_cell(3 * nat);
    std::vector<index_t> rel_cell(3 * nat);
    std::vector<real_t> rel_pos(3 * nat);
    std::vector<real_t> per_atom_s(per_atom_cutoff ? nat : 0);
    std::vector<index_t> types_s(types ? nat : 0);
    pyscal::parallel_for(nat, [&](index_t begin, index_t end) {
    for (index_t s = begin; s < end; s++) {
        const index_t a = cg.sorted_atom[s];
        const index_t c1 = raw[3 * a + 0], c2 = raw[3 * a + 1],
                      c3 = raw[3 * a + 2];
        raw_cell[3 * s + 0] = c1;
        raw_cell[3 * s + 1] = c2;
        raw_cell[3 * s + 2] = c3;

        const index_t r1 = pbc[0] ? c1 : bin_trunc(c1, n1);
        const index_t r2 = pbc[1] ? c2 : bin_trunc(c2, n2);
        const index_t r3 = pbc[2] ? c3 : bin_trunc(c3, n3);
        rel_cell[3 * s + 0] = r1;
        rel_cell[3 * s + 1] = r2;
        rel_cell[3 * s + 2] = r3;
        for (int k = 0; k < 3; k++) {
            rel_pos[3 * s + k] = r[3 * a + k] - r1 * cg.bin1[k] -
                                 r2 * cg.bin2[k] - r3 * cg.bin3[k];
        }
        if (per_atom_cutoff) per_atom_s[s] = per_atom_cutoff[a];
        if (types) types_s[s] = types[a];
    }
    }, 4096);

    NeighbourContext ctx;
    ctx.pbc0 = pbc[0];
    ctx.pbc1 = pbc[1];
    ctx.pbc2 = pbc[2];
    ctx.n1 = n1;
    ctx.n2 = n2;
    ctx.n3 = n3;
    ctx.nx = static_cast<index_t>(std::ceil(cutoff * n1 / cg.len[0]));
    ctx.ny = static_cast<index_t>(std::ceil(cutoff * n2 / cg.len[1]));
    ctx.nz = static_cast<index_t>(std::ceil(cutoff * n3 / cg.len[2]));
    for (int k = 0; k < 3; k++) {
        ctx.bin1[k] = cg.bin1[k];
        ctx.bin2[k] = cg.bin2[k];
        ctx.bin3[k] = cg.bin3[k];
    }
    ctx.cutoff_sq = cutoff * cutoff;
    ctx.per_type_cutoff_sq = per_type_cutoff_sq;
    ctx.ncutoffs = ncutoffs;
    ctx.sorted_atom = cg.sorted_atom.data();
    ctx.raw_cell = raw_cell.data();
    ctx.rel_cell = rel_cell.data();
    ctx.rel_pos = rel_pos.data();
    ctx.per_atom = per_atom_cutoff ? per_atom_s.data() : nullptr;
    ctx.types = types ? types_s.data() : nullptr;

    if (cg.sparse) {
        SparseQuery q{n1,           n2,           cg.hmask, cg.hkey.data(),
                      cg.hfirst.data(), cg.hcount.data()};
        build_pairs(ctx, q, nat, quantities, out);
    } else {
        DenseQuery q{n1, n2, cg.cell_first.data(), cg.cell_count.data()};
        build_pairs(ctx, q, nat, quantities, out);
    }

    return NL_SUCCESS;
}

error_t neighbour_list(int quantities, const real_t cell_origin[3],
                       const real_t cell[9], const real_t inv_cell[9],
                       const bool pbc[3], index_t nat, const real_t *r,
                       real_t cutoff, const real_t *per_atom_cutoff,
                       const real_t *per_type_cutoff_sq, index_t ncutoffs,
                       const index_t *types, NeighbourList &out,
                       CellOrder order) {
    clear_error();
    return catch_out_of_memory("Out of memory building the neighbour list.", [&] {
        return neighbour_list_body(quantities, cell_origin, cell, inv_cell, pbc,
                                   nat, r, cutoff, per_atom_cutoff,
                                   per_type_cutoff_sq, ncutoffs, types, out,
                                   order);
    });
}

static error_t neighbour_matrix_body(const real_t cell_origin[3],
                                     const real_t cell[9],
                                     const real_t inv_cell[9], const bool pbc[3],
                                     index_t nat, const real_t *positions,
                                     real_t cutoff, const real_t *per_atom_cutoff,
                                     const real_t *per_type_cutoff_sq,
                                     index_t ncutoffs, const index_t *types,
                                     index_t max_neighbours, NeighbourMatrix &out,
                                     CellOrder order, int quantities) {
    const index_t K = max_neighbours;
    const bool wD = quantities & QUANTITY_DISTVEC;
    const bool wS = quantities & QUANTITY_SHIFT;
    out.n = nat;
    out.max_neighbours = K;
    out.idx.clear();
    out.dist.clear();
    out.shift.clear();
    out.count.clear();
    out.overflow = false;

    error_t status = validate_neighbour_args(nat, cutoff, per_atom_cutoff,
                                             per_type_cutoff_sq, ncutoffs, types);
    if (status != NL_SUCCESS) return status;
    if (K < 0) {
        return set_invalid_argument("max_neighbours must be non-negative.");
    }

    /* std::vector value-initialises, so the unused slots happen to be 0 here;
       the contract leaves them unspecified (the GPU path does not clear). */
    out.idx.resize((size_t)nat * K);
    if (wD) out.dist.resize((size_t)nat * K * 3);
    if (wS) out.shift.resize((size_t)nat * K * 3);
    out.count.assign(nat, 0);
    if (nat <= 0) return NL_SUCCESS;

    /* The dense matrix is a reshape of the pair list: build the pairs, then
       scatter each into its atom's row. */
    NeighbourList nl;
    error_t e = neighbour_list_body(
        QUANTITY_FIRST | QUANTITY_SECOND | (wD ? QUANTITY_DISTVEC : 0) |
            (wS ? QUANTITY_SHIFT : 0),
        cell_origin, cell, inv_cell, pbc, nat, positions, cutoff, per_atom_cutoff,
        per_type_cutoff_sq, ncutoffs, types, nl, order);
    if (e != NL_SUCCESS) return e;

    /* CSR offsets from the i-sorted pair list, so rows fill in parallel. */
    std::vector<index_t> off(nat + 1, 0);
    for (index_t p = 0; p < nl.npairs; p++) off[nl.first[p] + 1]++;
    for (index_t a = 0; a < nat; a++) off[a + 1] += off[a];

    bool overflow = false;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(|| : overflow)
#endif
    for (index_t i = 0; i < nat; i++) {
        const index_t deg = off[i + 1] - off[i];
        out.count[i] = deg;
        if (deg > K) overflow = true;
        const index_t m = deg < K ? deg : K;
        for (index_t s = 0; s < m; s++) {
            const index_t p = off[i] + s;
            const size_t slot = (size_t)i * K + s;
            out.idx[slot] = nl.secnd[p];
            for (int k = 0; k < 3; k++) {
                if (wD) out.dist[slot * 3 + k] = nl.distvec[3 * p + k];
                if (wS) out.shift[slot * 3 + k] = nl.shift[3 * p + k];
            }
        }
    }
    out.overflow = overflow;
    return NL_SUCCESS;
}

error_t neighbour_matrix(const real_t cell_origin[3], const real_t cell[9],
                         const real_t inv_cell[9], const bool pbc[3], index_t nat,
                         const real_t *positions, real_t cutoff,
                         const real_t *per_atom_cutoff,
                         const real_t *per_type_cutoff_sq, index_t ncutoffs,
                         const index_t *types, index_t max_neighbours,
                         NeighbourMatrix &out, CellOrder order, int quantities) {
    clear_error();
    return catch_out_of_memory("Out of memory building the neighbour matrix.", [&] {
        return neighbour_matrix_body(cell_origin, cell, inv_cell, pbc, nat,
                                     positions, cutoff, per_atom_cutoff,
                                     per_type_cutoff_sq, ncutoffs, types,
                                     max_neighbours, out, order, quantities);
    });
}

}  // namespace matscipy
