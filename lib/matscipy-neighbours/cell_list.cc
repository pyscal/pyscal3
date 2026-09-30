/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#include "cell_list.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numeric>

#include "tools.hh"

namespace matscipy {

namespace {

/* Wrap (periodic) / clamp (non-periodic) a raw cell coordinate into [0, n). */
inline index_t fold(index_t c, index_t n, bool periodic) {
    return periodic ? bin_wrap(c, n) : bin_trunc(c, n);
}

/* Wrapped cell coordinates of every atom (used by both backends). */
void folded_coords(const index_t *raw_cell, const bool *pbc, index_t nat,
                   index_t n1, index_t n2, index_t n3,
                   std::vector<index_t> &cc) {
    cc.resize(3 * nat);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (index_t a = 0; a < nat; a++) {
        cc[3 * a + 0] = fold(raw_cell[3 * a + 0], n1, pbc[0]);
        cc[3 * a + 1] = fold(raw_cell[3 * a + 1], n2, pbc[1]);
        cc[3 * a + 2] = fold(raw_cell[3 * a + 2], n3, pbc[2]);
    }
}

void build_dense(const std::vector<index_t> &cc, index_t nat, index_t n1,
                 index_t n2, index_t n3, CellOrder order, CellGrid &cl) {
    const index_t ncells = n1 * n2 * n3;
    cl.cell_first.assign(ncells, 0);
    cl.cell_count.assign(ncells, 0);
    cl.sorted_atom.resize(nat);

    std::vector<index_t> lin(nat);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (index_t a = 0; a < nat; a++) {
        lin[a] = cc[3 * a] + n1 * (cc[3 * a + 1] + n2 * cc[3 * a + 2]);
    }

    /* Histogram of atoms per cell. */
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (index_t a = 0; a < nat; a++) {
#ifdef _OPENMP
#pragma omp atomic
#endif
        cl.cell_count[lin[a]]++;
    }

    if (order == CellOrder::Linear) {
        /* Prefix sum -> first index of each cell. */
        index_t acc = 0;
        for (index_t c = 0; c < ncells; c++) {
            cl.cell_first[c] = acc;
            acc += cl.cell_count[c];
        }
        /* Scatter atoms into their cell's slice. */
        std::vector<index_t> cursor(cl.cell_first);
        for (index_t a = 0; a < nat; a++) {
            cl.sorted_atom[cursor[lin[a]]++] = a;
        }
    } else {
        /* Morton order: sort atoms by Z-curve key, then record each cell's
           contiguous run. */
        std::vector<std::uint64_t> key(nat);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (index_t a = 0; a < nat; a++) {
            key[a] = morton3(cc[3 * a], cc[3 * a + 1], cc[3 * a + 2]);
        }
        std::iota(cl.sorted_atom.begin(), cl.sorted_atom.end(), 0);
        std::sort(cl.sorted_atom.begin(), cl.sorted_atom.end(),
                  [&](index_t x, index_t y) { return key[x] < key[y]; });
        /* Each cell is one contiguous run (the Morton key is unique per cell);
           record where each run starts. cell_count already holds the lengths. */
        for (index_t s = 0; s < nat; s++) {
            index_t c = lin[cl.sorted_atom[s]];
            if (s == 0 || lin[cl.sorted_atom[s - 1]] != c) {
                cl.cell_first[c] = s;
            }
        }
    }
}

void build_sparse(const std::vector<index_t> &cc, index_t nat, index_t n1,
                  index_t n2, CellGrid &cl) {
    /* Power-of-two capacity, load factor <= 0.5. */
    std::int64_t cap = 1;
    while (cap < 2 * static_cast<std::int64_t>(nat) + 1) cap <<= 1;
    cl.hmask = cap - 1;
    cl.hkey.assign(cap, -1);
    cl.hfirst.assign(cap, 0);
    cl.hcount.assign(cap, 0);
    cl.sorted_atom.resize(nat);

    auto key_of = [&](index_t a) -> std::int64_t {
        return static_cast<std::int64_t>(cc[3 * a]) +
               static_cast<std::int64_t>(n1) *
                   (static_cast<std::int64_t>(cc[3 * a + 1]) +
                    static_cast<std::int64_t>(n2) * cc[3 * a + 2]);
    };

    /* Insert keys and count occupants (serial; sparse is the rare path). */
    for (index_t a = 0; a < nat; a++) {
        std::int64_t key = key_of(a);
        std::int64_t h = cell_hash(key) & cl.hmask;
        while (cl.hkey[h] != -1 && cl.hkey[h] != key) h = (h + 1) & cl.hmask;
        cl.hkey[h] = key;
        cl.hcount[h]++;
    }

    /* Assign each occupied slot its slice start (prefix sum over slots). */
    index_t acc = 0;
    for (std::int64_t h = 0; h <= cl.hmask; h++) {
        if (cl.hkey[h] != -1) {
            cl.hfirst[h] = acc;
            acc += cl.hcount[h];
        }
    }

    /* Scatter atoms. */
    std::vector<index_t> cursor(cl.hfirst);
    for (index_t a = 0; a < nat; a++) {
        std::int64_t key = key_of(a);
        std::int64_t h = cell_hash(key) & cl.hmask;
        while (cl.hkey[h] != key) h = (h + 1) & cl.hmask;
        cl.sorted_atom[cursor[h]++] = a;
    }
}

}  // namespace

bool cell_grid_geometry(const real_t origin[3], const real_t cell[9],
                        const real_t inv_cell[9], const bool pbc[3],
                        real_t cutoff, CellGrid &cg) {
    const real_t *c1 = &cell[0], *c2 = &cell[3], *c3 = &cell[6];
    real_t nrm1[3], nrm2[3], nrm3[3];
    cross_product(c2, c3, nrm1);
    cross_product(c3, c1, nrm2);
    cross_product(c1, c2, nrm3);
    real_t volume = std::fabs(c3[0] * nrm3[0] + c3[1] * nrm3[1] + c3[2] * nrm3[2]);
    /* Degenerate if the volume is tiny relative to the product of the edge
       lengths (i.e. the lattice vectors are nearly coplanar), zero, or not
       finite. The ratio is scale-free, so the test is independent of units and
       of the cell's aspect ratio. */
    const real_t edges = norm(c1) * norm(c2) * norm(c3);
    if (!(volume > 1e-12 * edges) || !std::isfinite(volume)) return false;

    cg.len[0] = volume / norm(nrm1);
    cg.len[1] = volume / norm(nrm2);
    cg.len[2] = volume / norm(nrm3);
    /* Cells per direction: one cutoff wide, at least 1, and at most
       MAX_CELLS_PER_DIRECTION (clamped in floating point before the cast). */
    auto resolution = [cutoff](real_t len) -> index_t {
        const real_t n = std::floor(len / cutoff);
        if (!(n > 1.0)) return 1;
        if (n >= static_cast<real_t>(MAX_CELLS_PER_DIRECTION))
            return MAX_CELLS_PER_DIRECTION;
        return static_cast<index_t>(n);
    };
    cg.n1 = resolution(cg.len[0]);
    cg.n2 = resolution(cg.len[1]);
    cg.n3 = resolution(cg.len[2]);
    for (int k = 0; k < 3; k++) {
        cg.origin[k] = origin[k];
        cg.pbc[k] = pbc[k];
        cg.bin1[k] = c1[k] / cg.n1;
        cg.bin2[k] = c2[k] / cg.n2;
        cg.bin3[k] = c3[k] / cg.n3;
    }
    for (int k = 0; k < 9; k++) {
        cg.cell[k] = cell[k];
        cg.inv_cell[k] = inv_cell[k];
    }
    return true;
}

void build_cell_grid(const index_t *raw_cell, index_t nat, CellOrder order,
                     CellGrid &cg) {
    const index_t n1 = cg.n1, n2 = cg.n2, n3 = cg.n3;

    std::vector<index_t> cc;
    folded_coords(raw_cell, cg.pbc, nat, n1, n2, n3, cc);

    /* Choose the backend: dense unless the grid is far larger than the atom
       count (or would overflow a 32-bit cell index). */
    const std::int64_t ncells = static_cast<std::int64_t>(n1) * n2 * n3;
    const std::int64_t budget =
        std::max<std::int64_t>(1 << 20, 8 * static_cast<std::int64_t>(nat));
    cg.sparse = ncells > budget;

    if (!cg.sparse) {
        build_dense(cc, nat, n1, n2, n3, order, cg);
    } else {
        build_sparse(cc, nat, n1, n2, cg);
    }
}

}  // namespace matscipy
