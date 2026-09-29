/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 */

#ifndef MATSCIPY_CELL_LIST_HH
#define MATSCIPY_CELL_LIST_HH

#include <cstdint>
#include <vector>

#include "types.hh"

namespace matscipy {

/* Order in which atoms of different cells are laid out in `sorted_atom`.
   Linear is z-major (counting sort); Morton interleaves the cell coordinates
   (Z-curve) for tighter 3D spatial locality. Both group each cell contiguously;
   only the order of cells differs. */
enum class CellOrder { Linear, Morton };

/*
 * Atoms binned into an n1 x n2 x n3 grid. A CellGrid bundles two things:
 *
 *  - the *grid definition* — the spatial partition itself: cell matrix/origin,
 *    its inverse, the resolution n1/n2/n3, the per-cell edge vectors, the box
 *    widths `len`, and periodicity. This is what makes the cell indices and
 *    `sorted_atom` meaningful. The query *cutoff* is not part of it: the cutoff
 *    is a property of a neighbour search, and the grid only bounds the usable
 *    cutoff via its resolution.
 *  - the *binning* — atoms grouped per cell. Atoms of one cell form a contiguous
 *    slice of `sorted_atom`, via one of two backends:
 *      dense  — `cell_first[c]` / `cell_count[c]` indexed by linear cell index
 *               c = c1 + n1*(c2 + n2*c3). O(ncells) memory.
 *      sparse — an open-addressing hash table keyed by the 64-bit linear cell
 *               index. O(#non-empty cells) memory, for huge/sparse grids.
 */
struct CellGrid {
    /* Grid definition. */
    real_t origin[3] = {0, 0, 0};
    real_t cell[9] = {0};              /* row-major lattice vectors */
    real_t inv_cell[9] = {0};
    index_t n1 = 0, n2 = 0, n3 = 0;
    real_t bin1[3] = {0}, bin2[3] = {0}, bin3[3] = {0};  /* cell edge vectors */
    real_t len[3] = {0};               /* box width along each lattice direction */
    bool pbc[3] = {false, false, false};

    /* Binning. */
    bool sparse = false;
    std::vector<index_t> sorted_atom;  /* size nat */

    /* Dense backend. */
    std::vector<index_t> cell_first;   /* size ncells */
    std::vector<index_t> cell_count;   /* size ncells */

    /* Sparse backend (open addressing, power-of-two capacity). */
    std::vector<std::int64_t> hkey;    /* size cap, -1 = empty */
    std::vector<index_t> hfirst;       /* size cap */
    std::vector<index_t> hcount;       /* size cap */
    std::int64_t hmask = 0;            /* cap - 1 */
};

/* These leaf helpers (hash, Morton key, cell lookup) are host+device: one
   definition serves both the CPU loops and the GPU kernels. */

/* 64-bit hash used by both the sparse build and query. The table index is the
   low bits of this value, so every input bit must influence the low output
   bits: a bare multiply (Fibonacci hashing) only spreads the *high* bits, and
   cell keys that differ solely above the mask width (e.g. a wire along z in a
   power-of-two grid, key = c3 << 20) would all land in one probe chain. The
   xor-shift/multiply finaliser (splitmix64) mixes both ways. */
MATSCIPY_HD inline std::int64_t cell_hash(std::int64_t key) {
    std::uint64_t h = static_cast<std::uint64_t>(key) * 0x9E3779B97F4A7C15ull;
    h ^= h >> 32;
    h *= 0xBF58476D1CE4E5B9ull;
    h ^= h >> 29;
    return static_cast<std::int64_t>(h >> 1);  /* keep non-negative */
}

/* Upper bound on the cells per lattice direction. Keeps the cell coordinates
   within the 21 bits the Morton key can interleave (part1by2 below) and the
   linear cell index within 2^60, so neither can overflow. Bins are widened
   beyond the cutoff when the bound is hit, which stays correct. */
constexpr index_t MAX_CELLS_PER_DIRECTION = index_t(1) << 20;

/* Spread the low 21 bits of x so they occupy every third bit (for Morton). */
MATSCIPY_HD inline std::uint64_t part1by2(std::uint64_t x) {
    x &= 0x1fffffull;
    x = (x | (x << 32)) & 0x1f00000000ffffull;
    x = (x | (x << 16)) & 0x1f0000ff0000ffull;
    x = (x | (x << 8)) & 0x100f00f00f00f00full;
    x = (x | (x << 4)) & 0x10c30c30c30c30c3ull;
    x = (x | (x << 2)) & 0x1249249249249249ull;
    return x;
}

/* Morton (Z-curve) key interleaving three cell coordinates. */
MATSCIPY_HD inline std::uint64_t morton3(index_t a, index_t b, index_t c) {
    return part1by2(a) | (part1by2(b) << 1) | (part1by2(c) << 2);
}

/* Dense lookup: neighbour cell coords -> [b, e) slice of sorted_atom. */
struct DenseQuery {
    index_t n1, n2;
    const index_t *cell_first;
    const index_t *cell_count;
    MATSCIPY_HD inline void slice(index_t c1, index_t c2, index_t c3, index_t &b,
                                  index_t &e) const {
        index_t c = c1 + n1 * (c2 + n2 * c3);
        b = cell_first[c];
        e = b + cell_count[c];
    }
};

/* Sparse lookup: hash the 64-bit linear cell key, probe linearly. */
struct SparseQuery {
    index_t n1, n2;
    std::int64_t mask;
    const std::int64_t *hkey;
    const index_t *hfirst;
    const index_t *hcount;
    MATSCIPY_HD inline void slice(index_t c1, index_t c2, index_t c3, index_t &b,
                                  index_t &e) const {
        std::int64_t key = static_cast<std::int64_t>(c1) +
                           static_cast<std::int64_t>(n1) *
                               (static_cast<std::int64_t>(c2) +
                                static_cast<std::int64_t>(n2) * c3);
        std::int64_t h = cell_hash(key) & mask;
        while (hkey[h] != -1) {
            if (hkey[h] == key) {
                b = hfirst[h];
                e = b + hcount[h];
                return;
            }
            h = (h + 1) & mask;
        }
        b = 0;
        e = 0;  /* empty cell */
    }
};

/*
 * Fill the grid *definition* of `cg` from the cell geometry and a cutoff: the
 * resolution n1/n2/n3 (so each cell is ~one cutoff and only ±1 cell is scanned),
 * the edge vectors, box widths, and a copy of origin/cell/inv_cell/pbc. Computes
 * no binning. Shared by the CPU and GPU paths. Returns false on a degenerate
 * (zero-volume) cell. The grid bounds the usable query cutoff to <= `cutoff`.
 */
bool cell_grid_geometry(const real_t origin[3], const real_t cell[9],
                        const real_t inv_cell[9], const bool pbc[3],
                        real_t cutoff, CellGrid &cg);

/*
 * Bin atoms into the grid already described in `cg` (by cell_grid_geometry).
 * `raw_cell` holds the per-atom (un-wrapped) integer cell coordinates (3*nat,
 * original atom order) from position_to_cell_index(); periodic directions are
 * wrapped and non-periodic clamped here. The dense/sparse backend is chosen
 * automatically. Reads cg.n1/n2/n3 and cg.pbc; fills the binning fields.
 */
void build_cell_grid(const index_t *raw_cell, index_t nat, CellOrder order,
                     CellGrid &cg);

}  // namespace matscipy

#endif
