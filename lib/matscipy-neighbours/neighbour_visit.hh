/*
 * matscipy-neighbours — Neighbour list for particle simulations
 * https://github.com/libAtoms/matscipy-neighbours
 *
 * SPDX-License-Identifier: MIT
 * Copyright (2014-2026) James Kermode, University of Warwick
 *                       Lars Pastewka, University of Freiburg
 *                       and others (see toplevel AUTHORS file)
 *
 * The neighbour-cell traversal, shared by the CPU loops (neighbour_list.cc) and
 * the GPU kernels (neighbour_list_gpu.cc). It is one `__host__ __device__`
 * template so the distance/shift contract (D == r[j] - r[i] + S @ cell) lives in
 * exactly one place. The cell-lookup backend is a Query policy (dense or sparse,
 * from cell_list.hh); the per-pair sink is an Emit functor (a host lambda on the
 * CPU, a device functor in kernels).
 */

#ifndef MATSCIPY_NEIGHBOUR_VISIT_HH
#define MATSCIPY_NEIGHBOUR_VISIT_HH

#include "tools.hh"  /* bin_wrap / bin_trunc (host+device) */
#include "types.hh"

namespace matscipy {

/*
 * Per-atom data in cell-sorted order (index `s` is atom `sorted_atom[s]`) plus
 * the grid resolution and the cutoff-derived search range. Held by value (pbc as
 * three ints, not a pointer) so the struct is usable inside device kernels.
 */
struct NeighbourContext {
    index_t n1, n2, n3;
    index_t nx, ny, nz;          /* cells to scan each way (from the cutoff) */
    int pbc0, pbc1, pbc2;
    real_t bin1[3], bin2[3], bin3[3];
    real_t cutoff_sq;
    const real_t *per_type_cutoff_sq;
    index_t ncutoffs;
    const index_t *sorted_atom;  /* sorted position -> original atom index */
    const index_t *raw_cell;     /* 3*nat, sorted: raw (unwrapped) cell coord */
    const index_t *rel_cell;     /* 3*nat, sorted: trunc-if-non-periodic index */
    const real_t *rel_pos;       /* 3*nat, sorted: bin-relative position */
    const real_t *per_atom;      /* nat, sorted, or nullptr */
    const index_t *types;        /* nat, sorted, or nullptr */
};

/*
 * Visit every neighbour of the atom at sorted position `si`, invoking
 *     emit(s_j, dr[3], abs_dr_sq, shift[3]).
 */
template <typename Query, typename Emit>
MATSCIPY_HD inline void visit_neighbours(const NeighbourContext &c,
                                         const Query &q, index_t si,
                                         Emit &&emit) {
    const index_t n1 = c.n1, n2 = c.n2, n3 = c.n3;
    const index_t *raw_i = &c.raw_cell[3 * si];
    /* Cell index used to build dri (raw for periodic, truncated for
       non-periodic); the reference for the shift so it matches the distance
       vector. Using raw_i here would leak a spurious shift for a non-periodic
       atom outside the cell, where no boundary can be crossed. */
    const index_t *rel_i = &c.rel_cell[3 * si];
    const real_t *dri = &c.rel_pos[3 * si];

    const index_t ci1 = c.pbc0 ? bin_wrap(raw_i[0], n1) : bin_trunc(raw_i[0], n1);
    const index_t ci2 = c.pbc1 ? bin_wrap(raw_i[1], n2) : bin_trunc(raw_i[1], n2);
    const index_t ci3 = c.pbc2 ? bin_wrap(raw_i[2], n3) : bin_trunc(raw_i[2], n3);

    for (index_t z = -c.nz; z <= c.nz; z++) {
        index_t cj3 = ci3 + z;
        if (c.pbc2) cj3 = bin_wrap(cj3, n3);
        if (cj3 < 0 || cj3 >= n3) continue;
        real_t off3[3] = {z * c.bin3[0], z * c.bin3[1], z * c.bin3[2]};

        for (index_t y = -c.ny; y <= c.ny; y++) {
            index_t cj2 = ci2 + y;
            if (c.pbc1) cj2 = bin_wrap(cj2, n2);
            if (cj2 < 0 || cj2 >= n2) continue;
            real_t off2[3] = {off3[0] + y * c.bin2[0], off3[1] + y * c.bin2[1],
                              off3[2] + y * c.bin2[2]};

            for (index_t x = -c.nx; x <= c.nx; x++) {
                index_t cj1 = ci1 + x;
                if (c.pbc0) cj1 = bin_wrap(cj1, n1);
                if (cj1 < 0 || cj1 >= n1) continue;
                real_t off[3] = {off2[0] + x * c.bin1[0], off2[1] + x * c.bin1[1],
                                 off2[2] + x * c.bin1[2]};

                index_t begin, end;
                q.slice(cj1, cj2, cj3, begin, end);
                for (index_t sj = begin; sj < end; sj++) {
                    if (sj == si && x == 0 && y == 0 && z == 0) continue;

                    const real_t *drj = &c.rel_pos[3 * sj];
                    real_t dr[3] = {drj[0] - dri[0] + off[0],
                                    drj[1] - dri[1] + off[1],
                                    drj[2] - dri[2] + off[2]};
                    real_t abs_dr_sq =
                        dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2];
                    /* Written as !(d < cutoff) so a NaN distance (from a
                       non-finite position) is rejected, never accepted. */
                    if (!(abs_dr_sq < c.cutoff_sq)) continue;

                    bool inside_cutoff = true;
                    if (c.per_atom) {
                        real_t cc = c.per_atom[si] + c.per_atom[sj];
                        inside_cutoff = abs_dr_sq < cc * cc;
                    } else if (c.per_type_cutoff_sq && c.types) {
                        /* Types are validated against ncutoffs up front. */
                        index_t ti = c.types[si], tj = c.types[sj];
                        inside_cutoff =
                            abs_dr_sq < c.per_type_cutoff_sq[ti * c.ncutoffs + tj];
                    }
                    if (!inside_cutoff) continue;

                    const index_t *crj = &c.rel_cell[3 * sj];
                    index_t shift[3] = {(rel_i[0] - crj[0] + x) / n1,
                                        (rel_i[1] - crj[1] + y) / n2,
                                        (rel_i[2] - crj[2] + z) / n3};
                    emit(sj, dr, abs_dr_sq, shift);
                }
            }
        }
    }
}

}  // namespace matscipy

#endif
