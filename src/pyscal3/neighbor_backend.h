/*
Internal interface of neighbor_backend.cpp for other C++ routines (CNA):
the geometry of a structure and per-atom rows of neighbour candidates.
*/
#ifndef PYSCAL_NEIGHBOR_BACKEND_H
#define PYSCAL_NEIGHBOR_BACKEND_H

#include <cstdint>
#include <vector>
#include <pybind11/numpy.h>

#include "parallel.h"

namespace nlb {

namespace py = pybind11;
using idx = std::int64_t;
using darray = py::array_t<double, py::array::c_style | py::array::forcecast>;

template <typename T>
using buffer = pyscal::buffer<T>;

struct Geometry {
    idx n;
    const double *positions;
    double cell[9];
    double inv_cell[9];
    bool pbc[3];
    double volume;
};

// per-atom rows of selected pairs
struct Rows {
    std::vector<idx> offsets;
    buffer<idx> j;
    buffer<double> d, v;
    explicit Rows(idx n) : offsets(n + 1, 0) {}
    void add(idx jj, double dd, const double *vv) {
        j.push_back(jj);
        d.push_back(dd);
        v.insert(v.end(), vv, vv + 3);
    }
};

// positions (n, 3), row-major cell (3, 3) with the lattice vectors as rows,
// and the periodicity of the three directions
Geometry make_geometry(const darray &positions, const darray &cell,
                       const std::vector<bool> &pbc);

// prefactor * (V / N)^(1/3) of the cell
double guess_radius(const Geometry &g, double prefactor);

// candidates with d <= guess, each row sorted by distance (rounded to
// 1e-10) and then by neighbour index; v holds r_i - r_j
Rows candidates(const Geometry &g, double guess);

// The first nneed entries of every row are those of candidates(g, r_full),
// found with less work: atoms with at least nneed candidates within r_small
// get only those, the other atoms all candidates within r_full.
Rows nearest_candidates(const Geometry &g, double r_small, double r_full, int nneed);

}  // namespace nlb

#endif
