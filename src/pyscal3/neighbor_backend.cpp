/*
Neighbour search on top of the vendored matscipy-neighbours core
(lib/matscipy-neighbours).

Each function takes the positions (n, 3), the row-major cell (3, 3), whose
rows are the lattice vectors, and the periodicity of the three directions,
and returns a dict of numpy arrays in CSR form:

    offsets (n + 1), j, d, r, theta, phi, weight (per pair), diff (pairs, 3),
    cutoff (n), and for the candidate-based methods temp_offsets, temp_j,
    temp_d, plus a bool "finished".

The conventions are the ones of the old search in neighbor.cpp:
diff = r_i - r_j (minimum image), theta = acos(z/r), phi = atan2(y, x) of
diff, fixed cutoff d < rc, shell dmin <= d <= dmax, candidates d <= guess.
Candidate rows are sorted by distance, rounded to TIE_DISTANCE, and then by
neighbour index.
*/
#include "system.h"
#include "neighbour_list.hh"
#include "error.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <new>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace {

using idx = std::int64_t;
using darray = py::array_t<double, py::array::c_style | py::array::forcecast>;

// periodic images are searched up to this radius; the exact comparison of
// each method is applied afterwards
double query_radius(double r) { return r * (1.0 + 1e-9) + 1e-12; }

struct Geometry {
    idx n;
    const double *positions;
    double cell[9];
    double inv_cell[9];
    bool pbc[3];
    double volume;
};

Geometry make_geometry(const darray &positions, const darray &cell,
                       const vector<bool> &pbc) {
    if (positions.ndim() != 2 || positions.shape(1) != 3)
        throw std::invalid_argument("positions must have shape (n, 3)");
    if (cell.ndim() != 2 || cell.shape(0) != 3 || cell.shape(1) != 3)
        throw std::invalid_argument("cell must have shape (3, 3)");
    if (pbc.size() != 3)
        throw std::invalid_argument("pbc must have three entries");

    Geometry g;
    g.n = positions.shape(0);
    g.positions = positions.data();
    const double *c = cell.data();
    for (int k = 0; k < 9; k++) g.cell[k] = c[k];
    for (int k = 0; k < 3; k++) g.pbc[k] = pbc[k];

    // inverse of cell^T, so that inv_cell . r gives fractional coordinates
    const double a = c[0], b = c[3], cc = c[6];
    const double d = c[1], e = c[4], f = c[7];
    const double gg = c[2], h = c[5], i = c[8];
    const double det = a * (e * i - f * h) - b * (d * i - f * gg) + cc * (d * h - e * gg);
    if (det == 0.0)
        throw std::invalid_argument("the cell vectors are linearly dependent");
    const double m[9] = {e * i - f * h, cc * h - b * i, b * f - cc * e,
                         f * gg - d * i, a * i - cc * gg, cc * d - a * f,
                         d * h - e * gg, b * gg - a * h, a * e - b * d};
    for (int k = 0; k < 9; k++) g.inv_cell[k] = m[k] / det;
    g.volume = std::fabs(det);
    return g;
}

// all pairs (i, j) with |r_i - r_j| up to about rmax, sorted by i;
// v holds r_i - r_j and d = |v|
struct Pairs {
    vector<idx> i, j;
    vector<double> d, v;
};

Pairs query(const Geometry &g, double rmax) {
    Pairs p;
    if (g.n == 0 || !(rmax > 0.0)) return p;
    matscipy::NeighbourList nl;
    const double origin[3] = {0.0, 0.0, 0.0};
    const int quantities = matscipy::QUANTITY_FIRST | matscipy::QUANTITY_SECOND |
                           matscipy::QUANTITY_DISTVEC;
    const matscipy::error_t status = matscipy::neighbour_list(
        quantities, origin, g.cell, g.inv_cell, g.pbc, g.n, g.positions,
        query_radius(rmax), nullptr, nullptr, 0, nullptr, nl);
    if (status == matscipy::NL_INVALID_ARGUMENT)
        throw std::invalid_argument(matscipy::error_string);
    if (status == matscipy::NL_OUT_OF_MEMORY) throw std::bad_alloc();
    if (status != matscipy::NL_SUCCESS) throw std::runtime_error(matscipy::error_string);

    const idx m = nl.npairs;
    p.i = std::move(nl.first);
    p.j = std::move(nl.secnd);
    p.v.resize(3 * m);
    p.d.resize(m);
    for (idx k = 0; k < m; k++) {
        const double x = -nl.distvec[3 * k];
        const double y = -nl.distvec[3 * k + 1];
        const double z = -nl.distvec[3 * k + 2];
        p.v[3 * k] = x;
        p.v[3 * k + 1] = y;
        p.v[3 * k + 2] = z;
        p.d[k] = sqrt(x * x + y * y + z * z);
    }
    return p;
}

// per-atom rows of selected pairs
struct Rows {
    vector<idx> offsets, j;
    vector<double> d, v;
    explicit Rows(idx n) : offsets(n + 1, 0) {}
    void add(idx jj, double dd, const double *vv) {
        j.push_back(jj);
        d.push_back(dd);
        v.insert(v.end(), vv, vv + 3);
    }
};

// rows of the pairs of p that satisfy keep(d), in the order of p
template <typename Keep>
Rows select(const Pairs &p, idx n, Keep keep) {
    Rows rows(n);
    for (size_t k = 0; k < p.i.size(); k++) {
        if (keep(p.d[k])) {
            rows.add(p.j[k], p.d[k], &p.v[3 * k]);
            rows.offsets[p.i[k] + 1]++;
        }
    }
    for (idx a = 0; a < n; a++) rows.offsets[a + 1] += rows.offsets[a];
    return rows;
}

// distances closer than this count as equal when candidate rows are sorted,
// so that neighbours of one shell come out in index order although their
// computed distances differ in the last bits
constexpr double TIE_DISTANCE = 1e-10;

// candidates with d <= guess, each row sorted by distance (to TIE_DISTANCE)
// and then by neighbour index
Rows candidates(const Geometry &g, double guess) {
    const Pairs p = query(g, guess);
    Rows unsorted = select(p, g.n, [guess](double d) { return d <= guess; });
    Rows rows(g.n);
    rows.offsets = unsorted.offsets;
    rows.j.reserve(unsorted.j.size());
    rows.d.reserve(unsorted.d.size());
    rows.v.reserve(unsorted.v.size());
    vector<long long> key(unsorted.d.size());
    for (size_t k = 0; k < key.size(); k++) key[k] = std::llround(unsorted.d[k] / TIE_DISTANCE);
    vector<idx> order;
    for (idx a = 0; a < g.n; a++) {
        const idx lo = unsorted.offsets[a], hi = unsorted.offsets[a + 1];
        order.resize(hi - lo);
        std::iota(order.begin(), order.end(), lo);
        std::stable_sort(order.begin(), order.end(), [&](idx x, idx y) {
            if (key[x] != key[y]) return key[x] < key[y];
            return unsorted.j[x] < unsorted.j[y];
        });
        for (idx k : order) rows.add(unsorted.j[k], unsorted.d[k], &unsorted.v[3 * k]);
    }
    return rows;
}

// the first count[a] entries of each row of c
Rows prefix(const Rows &c, const vector<idx> &count) {
    const idx n = static_cast<idx>(count.size());
    Rows rows(n);
    for (idx a = 0; a < n; a++) {
        const idx lo = c.offsets[a];
        for (idx k = lo; k < lo + count[a]; k++) rows.add(c.j[k], c.d[k], &c.v[3 * k]);
        rows.offsets[a + 1] = rows.offsets[a] + count[a];
    }
    return rows;
}

double guess_radius(const Geometry &g, double prefactor) {
    return prefactor * cbrt(g.volume / double(g.n));
}

template <typename T>
py::array_t<T> to_array(vector<T> &&values, vector<py::ssize_t> shape) {
    auto *owned = new vector<T>(std::move(values));
    py::capsule release(owned, [](void *ptr) { delete reinterpret_cast<vector<T> *>(ptr); });
    return py::array_t<T>(shape, owned->data(), release);
}

py::dict to_dict(Rows &&rows, vector<double> &&cutoff) {
    const idx n = static_cast<idx>(rows.offsets.size()) - 1;
    const idx m = static_cast<idx>(rows.j.size());
    vector<double> r(m), theta(m), phi(m);
    for (idx k = 0; k < m; k++) {
        const double x = rows.v[3 * k], y = rows.v[3 * k + 1], z = rows.v[3 * k + 2];
        convert_to_spherical_coordinates(x, y, z, r[k], phi[k], theta[k]);
        // r is the bond length; keep it identical to d, whatever the compiler
        // does with the two sqrt expressions
        r[k] = rows.d[k];
    }
    py::dict out;
    out["offsets"] = to_array(std::move(rows.offsets), {py::ssize_t(n + 1)});
    out["j"] = to_array(std::move(rows.j), {py::ssize_t(m)});
    out["d"] = to_array(std::move(rows.d), {py::ssize_t(m)});
    out["diff"] = to_array(std::move(rows.v), {py::ssize_t(m), py::ssize_t(3)});
    out["r"] = to_array(std::move(r), {py::ssize_t(m)});
    out["theta"] = to_array(std::move(theta), {py::ssize_t(m)});
    out["phi"] = to_array(std::move(phi), {py::ssize_t(m)});
    out["weight"] = to_array(vector<double>(m, 1.0), {py::ssize_t(m)});
    out["cutoff"] = to_array(std::move(cutoff), {py::ssize_t(n)});
    out["finished"] = true;
    return out;
}

void add_candidates(py::dict &out, Rows &&c) {
    const idx n = static_cast<idx>(c.offsets.size()) - 1;
    const idx m = static_cast<idx>(c.j.size());
    out["temp_offsets"] = to_array(std::move(c.offsets), {py::ssize_t(n + 1)});
    out["temp_j"] = to_array(std::move(c.j), {py::ssize_t(m)});
    out["temp_d"] = to_array(std::move(c.d), {py::ssize_t(m)});
}

py::dict not_finished() {
    py::dict out;
    out["finished"] = false;
    return out;
}

// cutoff value for atoms that got at least one neighbour, 0 otherwise
vector<double> cutoff_if_any(const Rows &rows, double value) {
    const idx n = static_cast<idx>(rows.offsets.size()) - 1;
    vector<double> cutoff(n, 0.0);
    for (idx a = 0; a < n; a++)
        if (rows.offsets[a + 1] > rows.offsets[a]) cutoff[a] = value;
    return cutoff;
}

}  // namespace

py::dict nl_cutoff(const darray &positions, const darray &cell,
                   const vector<bool> &pbc, double rc) {
    const Geometry g = make_geometry(positions, cell, pbc);
    Rows rows = select(query(g, rc), g.n, [rc](double d) { return d < rc; });
    vector<double> cutoff = cutoff_if_any(rows, rc);
    return to_dict(std::move(rows), std::move(cutoff));
}

py::dict nl_shell(const darray &positions, const darray &cell,
                  const vector<bool> &pbc, double dmin, double dmax) {
    const Geometry g = make_geometry(positions, cell, pbc);
    Rows rows = select(query(g, dmax), g.n,
                       [dmin, dmax](double d) { return d >= dmin && d <= dmax; });
    vector<double> cutoff = cutoff_if_any(rows, dmax);
    return to_dict(std::move(rows), std::move(cutoff));
}

py::dict nl_candidates(const darray &positions, const darray &cell,
                       const vector<bool> &pbc, double prefactor, int nmin) {
    const Geometry g = make_geometry(positions, cell, pbc);
    Rows c = candidates(g, guess_radius(g, prefactor));
    bool finished = true;
    for (idx a = 0; a < g.n; a++)
        if (c.offsets[a + 1] - c.offsets[a] < nmin) finished = false;
    py::dict out;
    add_candidates(out, std::move(c));
    out["finished"] = finished;
    return out;
}

py::dict nl_number(const darray &positions, const darray &cell,
                   const vector<bool> &pbc, double prefactor, int nmax, bool assign) {
    const Geometry g = make_geometry(positions, cell, pbc);
    const double guess = guess_radius(g, prefactor);
    Rows c = candidates(g, guess);
    vector<idx> count(g.n, 0);
    vector<double> cutoff(g.n, 0.0);
    bool finished = true;
    for (idx a = 0; a < g.n; a++) {
        if (c.offsets[a + 1] - c.offsets[a] < nmax) {
            finished = false;
            continue;
        }
        if (assign) {
            count[a] = nmax;
            cutoff[a] = guess;
        }
    }
    py::dict out = to_dict(prefix(c, count), std::move(cutoff));
    add_candidates(out, std::move(c));
    out["finished"] = finished;
    return out;
}

py::dict nl_adaptive(const darray &positions, const darray &cell,
                     const vector<bool> &pbc, double prefactor, int nlimit, double padding) {
    const Geometry g = make_geometry(positions, cell, pbc);
    Rows c = candidates(g, guess_radius(g, prefactor));
    Rows rows(g.n);
    vector<double> cutoff(g.n, 0.0);
    for (idx a = 0; a < g.n; a++) {
        const idx lo = c.offsets[a], hi = c.offsets[a + 1];
        if (hi - lo < nlimit) return not_finished();
        double summ = 0.0;
        for (idx k = lo; k < lo + nlimit; k++) summ += c.d[k];
        const double dcut = padding * (1.0 / double(nlimit)) * summ;
        cutoff[a] = dcut;
        idx kept = 0;
        for (idx k = lo; k < hi; k++) {
            if (c.d[k] < dcut) {
                rows.add(c.j[k], c.d[k], &c.v[3 * k]);
                kept++;
            }
        }
        rows.offsets[a + 1] = rows.offsets[a] + kept;
    }
    py::dict out = to_dict(std::move(rows), std::move(cutoff));
    add_candidates(out, std::move(c));
    return out;
}

py::dict nl_sann(const darray &positions, const darray &cell,
                 const vector<bool> &pbc, double prefactor) {
    const Geometry g = make_geometry(positions, cell, pbc);
    Rows c = candidates(g, guess_radius(g, prefactor));
    vector<idx> count(g.n, 0);
    vector<double> cutoff(g.n, 0.0);
    for (idx a = 0; a < g.n; a++) {
        const idx lo = c.offsets[a], hi = c.offsets[a + 1];
        const idx maxneighs = hi - lo;
        if (maxneighs < 3) return not_finished();
        idx m = 3;
        double summ = c.d[lo] + c.d[lo + 1] + c.d[lo + 2];
        double dcut = summ / double(m - 2);
        while (m < maxneighs && dcut >= c.d[lo + m]) {
            summ += c.d[lo + m];
            m++;
            dcut = summ / double(m - 2);
        }
        if (m == maxneighs) return not_finished();
        count[a] = m;
        cutoff[a] = dcut;
    }
    py::dict out = to_dict(prefix(c, count), std::move(cutoff));
    add_candidates(out, std::move(c));
    return out;
}
