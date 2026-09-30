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
#include "neighbor_backend.h"
#include "parallel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <new>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace nlb {

using pyscal::parallel_for;

// block sizes for parallel_for: per-pair loops are cheap per element
constexpr idx PAIR_GRAIN = 1 << 14;
constexpr idx ATOM_GRAIN = 512;

// periodic images are searched up to this radius; the exact comparison of
// each method is applied afterwards
double query_radius(double r) { return r * (1.0 + 1e-9) + 1e-12; }


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

// radii, if given, are per-atom radii: the pair i, j is searched up to
// radii[i] + radii[j] (rmax must be the largest such sum)
Pairs query(const Geometry &g, double rmax, const double *radii = nullptr) {
    Pairs p;
    if (g.n == 0 || !(rmax > 0.0)) return p;
    matscipy::NeighbourList nl;
    const double origin[3] = {0.0, 0.0, 0.0};
    const int quantities = matscipy::QUANTITY_FIRST | matscipy::QUANTITY_SECOND |
                           matscipy::QUANTITY_DISTVEC;
    const matscipy::error_t status = matscipy::neighbour_list(
        quantities, origin, g.cell, g.inv_cell, g.pbc, g.n, g.positions,
        query_radius(rmax), radii, nullptr, 0, nullptr, nl);
    if (status == matscipy::NL_INVALID_ARGUMENT)
        throw std::invalid_argument(matscipy::error_string);
    if (status == matscipy::NL_OUT_OF_MEMORY) throw std::bad_alloc();
    if (status != matscipy::NL_SUCCESS) throw std::runtime_error(matscipy::error_string);

    const idx m = nl.npairs;
    p.i = std::move(nl.first);
    p.j = std::move(nl.secnd);
    p.v.resize(3 * m);
    p.d.resize(m);
    parallel_for(m, [&](idx begin, idx end) {
        for (idx k = begin; k < end; k++) {
            const double x = -nl.distvec[3 * k];
            const double y = -nl.distvec[3 * k + 1];
            const double z = -nl.distvec[3 * k + 2];
            p.v[3 * k] = x;
            p.v[3 * k + 1] = y;
            p.v[3 * k + 2] = z;
            p.d[k] = sqrt(x * x + y * y + z * z);
        }
    }, PAIR_GRAIN);
    return p;
}

// start of the pairs of each atom in p, which is sorted by i
vector<idx> pair_offsets(const Pairs &p, idx n) {
    vector<idx> off(n + 1, 0);
    for (size_t k = 0; k < p.i.size(); k++) off[p.i[k] + 1]++;
    for (idx a = 0; a < n; a++) off[a + 1] += off[a];
    return off;
}

// rows of the pairs k of each atom a for which keep(a, k) holds, in the
// order of p: count per atom, then fill each atom's slice
template <typename Keep>
Rows filter(const Pairs &p, idx n, Keep keep) {
    const vector<idx> off = pair_offsets(p, n);
    Rows rows(n);
    parallel_for(n, [&](idx begin, idx end) {
        for (idx a = begin; a < end; a++) {
            idx c = 0;
            for (idx k = off[a]; k < off[a + 1]; k++)
                if (keep(a, k)) c++;
            rows.offsets[a + 1] = c;
        }
    }, ATOM_GRAIN);
    for (idx a = 0; a < n; a++) rows.offsets[a + 1] += rows.offsets[a];
    const idx total = rows.offsets[n];
    rows.j.resize(total);
    rows.d.resize(total);
    rows.v.resize(3 * total);
    parallel_for(n, [&](idx begin, idx end) {
        for (idx a = begin; a < end; a++) {
            idx w = rows.offsets[a];
            for (idx k = off[a]; k < off[a + 1]; k++) {
                if (!keep(a, k)) continue;
                rows.j[w] = p.j[k];
                rows.d[w] = p.d[k];
                std::copy(&p.v[3 * k], &p.v[3 * k] + 3, &rows.v[3 * w]);
                w++;
            }
        }
    }, ATOM_GRAIN);
    return rows;
}


// rows of the pairs of p that satisfy keep(d), in the order of p
template <typename Keep>
Rows select(const Pairs &p, idx n, Keep keep) {
    return filter(p, n, [&](idx, idx k) { return keep(p.d[k]); });
}

// distances closer than this count as equal when candidate rows are sorted,
// so that neighbours of one shell come out in index order although their
// computed distances differ in the last bits
constexpr double TIE_DISTANCE = 1e-10;

// each row sorted by distance (to TIE_DISTANCE) and then by neighbour index
Rows sort_rows(const Rows &unsorted, idx n) {
    Rows rows(n);
    rows.offsets = unsorted.offsets;
    rows.j.resize(unsorted.j.size());
    rows.d.resize(unsorted.d.size());
    rows.v.resize(unsorted.v.size());
    parallel_for(n, [&](idx begin, idx end) {
        vector<idx> order;
        vector<long long> key;
        for (idx a = begin; a < end; a++) {
            const idx lo = unsorted.offsets[a], hi = unsorted.offsets[a + 1];
            order.resize(hi - lo);
            key.resize(hi - lo);
            std::iota(order.begin(), order.end(), lo);
            for (idx k = lo; k < hi; k++) key[k - lo] = std::llround(unsorted.d[k] / TIE_DISTANCE);
            std::stable_sort(order.begin(), order.end(), [&](idx x, idx y) {
                if (key[x - lo] != key[y - lo]) return key[x - lo] < key[y - lo];
                return unsorted.j[x] < unsorted.j[y];
            });
            idx w = lo;
            for (idx k : order) {
                rows.j[w] = unsorted.j[k];
                rows.d[w] = unsorted.d[k];
                std::copy(&unsorted.v[3 * k], &unsorted.v[3 * k] + 3, &rows.v[3 * w]);
                w++;
            }
        }
    }, ATOM_GRAIN);
    return rows;
}

Rows candidates(const Geometry &g, double guess) {
    const Pairs p = query(g, guess);
    return sort_rows(select(p, g.n, [guess](double d) { return d <= guess; }), g.n);
}

Rows nearest_candidates(const Geometry &g, double r_small, double r_full, int nneed) {
    Rows c = candidates(g, r_small);
    vector<char> full(g.n, 0);
    bool any = false;
    for (idx a = 0; a < g.n; a++) {
        if (c.offsets[a + 1] - c.offsets[a] < nneed) {
            full[a] = 1;
            any = true;
        }
    }
    if (!any) return c;
    c = Rows(0);
    // one more search: atom i needs neighbours up to r_i (r_full for the
    // atoms above, r_small otherwise); with radii R_i = r_full - r_small / 2
    // and r_small / 2, R_i + R_j is at least r_i for every pair
    vector<double> radii(g.n);
    for (idx a = 0; a < g.n; a++) radii[a] = full[a] ? r_full - 0.5 * r_small : 0.5 * r_small;
    const Pairs p = query(g, 2.0 * r_full - r_small, radii.data());
    const Rows rows = filter(p, g.n, [&](idx a, idx k) {
        return p.d[k] <= (full[a] ? r_full : r_small);
    });
    return sort_rows(rows, g.n);
}

// the first count[a] entries of each row of c
Rows prefix(const Rows &c, const vector<idx> &count) {
    const idx n = static_cast<idx>(count.size());
    Rows rows(n);
    for (idx a = 0; a < n; a++) rows.offsets[a + 1] = rows.offsets[a] + count[a];
    const idx total = rows.offsets[n];
    rows.j.resize(total);
    rows.d.resize(total);
    rows.v.resize(3 * total);
    parallel_for(n, [&](idx begin, idx end) {
        for (idx a = begin; a < end; a++) {
            const idx lo = c.offsets[a], w = rows.offsets[a];
            std::copy(&c.j[lo], &c.j[lo] + count[a], &rows.j[w]);
            std::copy(&c.d[lo], &c.d[lo] + count[a], &rows.d[w]);
            std::copy(&c.v[3 * lo], &c.v[3 * lo] + 3 * count[a], &rows.v[3 * w]);
        }
    }, ATOM_GRAIN);
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

// the neighbour data of find_neighbors, computed without the GIL
struct Result {
    Rows rows{0};
    vector<double> cutoff, r, theta, phi;
    Rows candidates{0};
    bool has_candidates = false;
    bool finished = true;
};

void add_angles(Result &res) {
    const Rows &rows = res.rows;
    const idx m = static_cast<idx>(rows.j.size());
    res.r.resize(m);
    res.theta.resize(m);
    res.phi.resize(m);
    parallel_for(m, [&](idx begin, idx end) {
        for (idx k = begin; k < end; k++) {
            const double x = rows.v[3 * k], y = rows.v[3 * k + 1], z = rows.v[3 * k + 2];
            convert_to_spherical_coordinates(x, y, z, res.r[k], res.phi[k], res.theta[k]);
            // r is the bond length; keep it identical to d, whatever the compiler
            // does with the two sqrt expressions
            res.r[k] = rows.d[k];
        }
    }, PAIR_GRAIN);
}

py::dict to_dict(Rows &&rows, vector<double> &&cutoff, vector<double> &&r,
                 vector<double> &&theta, vector<double> &&phi) {
    const idx n = static_cast<idx>(rows.offsets.size()) - 1;
    const idx m = static_cast<idx>(rows.j.size());
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

// cutoff value for atoms that got at least one neighbour, 0 otherwise
vector<double> cutoff_if_any(const Rows &rows, double value) {
    const idx n = static_cast<idx>(rows.offsets.size()) - 1;
    vector<double> cutoff(n, 0.0);
    for (idx a = 0; a < n; a++)
        if (rows.offsets[a + 1] > rows.offsets[a]) cutoff[a] = value;
    return cutoff;
}

// run compute(res) without the GIL, then build the returned dict
template <typename Compute>
py::dict run(Compute compute) {
    Result res;
    {
        py::gil_scoped_release release;
        compute(res);
        if (res.finished) add_angles(res);
    }
    if (!res.finished) {
        py::dict out;
        out["finished"] = false;
        return out;
    }
    py::dict out = to_dict(std::move(res.rows), std::move(res.cutoff), std::move(res.r),
                           std::move(res.theta), std::move(res.phi));
    if (res.has_candidates) add_candidates(out, std::move(res.candidates));
    return out;
}

}  // namespace nlb

using namespace nlb;

py::dict nl_cutoff(const darray &positions, const darray &cell,
                   const vector<bool> &pbc, double rc) {
    const Geometry g = make_geometry(positions, cell, pbc);
    return run([&](Result &res) {
        res.rows = select(query(g, rc), g.n, [rc](double d) { return d < rc; });
        res.cutoff = cutoff_if_any(res.rows, rc);
    });
}

py::dict nl_shell(const darray &positions, const darray &cell,
                  const vector<bool> &pbc, double dmin, double dmax) {
    const Geometry g = make_geometry(positions, cell, pbc);
    return run([&](Result &res) {
        res.rows = select(query(g, dmax), g.n,
                          [dmin, dmax](double d) { return d >= dmin && d <= dmax; });
        res.cutoff = cutoff_if_any(res.rows, dmax);
    });
}

py::dict nl_number(const darray &positions, const darray &cell,
                   const vector<bool> &pbc, double prefactor, int nmax, bool assign) {
    const Geometry g = make_geometry(positions, cell, pbc);
    bool enough = true;
    py::dict out = run([&](Result &res) {
        const double guess = guess_radius(g, prefactor);
        res.candidates = candidates(g, guess);
        res.has_candidates = true;
        const Rows &c = res.candidates;
        vector<idx> count(g.n, 0);
        res.cutoff.assign(g.n, 0.0);
        for (idx a = 0; a < g.n; a++) {
            if (c.offsets[a + 1] - c.offsets[a] < nmax) {
                enough = false;
                continue;
            }
            if (assign) {
                count[a] = nmax;
                res.cutoff[a] = guess;
            }
        }
        res.rows = prefix(c, count);
    });
    out["finished"] = enough;
    return out;
}

py::dict nl_adaptive(const darray &positions, const darray &cell,
                     const vector<bool> &pbc, double prefactor, int nlimit, double padding) {
    const Geometry g = make_geometry(positions, cell, pbc);
    return run([&](Result &res) {
        res.candidates = candidates(g, guess_radius(g, prefactor));
        res.has_candidates = true;
        const Rows &c = res.candidates;
        for (idx a = 0; a < g.n; a++) {
            if (c.offsets[a + 1] - c.offsets[a] < nlimit) {
                res.finished = false;
                return;
            }
        }
        res.cutoff.assign(g.n, 0.0);
        vector<idx> kept(g.n, 0);
        parallel_for(g.n, [&](idx begin, idx end) {
            for (idx a = begin; a < end; a++) {
                const idx lo = c.offsets[a], hi = c.offsets[a + 1];
                double summ = 0.0;
                for (idx k = lo; k < lo + nlimit; k++) summ += c.d[k];
                const double dcut = padding * (1.0 / double(nlimit)) * summ;
                res.cutoff[a] = dcut;
                for (idx k = lo; k < hi; k++)
                    if (c.d[k] < dcut) kept[a]++;
            }
        }, ATOM_GRAIN);
        // the kept candidates of a row are not always its first entries, so
        // copy them one by one
        Rows rows(g.n);
        for (idx a = 0; a < g.n; a++) rows.offsets[a + 1] = rows.offsets[a] + kept[a];
        const idx total = rows.offsets[g.n];
        rows.j.resize(total);
        rows.d.resize(total);
        rows.v.resize(3 * total);
        parallel_for(g.n, [&](idx begin, idx end) {
            for (idx a = begin; a < end; a++) {
                idx w = rows.offsets[a];
                for (idx k = c.offsets[a]; k < c.offsets[a + 1]; k++) {
                    if (!(c.d[k] < res.cutoff[a])) continue;
                    rows.j[w] = c.j[k];
                    rows.d[w] = c.d[k];
                    std::copy(&c.v[3 * k], &c.v[3 * k] + 3, &rows.v[3 * w]);
                    w++;
                }
            }
        }, ATOM_GRAIN);
        res.rows = std::move(rows);
    });
}

py::dict nl_sann(const darray &positions, const darray &cell,
                 const vector<bool> &pbc, double prefactor) {
    const Geometry g = make_geometry(positions, cell, pbc);
    return run([&](Result &res) {
        res.candidates = candidates(g, guess_radius(g, prefactor));
        res.has_candidates = true;
        const Rows &c = res.candidates;
        vector<idx> count(g.n, 0);
        vector<char> failed(g.n, 0);
        res.cutoff.assign(g.n, 0.0);
        parallel_for(g.n, [&](idx begin, idx end) {
            for (idx a = begin; a < end; a++) {
                const idx lo = c.offsets[a], hi = c.offsets[a + 1];
                const idx maxneighs = hi - lo;
                if (maxneighs < 3) {
                    failed[a] = 1;
                    continue;
                }
                idx m = 3;
                double summ = c.d[lo] + c.d[lo + 1] + c.d[lo + 2];
                double dcut = summ / double(m - 2);
                while (m < maxneighs && dcut >= c.d[lo + m]) {
                    summ += c.d[lo + m];
                    m++;
                    dcut = summ / double(m - 2);
                }
                if (m == maxneighs) {
                    failed[a] = 1;
                    continue;
                }
                count[a] = m;
                res.cutoff[a] = dcut;
            }
        }, ATOM_GRAIN);
        for (idx a = 0; a < g.n; a++) {
            if (failed[a]) {
                res.finished = false;
                return;
            }
        }
        res.rows = prefix(c, count);
    });
}
