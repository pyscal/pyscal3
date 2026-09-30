/*
 * puremath.cpp – C++ implementations of formerly pure-Python descriptors:
 *   1. calculate_chi_params      (chi parameter vectors)
 *   2. calculate_angular_criteria (angular parameter A)
 *   3. calculate_voronoi_vector   (Voronoi n3,n4,n5,n6 vector)
 *   4. calculate_short_range_order (Warren-Cowley SRO)
 */

#include "system.h"
#include "parallel.h"
#include <limits>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <map>
#include <string>

namespace py = pybind11;
using namespace std;

/* -----------------------------------------------------------------------
 *  Helper: cosine of angle between two 3-vectors stored as flat sub-arrays
 *  in diff[atom][neighbor_idx] = {dx, dy, dz}
 * ----------------------------------------------------------------------- */
static inline double cosine_angle(const double* v1,
                                  const double* v2) {
    double dot = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
    double m1  = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
    double m2  = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
    if (m1 > 0.0 && m2 > 0.0) {
        double ct = dot / (m1 * m2);
        // Clamp to [-1, 1] to handle floating-point rounding
        if (ct < -1.0) ct = -1.0;
        if (ct >  1.0) ct =  1.0;
        return ct;
    }
    return 0.0;
}


/* =======================================================================
 *  1. Chi Parameters
 *
 *  For each atom, compute pairwise cosine angles between ALL neighbor
 *  displacement vectors, then histogram into 9 bins.
 *  bins = [-1.0, -0.945, -0.915, -0.755, -0.705, -0.195,
 *           0.195, 0.245, 0.795, 1.0]
 * ======================================================================= */
py::tuple calculate_chi_params(const nl_index& offsets, const nl_values& diff) {

    // diff holds the neighbour vectors, (bonds, 3); bonds of atom i are
    // offsets[i] .. offsets[i + 1]. Returns the (n, 9) chi counts, the
    // cosines of all neighbour pairs of each atom, and offsets into them.
    const std::int64_t* off = offsets.data();
    const double* v = diff.data();
    const py::ssize_t nop = offsets.shape(0) - 1;

    // 10 bin edges → 9 bins
    const double bins[10] = {-1.0, -0.945, -0.915, -0.755, -0.705,
                             -0.195, 0.195, 0.245, 0.795, 1.0};

    py::array_t<std::int64_t> cos_offsets(nop + 1);
    std::int64_t* co = cos_offsets.mutable_data();
    co[0] = 0;
    for (py::ssize_t ti = 0; ti < nop; ti++) {
        const std::int64_t nn = off[ti + 1] - off[ti];
        co[ti + 1] = co[ti] + nn * (nn - 1) / 2;
    }

    py::array_t<std::int64_t> chiparams(vector<py::ssize_t>{nop, 9});
    py::array_t<double> cosines(static_cast<py::ssize_t>(co[nop]));
    std::int64_t* chi = chiparams.mutable_data();
    double* cs = cosines.mutable_data();
    fill(chi, chi + 9 * nop, 0);

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        // all pairwise cosines
        std::int64_t k = co[ti];
        for (std::int64_t i = off[ti]; i < off[ti + 1]; i++) {
            for (std::int64_t j = i + 1; j < off[ti + 1]; j++) {
                cs[k++] = cosine_angle(v + 3 * i, v + 3 * j);
            }
        }

        // Histogram into bins
        // Match numpy behavior: bins 0..7 are [left, right), bin 8 is [left, right]
        for (std::int64_t p = co[ti]; p < co[ti + 1]; p++) {
            const double ct = cs[p];
            for (int b = 0; b < 9; b++) {
                if (b < 8) {
                    if (ct >= bins[b] && ct < bins[b + 1]) {
                        chi[9 * ti + b]++;
                        break;
                    }
                } else {
                    // Last bin: closed on both sides [0.795, 1.0]
                    if (ct >= bins[b] && ct <= bins[b + 1]) {
                        chi[9 * ti + b]++;
                        break;
                    }
                }
            }
        }
    }
    }, 256);
    }
    return py::make_tuple(chiparams, cosines, cos_offsets);
}


/* =======================================================================
 *  2. Angular Criteria
 *
 *  For each atom, take the 4 closest neighbors (by distance), compute
 *  C(4,2)=6 pairwise cosines, and sum (cos(theta) + 1/3)^2.
 * ======================================================================= */
py::array_t<double> calculate_angular_criteria(const nl_index& offsets,
    const nl_values& neighbordist,
    const nl_values& diff) {

    const std::int64_t* off = offsets.data();
    const double* dist = neighbordist.data();
    const double* v = diff.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    py::array_t<double> angular(nop);
    double* out = angular.mutable_data();

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        const int nn = (int)(off[ti + 1] - off[ti]);
        if (nn < 4) {
            out[ti] = 0.0;
            continue;
        }

        // Find indices of 4 closest neighbors using stable sort
        // (matches numpy's argsort which is stable)
        vector<datom> dv(nn);
        for (int i = 0; i < nn; i++) {
            dv[i].dist  = dist[off[ti] + i];
            dv[i].index = i;
        }
        stable_sort(dv.begin(), dv.end(), by_dist());

        const std::int64_t top4[4] = { off[ti] + dv[0].index, off[ti] + dv[1].index,
                                       off[ti] + dv[2].index, off[ti] + dv[3].index };

        double costhetasum = 0.0;
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                double ct = cosine_angle(v + 3 * top4[i], v + 3 * top4[j]);
                double term = ct + 1.0 / 3.0;
                costhetasum += term * term;
            }
        }
        out[ti] = costhetasum;
    }
    }, 256);
    }
    return angular;
}


/* =======================================================================
 *  3. Voronoi Vector
 *
 *  For each atom, walk through Voronoi face data:
 *    face_vertices[atom]   – list of vertex counts per face
 *    vertex_numbers[atom]  – flat list [skip, v0,v1,..., skip, v0,v1,..., ...]
 *    vertex_vectors[atom]  – flat list of vertex coordinates [x,y,z, x,y,z, ...]
 *    neighborweight[atom]  – face area fractions
 *
 *  For each face, compute edge lengths, normalise, count edges passing
 *  edge_cutoff, filter faces by area_cutoff, bin edge-counts into [n3,n4,n5,n6].
 * ======================================================================= */
void calculate_voronoi_vector(py::dict& atoms,
                              double edge_cutoff,
                              double area_cutoff) {

    vector<vector<int>> face_vertices =
        atoms[py::str("face_vertices")].cast<vector<vector<int>>>();
    vector<vector<int>> vertex_numbers =
        atoms[py::str("vertex_numbers")].cast<vector<vector<int>>>();
    vector<vector<double>> vertex_vectors =
        atoms[py::str("vertex_vectors")].cast<vector<vector<double>>>();
    vector<vector<double>> neighborweight =
        atoms[py::str("neighborweight")].cast<vector<vector<double>>>();

    int nop = (int)face_vertices.size();
    vector<vector<int>> vorovectors(nop, vector<int>(4, 0));

    for (int x = 0; x < nop; x++) {
        int st = 1;   // starting index into vertex_numbers (skip first entry)

        for (int fi = 0; fi < (int)face_vertices[x].size(); fi++) {
            int vno = face_vertices[x][fi];

            // Collect vertex indices for this face
            // vphase = vertex_numbers[x][st : st+vno]
            vector<int> vphase(vno);
            for (int k = 0; k < vno; k++) {
                vphase[k] = vertex_numbers[x][st + k];
            }

            // Edge lengths: vno edges between consecutive vertices
            // (the last vertex is paired with the first).
            double edge_sum = 0.0;
            vector<double> edge_lengths(vno);
            for (int i = -1; i < vno - 1; i++) {
                int idx_i = (i < 0) ? vno - 1 : i;
                int idx_j = i + 1;
                int vi3 = vphase[idx_i] * 3;
                int vj3 = vphase[idx_j] * 3;
                double dx = vertex_vectors[x][vi3]     - vertex_vectors[x][vj3];
                double dy = vertex_vectors[x][vi3 + 1] - vertex_vectors[x][vj3 + 1];
                double dz = vertex_vectors[x][vi3 + 2] - vertex_vectors[x][vj3 + 2];
                double elen = sqrt(dx*dx + dy*dy + dz*dz);
                edge_lengths[i + 1] = elen;  // i=-1→0, i=0→1, etc.
                edge_sum += elen;
            }

            st += (vno + 1);

            // Skip faces whose (relative) area is below area_cutoff.
            // neighborweight holds the area fraction of face fi.
            if (fi < (int)neighborweight[x].size() &&
                neighborweight[x][fi] > area_cutoff) {
                // Count edges passing edge_cutoff
                int edgecount = 0;
                if (edge_sum > 0.0) {
                    for (int e = 0; e < vno; e++) {
                        if (edge_lengths[e] / edge_sum > edge_cutoff)
                            edgecount++;
                    }
                }
                // Bin into n3,n4,n5,n6
                if (edgecount >= 3 && edgecount <= 6) {
                    vorovectors[x][edgecount - 3]++;
                }
            }
        }
    }

    atoms[py::str("vorovector")] = vorovectors;
}


/* =======================================================================
 *  4. Short-Range Order (Warren-Cowley)
 *
 *  alpha_AB(i) = 1 - p_AB(i) / c_B  for atoms i of type A (reference_type),
 *  where p_AB(i) is the fraction of neighbours of i that are of type B
 *  (compare_type) and c_B is the global concentration of B.  Atoms that are
 *  not of the reference type (or have no neighbours) get NaN.
 * ======================================================================= */
py::array_t<double> calculate_short_range_order(const nl_index& offsets,
                                 const nl_index& neighbors,
                                 const nl_index& types,
                                 int reference_type,
                                 int compare_type) {

    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const std::int64_t* ty = types.data();
    const int nop = (int)(offsets.shape(0) - 1);
    const double nan = std::numeric_limits<double>::quiet_NaN();

    // global concentration of the compare type
    int n_compare = 0;
    for (int i = 0; i < nop; i++) {
        if (ty[i] == compare_type) n_compare++;
    }
    double c_compare = (nop > 0) ? (double)n_compare / (double)nop : 0.0;

    py::array_t<double> sro_array(nop);
    double* sro = sro_array.mutable_data();
    fill(sro, sro + nop, nan);

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    for (py::ssize_t i = begin_; i < end_; i++) {
        if (ty[i] != reference_type) continue;
        int nn = (int)(off[i + 1] - off[i]);
        if (nn == 0 || c_compare <= 0.0) continue;

        int cmp_count = 0;
        for (std::int64_t j = off[i]; j < off[i + 1]; j++) {
            if (ty[nb[j]] == compare_type)
                cmp_count++;
        }
        double p = (double)cmp_count / (double)nn;
        sro[i] = 1.0 - p / c_compare;
    }
    }, 256);
    }
    return sro_array;
}


/* =======================================================================
 *  5. Average Disorder over Neighbors
 *
 *  Same pattern as calculate_average_entropy:
 *  avg_disorder[i] = mean( disorder[i], disorder[j] for j in neighbors[i] )
 * ======================================================================= */
py::array_t<double> calculate_average_disorder(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& disorder) {
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* dis = disorder.data();
    const py::ssize_t nop = offsets.shape(0) - 1;

    py::array_t<double> avg_disorder(nop);
    double* out = avg_disorder.mutable_data();
    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        double sum = dis[ti];
        for (std::int64_t ci = off[ti]; ci < off[ti+1]; ci++) {
            sum += dis[nb[ci]];
        }
        out[ti] = sum / (double)(off[ti+1] - off[ti] + 1);
    }
    }, 256);
    }
    return avg_disorder;
}


/* =======================================================================
 *  6. Generic Average-Over-Neighbors for 1-D per-atom values
 *
 *  Returns a py::list of averaged values rather than writing to a fixed key,
 *  so the Python caller can store under whatever key it wants.
 * ======================================================================= */
py::array_t<double> calculate_average_over_neighbors(const nl_index& offsets,
                                          const nl_index& neighbors,
                                          const nl_values& values,
                                          bool include_self) {
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* val = values.data();
    const py::ssize_t nop = offsets.shape(0) - 1;

    py::array_t<double> result(nop);
    double* out = result.mutable_data();
    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        double sum = include_self ? val[ti] : 0.0;
        int count  = include_self ? 1 : 0;
        for (std::int64_t j = off[ti]; j < off[ti + 1]; j++) {
            sum += val[nb[j]];
            count++;
        }
        out[ti] = count > 0 ? sum / (double)count : 0.0;
    }
    }, 256);
    }
    return result;
}


/* =======================================================================
 *  Local deformation: atomic strain, D^2_min and slip vector
 *
 *  For each atom, the neighbours that appear in both the current and the
 *  reference neighbour list are paired (first occurrence of each index in
 *  either list). With X the reference and Y the current neighbour vectors
 *  of the pairs, the affine fit is F^T = (X^T X)^-1 X^T Y. Returns
 *    strain (n, 3, 3)  E = (F^T F - I) / 2, NaN for fewer than 3 pairs or
 *                      a singular X^T X
 *    d2min  (n,)       mean |Y - X F^T|^2 over the pairs, NaN likewise
 *    slip   (n, 3)     mean of Y - X over the pairs, NaN without pairs
 * ======================================================================= */
static void first_occurrences(const std::int64_t* nb, std::int64_t lo, std::int64_t hi,
                              vector<pair<std::int64_t, std::int64_t>>& out) {
    // (neighbour index, bond index) sorted by neighbour index, first bond only
    out.clear();
    for (std::int64_t c = lo; c < hi; c++) out.emplace_back(nb[c], c);
    stable_sort(out.begin(), out.end(),
                [](const pair<std::int64_t, std::int64_t>& a,
                   const pair<std::int64_t, std::int64_t>& b) { return a.first < b.first; });
    out.erase(unique(out.begin(), out.end(),
                     [](const pair<std::int64_t, std::int64_t>& a,
                        const pair<std::int64_t, std::int64_t>& b) { return a.first == b.first; }),
              out.end());
}

static bool invert3(const double m[9], double inv[9]) {
    const double c00 = m[4]*m[8] - m[5]*m[7];
    const double c01 = m[5]*m[6] - m[3]*m[8];
    const double c02 = m[3]*m[7] - m[4]*m[6];
    const double det = m[0]*c00 + m[1]*c01 + m[2]*c02;
    if (det == 0.0 || !std::isfinite(det)) return false;
    inv[0] = c00 / det;
    inv[1] = (m[2]*m[7] - m[1]*m[8]) / det;
    inv[2] = (m[1]*m[5] - m[2]*m[4]) / det;
    inv[3] = c01 / det;
    inv[4] = (m[0]*m[8] - m[2]*m[6]) / det;
    inv[5] = (m[2]*m[3] - m[0]*m[5]) / det;
    inv[6] = c02 / det;
    inv[7] = (m[1]*m[6] - m[0]*m[7]) / det;
    inv[8] = (m[0]*m[4] - m[1]*m[3]) / det;
    return true;
}

py::tuple calculate_local_deformation(const nl_index& offsets_cur,
                                      const nl_index& neighbors_cur,
                                      const nl_values& diff_cur,
                                      const nl_index& offsets_ref,
                                      const nl_index& neighbors_ref,
                                      const nl_values& diff_ref) {
    const std::int64_t* oc = offsets_cur.data();
    const std::int64_t* nc = neighbors_cur.data();
    const double* vc = diff_cur.data();
    const std::int64_t* orf = offsets_ref.data();
    const std::int64_t* nr = neighbors_ref.data();
    const double* vr = diff_ref.data();
    const py::ssize_t nop = offsets_cur.shape(0) - 1;
    const double nan = std::numeric_limits<double>::quiet_NaN();

    py::array_t<double> strain(vector<py::ssize_t>{nop, 3, 3});
    py::array_t<double> d2min(nop);
    py::array_t<double> slip(vector<py::ssize_t>{nop, 3});
    double* E = strain.mutable_data();
    double* D = d2min.mutable_data();
    double* S = slip.mutable_data();


    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    vector<pair<std::int64_t, std::int64_t>> cur, ref;
    vector<pair<std::int64_t, std::int64_t>> pairs;   // (bond in cur, bond in ref)
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        first_occurrences(nc, oc[ti], oc[ti + 1], cur);
        first_occurrences(nr, orf[ti], orf[ti + 1], ref);
        pairs.clear();
        size_t a = 0, b = 0;
        while (a < cur.size() && b < ref.size()) {
            if (cur[a].first < ref[b].first) a++;
            else if (ref[b].first < cur[a].first) b++;
            else { pairs.emplace_back(cur[a].second, ref[b].second); a++; b++; }
        }
        const size_t np = pairs.size();

        // slip vector
        if (np == 0) {
            for (int k = 0; k < 3; k++) S[3*ti + k] = nan;
        } else {
            double s[3] = {0.0, 0.0, 0.0};
            for (const auto& p : pairs)
                for (int k = 0; k < 3; k++) s[k] += vc[3*p.first + k] - vr[3*p.second + k];
            for (int k = 0; k < 3; k++) S[3*ti + k] = s[k] / double(np);
        }

        // affine fit
        double xtx[9] = {0}, xty[9] = {0}, inv[9];
        for (const auto& p : pairs) {
            const double* x = vr + 3*p.second;
            const double* y = vc + 3*p.first;
            for (int r = 0; r < 3; r++)
                for (int c = 0; c < 3; c++) {
                    xtx[3*r + c] += x[r]*x[c];
                    xty[3*r + c] += x[r]*y[c];
                }
        }
        if (np < 3 || !invert3(xtx, inv)) {
            for (int k = 0; k < 9; k++) E[9*ti + k] = nan;
            D[ti] = nan;
            continue;
        }
        double ft[9];   // F^T = (X^T X)^-1 X^T Y
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) {
                double v = 0.0;
                for (int k = 0; k < 3; k++) v += inv[3*r + k]*xty[3*k + c];
                ft[3*r + c] = v;
            }
        // C = F^T F, with F = (F^T)^T: C_rc = sum_k F_kr F_kc = sum_k ft[r][k] ft[c][k]
        for (int r = 0; r < 3; r++)
            for (int c = 0; c < 3; c++) {
                double v = 0.0;
                for (int k = 0; k < 3; k++) v += ft[3*r + k]*ft[3*c + k];
                E[9*ti + 3*r + c] = 0.5*(v - (r == c ? 1.0 : 0.0));
            }
        // D^2_min: residuals Y - X F^T
        double d2 = 0.0;
        for (const auto& p : pairs) {
            const double* x = vr + 3*p.second;
            const double* y = vc + 3*p.first;
            for (int c = 0; c < 3; c++) {
                double pred = x[0]*ft[c] + x[1]*ft[3 + c] + x[2]*ft[6 + c];
                double res = y[c] - pred;
                d2 += res*res;
            }
        }
        D[ti] = d2 / double(np);
    }
    }, 256);
    }
    return py::make_tuple(strain, d2min, slip);
}
