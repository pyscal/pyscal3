/*
 * puremath.cpp – C++ implementations of formerly pure-Python descriptors:
 *   1. calculate_chi_params      (chi parameter vectors)
 *   2. calculate_angular_criteria (angular parameter A)
 *   3. calculate_voronoi_vector   (Voronoi n3,n4,n5,n6 vector)
 *   4. calculate_short_range_order (Warren-Cowley SRO)
 */

#include "system.h"
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
static inline double cosine_angle(const vector<double>& v1,
                                  const vector<double>& v2) {
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
void calculate_chi_params(py::dict& atoms) {

    vector<vector<int>> neighbors =
        atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    vector<vector<vector<double>>> diff =
        atoms[py::str("diff")].cast<vector<vector<vector<double>>>>();

    int nop = (int)neighbors.size();

    // 10 bin edges → 9 bins
    const double bins[10] = {-1.0, -0.945, -0.915, -0.755, -0.705,
                             -0.195, 0.195, 0.245, 0.795, 1.0};

    vector<vector<int>> chiparams(nop, vector<int>(9, 0));
    vector<vector<double>> cosines(nop);

    for (int ti = 0; ti < nop; ti++) {
        int nn = (int)diff[ti].size();
        // Reserve space for C(nn,2) cosine values
        int npairs = nn * (nn - 1) / 2;
        cosines[ti].reserve(npairs);

        // Compute all pairwise cosines
        for (int i = 0; i < nn; i++) {
            for (int j = i + 1; j < nn; j++) {
                double ct = cosine_angle(diff[ti][i], diff[ti][j]);
                cosines[ti].push_back(ct);
            }
        }

        // Histogram into bins
        // Match numpy behavior: bins 0..7 are [left, right), bin 8 is [left, right]
        for (double ct : cosines[ti]) {
            for (int b = 0; b < 9; b++) {
                if (b < 8) {
                    if (ct >= bins[b] && ct < bins[b + 1]) {
                        chiparams[ti][b]++;
                        break;
                    }
                } else {
                    // Last bin: closed on both sides [0.795, 1.0]
                    if (ct >= bins[b] && ct <= bins[b + 1]) {
                        chiparams[ti][b]++;
                        break;
                    }
                }
            }
        }
    }

    atoms[py::str("chiparams")] = chiparams;
    atoms[py::str("cosines")]   = cosines;
}


/* =======================================================================
 *  2. Angular Criteria
 *
 *  For each atom, take the 4 closest neighbors (by distance), compute
 *  C(4,2)=6 pairwise cosines, and sum (cos(theta) + 1/3)^2.
 * ======================================================================= */
void calculate_angular_criteria(py::dict& atoms) {

    vector<vector<double>> neighbordist =
        atoms[py::str("neighbordist")].cast<vector<vector<double>>>();
    vector<vector<vector<double>>> diff =
        atoms[py::str("diff")].cast<vector<vector<vector<double>>>>();

    int nop = (int)neighbordist.size();
    vector<double> angular(nop, 0.0);

    for (int ti = 0; ti < nop; ti++) {
        int nn = (int)neighbordist[ti].size();
        if (nn < 4) {
            angular[ti] = 0.0;
            continue;
        }

        // Find indices of 4 closest neighbors using stable sort
        // (matches numpy's argsort which is stable)
        vector<datom> dv(nn);
        for (int i = 0; i < nn; i++) {
            dv[i].dist  = neighbordist[ti][i];
            dv[i].index = i;
        }
        stable_sort(dv.begin(), dv.end(), by_dist());

        int top4[4] = { dv[0].index, dv[1].index, dv[2].index, dv[3].index };

        double costhetasum = 0.0;
        for (int i = 0; i < 4; i++) {
            for (int j = i + 1; j < 4; j++) {
                double ct = cosine_angle(diff[ti][top4[i]], diff[ti][top4[j]]);
                double term = ct + 1.0 / 3.0;
                costhetasum += term * term;
            }
        }
        angular[ti] = costhetasum;
    }

    atoms[py::str("angular")] = angular;
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
        int refined_edge_count = 0;

        for (int fi = 0; fi < (int)face_vertices[x].size(); fi++) {
            int vno = face_vertices[x][fi];

            // Collect vertex indices for this face
            // vphase = vertex_numbers[x][st : st+vno]
            vector<int> vphase(vno);
            for (int k = 0; k < vno; k++) {
                vphase[k] = vertex_numbers[x][st + k];
            }

            // Compute edge lengths
            double edge_sum = 0.0;
            vector<double> edge_lengths(vno);
            for (int i = -1; i < vno - 1; i++) {
                // Wrap: i=-1 → last vertex paired with first
                int ii = (i < 0) ? vno - 1 : i;
                int jj = i + 1;
                int vi_idx = vphase[ii] * 3;
                int vj_idx = vphase[jj] * 3;
                double dx = vertex_vectors[x][vi_idx]     - vertex_vectors[x][vj_idx];
                double dy = vertex_vectors[x][vi_idx + 1] - vertex_vectors[x][vj_idx + 1];
                double dz = vertex_vectors[x][vi_idx + 2] - vertex_vectors[x][vj_idx + 2];
                double elen = sqrt(dx*dx + dy*dy + dz*dz);
                edge_lengths[ii == vno - 1 ? 0 : ii + 1] = elen;
                // Actually, Python iterates i in range(-1, len(vphase)-1)
                // storing sequentially. Let's just store sequentially:
            }
            // Re-do more carefully to exactly match Python:
            // Python: for i in range(-1, len(vphase)-1):
            //   edgeln between vphase[i] and vphase[i+1]
            // i=-1: vphase[-1] (last) vs vphase[0]
            // i=0:  vphase[0] vs vphase[1]
            // ...
            // i=vno-2: vphase[vno-2] vs vphase[vno-1]
            // Total: vno edges
            edge_sum = 0.0;
            edge_lengths.resize(vno);
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

            // Normalise and check area cutoff
            if (refined_edge_count < (int)neighborweight[x].size() &&
                neighborweight[x][refined_edge_count] > area_cutoff) {
                // Count edges passing edge_cutoff
                int edgecount = 0;
                if (edge_sum > 0.0) {
                    for (int e = 0; e < vno; e++) {
                        if (edge_lengths[e] / edge_sum > edge_cutoff)
                            edgecount++;
                    }
                }
                refined_edge_count++;
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
 *  For each atom, count how many neighbors are of the reference type,
 *  compute local composition fraction, then SRO parameter.
 * ======================================================================= */
void calculate_short_range_order(py::dict& atoms,
                                 int reference_type,
                                 int compare_type) {

    vector<int> types =
        atoms[py::str("types")].cast<vector<int>>();
    vector<vector<int>> neighbors =
        atoms[py::str("neighbors")].cast<vector<vector<int>>>();

    int nop = (int)types.size();

    // Compute global composition
    map<int, int> type_counts;
    for (int i = 0; i < nop; i++) {
        type_counts[types[i]]++;
    }
    double total = (double)nop;
    double global_comp = 0.0;
    if (type_counts.count(reference_type)) {
        global_comp = type_counts[reference_type] / total;
    }

    vector<double> sro(nop, 0.0);

    for (int i = 0; i < nop; i++) {
        int nn = (int)neighbors[i].size();
        if (nn == 0) {
            sro[i] = 0.0;
            continue;
        }

        // Count neighbors of reference type
        int ref_count = 0;
        for (int j = 0; j < nn; j++) {
            if (types[neighbors[i][j]] == reference_type)
                ref_count++;
        }
        double lc = (double)ref_count / (double)nn;

        if (reference_type == compare_type) {
            if (global_comp < 1.0)
                sro[i] = (lc - global_comp) / (1.0 - global_comp);
            else
                sro[i] = 0.0;
        } else {
            if (global_comp > 0.0)
                sro[i] = 1.0 - (lc / global_comp);
            else
                sro[i] = 0.0;
        }
    }

    atoms[py::str("sro")] = sro;
}


/* =======================================================================
 *  5. Average Disorder over Neighbors
 *
 *  Same pattern as calculate_average_entropy:
 *  avg_disorder[i] = mean( disorder[i], disorder[j] for j in neighbors[i] )
 * ======================================================================= */
void calculate_average_disorder(py::dict& atoms) {
    vector<double> disorder =
        atoms[py::str("disorder")].cast<vector<double>>();
    vector<vector<int>> neighbors =
        atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    int nop = (int)neighbors.size();

    vector<double> avg_disorder(nop);
    for (int ti = 0; ti < nop; ti++) {
        double sum = disorder[ti];
        int nn = (int)neighbors[ti].size();
        for (int j = 0; j < nn; j++) {
            sum += disorder[neighbors[ti][j]];
        }
        avg_disorder[ti] = sum / (double)(nn + 1);
    }
    atoms[py::str("avg_disorder")] = avg_disorder;
}


/* =======================================================================
 *  6. Generic Average-Over-Neighbors for 1-D per-atom values
 *
 *  Returns a py::list of averaged values rather than writing to a fixed key,
 *  so the Python caller can store under whatever key it wants.
 * ======================================================================= */
py::list calculate_average_over_neighbors(py::dict& atoms,
                                          const vector<double>& values,
                                          bool include_self) {
    vector<vector<int>> neighbors =
        atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    int nop = (int)neighbors.size();

    py::list result(nop);
    for (int ti = 0; ti < nop; ti++) {
        double sum = include_self ? values[ti] : 0.0;
        int count  = include_self ? 1 : 0;
        int nn = (int)neighbors[ti].size();
        for (int j = 0; j < nn; j++) {
            sum += values[neighbors[ti][j]];
            count++;
        }
        result[ti] = count > 0 ? sum / (double)count : 0.0;
    }
    return result;
}
