/*
 * PTM wrapper for pyscal3.
 *
 * Provides calculate_ptm() which takes the standard pyscal3 atom dict
 * and runs Polyhedral Template Matching (Larsen et al. 2016) on every atom.
 *
 * The vendored PTM library lives in lib/ptm/ (MIT license).
 */

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>

#include "ptm_functions.h"
#include "ptm_constants.h"

namespace py = pybind11;

// ------------------------------------------------------------------
// Neighbor callback data structure
// ------------------------------------------------------------------
struct PyscalNbrData {
    const double *positions;   // (nop, 3) flat, absolute positions
    const int    *neighbors;   // (nop, max_nn) neighbor indices, -1 = unused
    const double *diffs;       // (nop, max_nn, 3) displacement vectors
    int nop;
    int max_nn;
};

// ------------------------------------------------------------------
// PTM neighbor callback
//
// The PTM library calls this to fetch neighbor positions for a given
// atom.  The first entry (index 0) must be the central atom itself
// (at the origin), and the rest are relative displacement vectors.
// Neighbors are returned sorted by distance (ascending).
// ------------------------------------------------------------------
struct NbrEntry {
    int index;
    double dx, dy, dz;
    double r2;
};

static bool cmp_nbr(const NbrEntry &a, const NbrEntry &b) {
    return a.r2 < b.r2;
}

static int pyscal_get_neighbours(
    void *vdata,
    size_t /* unused */,
    size_t atom_index,
    int num,            // max points to return (including central)
    int *ordering,      // not used for single-shell, can be NULL
    size_t *nbr_indices,
    int32_t *numbers,
    double (*nbr_pos)[3]
) {
    auto *data = static_cast<PyscalNbrData *>(vdata);
    int ai = static_cast<int>(atom_index);
    int max_nn = data->max_nn;

    // Central atom at the origin
    nbr_pos[0][0] = nbr_pos[0][1] = nbr_pos[0][2] = 0.0;
    nbr_indices[0] = atom_index;
    numbers[0] = 0;

    // Collect valid neighbors and sort by distance
    std::vector<NbrEntry> nbrs;
    nbrs.reserve(max_nn);

    for (int j = 0; j < max_nn; ++j) {
        int nidx = data->neighbors[ai * max_nn + j];
        if (nidx < 0) continue;

        NbrEntry e;
        e.index = nidx;
        e.dx = data->diffs[(ai * max_nn + j) * 3 + 0];
        e.dy = data->diffs[(ai * max_nn + j) * 3 + 1];
        e.dz = data->diffs[(ai * max_nn + j) * 3 + 2];
        e.r2 = e.dx * e.dx + e.dy * e.dy + e.dz * e.dz;
        nbrs.push_back(e);
    }

    std::sort(nbrs.begin(), nbrs.end(), cmp_nbr);

    int num_nbrs = std::min(num - 1, static_cast<int>(nbrs.size()));
    for (int j = 0; j < num_nbrs; ++j) {
        nbr_pos[j + 1][0] = nbrs[j].dx;
        nbr_pos[j + 1][1] = nbrs[j].dy;
        nbr_pos[j + 1][2] = nbrs[j].dz;
        nbr_indices[j + 1] = static_cast<size_t>(nbrs[j].index);
        numbers[j + 1] = 0;
    }

    return num_nbrs + 1;
}


// ------------------------------------------------------------------
// Main PTM analysis function exposed to Python
// ------------------------------------------------------------------
void calculate_ptm(py::dict &atoms, int32_t flags, double rmsd_cutoff) {
    // Get arrays from the pyscal3 dict
    auto positions_arr = atoms["positions"].cast<py::array_t<double>>();
    auto neighbors_arr = atoms["neighbors"].cast<py::array_t<int>>();
    auto diff_arr      = atoms["diff"].cast<py::array_t<double>>();

    auto pos_info = positions_arr.request();
    auto nbr_info = neighbors_arr.request();
    auto dif_info = diff_arr.request();

    int nop = static_cast<int>(pos_info.shape[0]);
    int max_nn = static_cast<int>(nbr_info.shape[1]);

    const double *positions = static_cast<const double *>(pos_info.ptr);
    const int    *neighbors = static_cast<const int *>(nbr_info.ptr);
    const double *diffs     = static_cast<const double *>(dif_info.ptr);

    // Prepare output arrays
    std::vector<int32_t> types(nop, 0);
    std::vector<double>  rmsds(nop, INFINITY);
    std::vector<double>  quats(nop * 4, 0.0);
    std::vector<double>  interatomic(nop, 0.0);
    std::vector<double>  lattice_constants(nop, 0.0);
    std::vector<int32_t> alloy_types(nop, 0);

    // Setup callback data
    PyscalNbrData nbrdata;
    nbrdata.positions = positions;
    nbrdata.neighbors = neighbors;
    nbrdata.diffs     = diffs;
    nbrdata.nop       = nop;
    nbrdata.max_nn    = max_nn;

    // Initialize PTM
    ptm_initialize_global();
    ptm_local_handle_t local_handle = ptm_initialize_local();

    double cutoff = (rmsd_cutoff <= 0) ? INFINITY : rmsd_cutoff;

    for (int i = 0; i < nop; ++i) {
        double scale, rmsd;
        double q[4] = {0, 0, 0, 0};
        int32_t type = PTM_MATCH_NONE, alloy = PTM_ALLOY_NONE;
        double iad = 0.0, lc = 0.0;

        ptm_index(local_handle,
                  static_cast<size_t>(i),
                  pyscal_get_neighbours,
                  static_cast<void *>(&nbrdata),
                  flags,
                  false,  // output_conventional_orientation
                  &type, &alloy, &scale, &rmsd, q,
                  NULL, NULL, NULL, NULL,  // F, F_res, U, P
                  &iad, &lc,
                  NULL, NULL,  // best_template_index, best_template
                  NULL);       // output_indices

        if (rmsd > cutoff) {
            type = PTM_MATCH_NONE;
            rmsd = INFINITY;
        }

        types[i] = type;
        rmsds[i] = rmsd;
        quats[i * 4 + 0] = q[0];
        quats[i * 4 + 1] = q[1];
        quats[i * 4 + 2] = q[2];
        quats[i * 4 + 3] = q[3];
        interatomic[i] = iad;
        lattice_constants[i] = lc;
        alloy_types[i] = alloy;
    }

    ptm_uninitialize_local(local_handle);

    // Store results back in the dict
    atoms[py::str("ptm_type")]       = types;
    atoms[py::str("ptm_rmsd")]       = rmsds;
    atoms[py::str("ptm_quat")]       = quats;
    atoms[py::str("ptm_interatomic_distance")] = interatomic;
    atoms[py::str("ptm_lattice_constant")]     = lattice_constants;
    atoms[py::str("ptm_alloy_type")] = alloy_types;
}
