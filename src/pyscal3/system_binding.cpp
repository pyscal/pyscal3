#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>
#include <string>
#include "system.h"
#include <map>
#include <string>
#include <any>

namespace py = pybind11;
using namespace std;


PYBIND11_MODULE(csystem, m) {
    // Only functions that are called from the Python layer are bound.
    // Functions with reference out-parameters cannot return their results
    // to Python and are used internally by the C++ code only.
    py::options options;
    options.disable_function_signatures();
    m.def("get_distance_vector", &get_distance_vector);
    m.def("get_all_neighbors_normal", &get_all_neighbors_normal);
    m.def("get_all_neighbors_shell_normal", &get_all_neighbors_shell_normal);
    m.def("get_all_neighbors_cells", &get_all_neighbors_cells);
    m.def("get_all_neighbors_shell_cells", &get_all_neighbors_shell_cells);
    m.def("get_all_neighbors_bynumber", &get_all_neighbors_bynumber);
    m.def("get_all_neighbors_sann", &get_all_neighbors_sann);
    m.def("get_all_neighbors_adaptive", &get_all_neighbors_adaptive);
    m.def("calculate_q_single", &calculate_q_single);
    m.def("calculate_aq_single", &calculate_aq_single);
    m.def("calculate_w_single", &calculate_w_single);
    m.def("calculate_aw_single", &calculate_aw_single);
    m.def("calculate_disorder", &calculate_disorder);
    m.def("calculate_bonds", &calculate_bonds);
    m.def("find_clusters", &find_clusters);
    m.def("get_all_neighbors_voronoi", &get_all_neighbors_voronoi);
    m.def("clean_voronoi_vertices", &clean_voronoi_vertices);
    m.def("get_cna_neighbors", &get_cna_neighbors);
    m.def("get_acna_neighbors_cn12", &get_acna_neighbors_cn12);
    m.def("get_acna_neighbors_cn14", &get_acna_neighbors_cn14);
    m.def("identify_cn12", &identify_cn12);
    m.def("identify_cn14", &identify_cn14);
    m.def("identify_diamond_cna", &identify_diamond_cna);
    m.def("calculate_centrosymmetry", &calculate_centrosymmetry);
    m.def("calculate_entropy", &calculate_entropy);
    m.def("calculate_average_entropy", &calculate_average_entropy);
    m.def("calculate_chi_params", &calculate_chi_params);
    m.def("calculate_angular_criteria", &calculate_angular_criteria);
    m.def("calculate_voronoi_vector", &calculate_voronoi_vector);
    m.def("calculate_short_range_order", &calculate_short_range_order);
    m.def("calculate_average_disorder", &calculate_average_disorder);
    m.def("calculate_average_over_neighbors", &calculate_average_over_neighbors);

#ifdef VERSION_INFO
    m.attr("__version__") = VERSION_INFO;
#else
    m.attr("__version__") = "dev";
#endif
}