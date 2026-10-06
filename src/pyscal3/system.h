#include <iostream>
#include <iostream>
#include <exception>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <cstdint>
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <string>
#include <any>

const double PI = 3.141592653589793;

namespace py = pybind11;
using namespace std;

// Flat neighbour data passed from Python: offsets (n + 1) into per-bond
// arrays, neighbour indices, and per-bond or per-atom values.
using nl_index = py::array_t<std::int64_t, py::array::c_style | py::array::forcecast>;
using nl_values = py::array_t<double, py::array::c_style | py::array::forcecast>;

/*-----------------------------------------------------
    Some utility objects
-----------------------------------------------------*/

struct datom{
    double dist;
    int  index;
};

//create another for the sorting algorithm
struct by_dist{
    bool operator()(datom const &datom1, datom const &datom2){
        return (datom1.dist < datom2.dist);
    }
};
/*-----------------------------------------------------
    Neighbor Methods
-----------------------------------------------------*/

double get_abs_distance(vector<double>, vector<double>,
    const int&, 
    const vector<vector<double>>&, 
    const vector<vector<double>>&, 
    const vector<double>&, 
    double&, double&, double&);

vector<double> get_distance_vector(vector<double> pos1, 
    vector<double> pos2, 
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box);

vector<double> remap_atom_into_box(vector<double> pos, 
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box);

vector<double> remap_and_displace_atom(vector<double> pos, 
    const int& triclinic, 
    const vector<vector<double>>& rot, 
    const vector<vector<double>>& rotinv,
    const vector<double>& box,
    const vector<double>& perturbation);

void convert_to_spherical_coordinates(double, double, double, 
    double&, double&, double&);

void get_all_neighbors_voronoi(py::dict& atoms,
    const double neighbordistance,
    const int triclinic,
    const vector<vector<double>> rot, 
    const vector<vector<double>> rotinv,
    const vector<double> box,
    const double face_area_exponent);

/*-----------------------------------------------------
    Neighbor search on matscipy-neighbours
    (neighbor_backend.cpp)
-----------------------------------------------------*/
using nl_positions = py::array_t<double, py::array::c_style | py::array::forcecast>;

py::dict nl_cutoff(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double rc);
py::dict nl_shell(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double dmin, double dmax);
py::dict nl_number(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double prefactor, int nmax, bool assign);
py::dict nl_adaptive(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double prefactor, int nlimit, double padding);
py::dict nl_sann(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double prefactor);


/*-----------------------------------------------------
    Steinhardt Methods
-----------------------------------------------------*/

void calculate_factors(const int lm, 
    vector<vector<double>>& alm,
    vector<vector<double>>& blm, 
    vector<vector<double>>& clm,
    vector<double>& dl,
    vector<double>& el);

vector<vector<double>> calculate_plm(const int lm,
    const double costheta, 
    const double sintheta);

double dfactorial(int l,
    int m);

vector<vector<vector<double>>> calculate_ylm(const int lm,
    const double costheta,
    const double sintheta,
    const double cosphi,
    const double sinphi);

vector<vector<vector<vector<double>>>> calculate_q_atom(const int lm,
    const vector<double>& theta,
    const vector<double>& phi);

void calculate_q(py::dict& atoms,
    const int lm);

vector<double> ylm_norms(const int l);

void ylm_all_m(const int l,
    const double theta,
    const double phi,
    const vector<double>& norm,
    double* ylm_real,
    double* ylm_imag);

py::tuple calculate_q_single(const nl_index& offsets,
    const nl_values& theta,
    const nl_values& phi,
    const nl_values& weights,
    const int lm);

py::array_t<double> calculate_aq_single(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm);

py::tuple calculate_w_single(const nl_values& q_real,
    const nl_values& q_imag,
    const int lm);

py::tuple calculate_aw_single(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm);

py::array_t<double> calculate_disorder(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm);

py::tuple calculate_bonds(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& q_real,
    const nl_values& q_imag,
    const int lm,
    const double threshold,
    const double avgthreshold,
    const double minbonds,
    const int comparecriteria,
    const int criteria);

py::array_t<std::int64_t> find_clusters(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& neighbordist,
    const nl_values& atom_cutoff,
    const py::array_t<bool, py::array::c_style | py::array::forcecast>& condition,
    double clustercutoff);

/*-----------------------------------------------------
    CNA Methods (cna.cpp)
-----------------------------------------------------*/
py::tuple cna_structure(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double prefactor, double lattice_constant, int nmin);
py::tuple diamond_structure_cna(const nl_positions& positions, const nl_positions& cell,
    const vector<bool>& pbc, double prefactor);
py::tuple ackland_jones_structure(const nl_positions& positions, const nl_positions& cell,
                                  const vector<bool>& pbc, double prefactor);

/*-----------------------------------------------------
    Other Methods
-----------------------------------------------------*/
py::array_t<double> calculate_centrosymmetry(const nl_index& offsets,
    const nl_values& diff, const int nmax);

/*-----------------------------------------------------
    Pure-math descriptors (chi, angular, voronoi-vec, SRO)
-----------------------------------------------------*/
py::tuple calculate_chi_params(const nl_index& offsets, const nl_values& diff);

py::array_t<double> calculate_angular_criteria(const nl_index& offsets,
    const nl_values& neighbordist, const nl_values& diff);

void calculate_voronoi_vector(py::dict& atoms,
    double edge_cutoff,
    double area_cutoff);

py::array_t<double> calculate_short_range_order(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_index& types,
    int reference_type,
    int compare_type);

/*-----------------------------------------------------
    Entropy Methods
-----------------------------------------------------*/
py::array_t<double> calculate_entropy(const nl_index& offsets,
    const nl_values& neighbordist,
    const nl_values& cutoff,
    double sigma,
    double rho,
    double rstart,
    double rstop,
    double h,
    double kb);

py::array_t<double> calculate_average_entropy(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& entropy);

/*-----------------------------------------------------
    Neighbor-averaging helpers (puremath.cpp)
-----------------------------------------------------*/
py::array_t<double> calculate_average_disorder(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& disorder);

py::tuple calculate_local_deformation(const nl_index& offsets_cur,
    const nl_index& neighbors_cur,
    const nl_values& diff_cur,
    const nl_index& offsets_ref,
    const nl_index& neighbors_ref,
    const nl_values& diff_ref);

py::array_t<double> calculate_average_over_neighbors(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& values,
    bool include_self);