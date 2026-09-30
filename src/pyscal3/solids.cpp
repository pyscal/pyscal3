#include "system.h"
#include "parallel.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "string.h"
#include <chrono>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <map>
#include <string>
#include <any>

double get_number_from_bond(const int lm,
	const double* real_qi,
	const double* imag_qi,
	const double* real_qj,
	const double* imag_qj){

    double sum2ti,sum2tj;
    double realdotproduct,imgdotproduct;
    double connection;

    sum2ti = 0.0;
    sum2tj = 0.0;
    realdotproduct = 0.0;
    imgdotproduct = 0.0;

    for (int mi = 0; mi<2*lm+1 ; mi++){
        sum2ti += real_qi[mi]*real_qi[mi] + imag_qi[mi]*imag_qi[mi];
        sum2tj += real_qj[mi]*real_qj[mi] + imag_qj[mi]*imag_qj[mi];
        realdotproduct += real_qi[mi]*real_qj[mi];
        imgdotproduct  += imag_qi[mi]*imag_qj[mi];
    }

    connection = (realdotproduct+imgdotproduct)/(sqrt(sum2tj)*sqrt(sum2ti));
    return connection;
}

py::tuple calculate_bonds(const nl_index& offsets,
	const nl_index& neighbors,
	const nl_values& q_real,
	const nl_values& q_imag,
	const int lm,
	const double threshold,
	const double avgthreshold,
	const double minbonds,
	const int comparecriteria,
	const int criteria){

    // returns bonds, sij (per bond, in neighbour order), avg_sij and solid
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* qr = q_real.data();
    const double* qi = q_imag.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    const int nm = 2*lm + 1;

    py::array_t<double> bonds(nop), avg_sij(nop), solid(nop);
    py::array_t<double> sij(static_cast<py::ssize_t>(off[nop]));
    double* bo = bonds.mutable_data();
    double* av = avg_sij.mutable_data();
    double* so = solid.mutable_data();
    double* sj = sij.mutable_data();


    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    int frenkelcons;
    double scalar, tempsij;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        frenkelcons = 0;
        tempsij = 0.0;
        for (std::int64_t ci=off[ti]; ci<off[ti+1]; ci++){
            const std::int64_t tj = nb[ci];
            scalar = get_number_from_bond(lm, qr + ti*nm, qi + ti*nm, qr + tj*nm, qi + tj*nm);
            sj[ci] = scalar;
            if (comparecriteria == 0){
                if (scalar > threshold) frenkelcons += 1;
            }
            else{
                if (scalar < threshold) frenkelcons += 1;
            }
            tempsij += scalar;
        }
        bo[ti] = frenkelcons;
        av[ti] = tempsij/double(off[ti+1] - off[ti]);
    }
    }, 256);
    }

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    int issolid;
    double tfrac;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        if (criteria == 0){
            if (comparecriteria==0){
                issolid = ((bo[ti] > minbonds) && (av[ti] > avgthreshold));
            }
            else{
                issolid = ((bo[ti] > minbonds) && (av[ti] < avgthreshold));
            }
        }
        else {
            tfrac = (bo[ti]/double(off[ti+1] - off[ti]) > minbonds);
            if (comparecriteria==0){
                issolid = (tfrac && (av[ti] > avgthreshold));
            }
            else{
                issolid = (tfrac && (av[ti] < avgthreshold));
            }
        }
        so[ti] = issolid;
    }
    }, 256);
    }
    return py::make_tuple(bonds, sij, avg_sij, solid);
}

static void extract_cluster(std::int64_t seed,
	int clusterindex,
	const bool* condition,
	const std::int64_t* off,
	const std::int64_t* nb,
	const double* dist,
	const vector<double>& cutoff,
	std::int64_t* cluster){

	// depth-first search with an explicit stack: a recursive search
	// overflowed the call stack for clusters of a few 100 000 atoms
	vector<std::int64_t> stack(1, seed);
	while (!stack.empty()){
		const std::int64_t ti = stack.back();
		stack.pop_back();
		for (std::int64_t ci=off[ti]; ci<off[ti+1]; ci++){
			const std::int64_t cc = nb[ci];
			if (!condition[cc]) continue;
			if (!(dist[ci] <= cutoff[ti])) continue;
			if (cluster[cc] == -1){
				cluster[cc] = clusterindex;
				stack.push_back(cc);
			}
		}
	}
}

py::array_t<std::int64_t> find_clusters(const nl_index& offsets,
	const nl_index& neighbors,
	const nl_values& neighbordist,
	const nl_values& atom_cutoff,
	const py::array_t<bool, py::array::c_style | py::array::forcecast>& condition,
	double clustercutoff){

    // cluster id of every atom that satisfies condition, -1 otherwise;
    // a bond is followed if its length is at most the cutoff of the atom
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* dist = neighbordist.data();
    const bool* cond = condition.data();
    const py::ssize_t nop = offsets.shape(0) - 1;

    vector<double> cutoff(atom_cutoff.data(), atom_cutoff.data() + nop);
    if (clustercutoff != 0){
        fill(cutoff.begin(), cutoff.end(), clustercutoff);
    }

    py::array_t<std::int64_t> cluster(nop);
    std::int64_t* cl = cluster.mutable_data();
    fill(cl, cl + nop, -1);
    int clusterindex = 0;

    for (py::ssize_t ti=0; ti<nop; ti++){
        if (!cond[ti]) continue;
        if (cl[ti] == -1){
            clusterindex += 1;
            cl[ti] = clusterindex;
            extract_cluster(ti, clusterindex, cond, off, nb, dist, cutoff, cl);
        }
    }
    return cluster;
}
