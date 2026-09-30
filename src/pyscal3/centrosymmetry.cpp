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

py::array_t<double> calculate_centrosymmetry(const nl_index& offsets,
    const nl_values& diff,
    const int nmax){

    // sum of the nmax/2 smallest |r_ij + r_ik| squared over neighbour pairs
    const std::int64_t* off = offsets.data();
    const double* v = diff.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    py::array_t<double> centrosymmetry(nop);
    double* out = centrosymmetry.mutable_data();

    {
    py::gil_scoped_release release_gil;
    pyscal::parallel_for(nop, [&](std::int64_t begin_, std::int64_t end_) {
    double dx, dy, dz, weight;
    vector<datom> temp;
    for (py::ssize_t ti = begin_; ti < end_; ti++) {
        temp.clear();
        int count = 0;
        for (std::int64_t i=off[ti]; i<off[ti+1]; i++){
            for (std::int64_t j=i+1; j<off[ti+1]; j++){
                dx = v[3*i] + v[3*j];
                dy = v[3*i + 1] + v[3*j + 1];
                dz = v[3*i + 2] + v[3*j + 2];
                weight = sqrt(dx*dx+dy*dy+dz*dz);
                datom x = {weight, count};
                temp.emplace_back(x);
                count++;
            }
        }
        sort(temp.begin(), temp.end(), by_dist());

        double csym = 0;
        for(int i=0; i<nmax/2; i++){
            csym += temp[i].dist*temp[i].dist;
        }
        out[ti] = csym;
    }
    }, 256);
    }
    return centrosymmetry;
}