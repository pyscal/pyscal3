#include "system.h"
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

void calculate_centrosymmetry(py::dict& atoms,
    const int nmax){
    
    double dx, dy, dz, weight;
    vector<datom> temp;

    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();;
    vector<vector<vector<double>>> diff = atoms[py::str("diff")].cast<vector<vector<vector<double>>>>();;

    int nop = neighbors.size();
    vector<double> centrosymmetry(nop);
    
    for (int ti=0; ti<nop; ti++){
    	int count = 0;
    	for (int i=0; i<neighbors[ti].size(); i++){
    		for (int j=i+1; j<neighbors[ti].size(); j++){
            	dx = diff[ti][i][0] + diff[ti][j][0];
            	dy = diff[ti][i][1] + diff[ti][j][1];
            	dz = diff[ti][i][2] + diff[ti][j][2];
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

	    centrosymmetry[ti] = csym;
    }
    atoms[py::str("centrosymmetry")] = centrosymmetry; 
}