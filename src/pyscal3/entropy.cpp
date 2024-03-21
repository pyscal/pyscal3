#include "system.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
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

double gmr(double r, 
	double sigma,
    double rho,
	int n_neighbors,
	vector<double> &neighbordist){

    double g = 0.00;
    double rij,r2;
    double sigma2 = sigma*sigma;
    double frho = 4.00*PI*rho*r*r;
    double fsigma = sqrt(2.00*PI*sigma2);
    double factor = (1.00/frho)*(1.00/fsigma);

    for(int i=0; i<n_neighbors; i++)
    {
            rij = neighbordist[i];
            r2 = (r-rij)*(r-rij);
            g+=exp((-1.00*r2)/(2.00*sigma2));
    }
    return factor*g;
}

//function which is to be integrated
double entropy_integrand(double r,
	double sigma,
    double rho,
	int n_neighbors,
	vector<double> &neighbordist){
        
    double g = gmr(r, sigma, rho, n_neighbors, neighbordist);
    return ((g*log(g)-g +1.00)*r*r);
}


double trapezoid_integration(const double rstart,
	const double rstop,
	const double h,
	double sigma,
    double rho,
	int n_neighbors,
	vector<double> &neighbordist,
	double kb){

    int nsteps = (rstop - rstart)/h;
    double summ=0.00;
    double xstart, xend;
    
    double rloop;
    double integral;

    xstart = entropy_integrand(rstart, sigma, rho, n_neighbors, neighbordist);

    for(int j=1; j<nsteps-1; j++)
    {
        rloop = rstart + j*h;
        summ += entropy_integrand(rloop, sigma, rho, n_neighbors, neighbordist);
        //cout<<'integrand'<<entropy_integrand(rloop, sigma, rho, n_neighbors, neighbordist)<<endl;
    }

    xend = entropy_integrand(rstart + nsteps*h, sigma, rho, n_neighbors, neighbordist);
    integral = (h/2.00)*(xstart + 2.00*summ + xend);

    integral = -1.*rho*kb*integral;
    return integral;
}

void calculate_entropy(py::dict& atoms, 
	double sigma, 
	double rho, 
	double rstart, 
	double rstop, 
	double h, 
	double kb){

    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    vector<double> cutoff = atoms[py::str("cutoff")].cast<vector<double>>();
    vector<vector<double>> neighbordist = atoms[py::str("neighbordist")].cast<vector<vector<double>>>();
	int nop = neighbors.size();

    vector<double> entropy(nop); 

    for(int ti=0; ti<nop; ti++){    
        if (rho == 0){
            rho = neighbors[ti].size()/(4.1887902047863905*pow(cutoff[ti], 3));
        }
        entropy[ti] = trapezoid_integration(rstart, rstop, h, sigma, rho, neighbors[ti].size(), neighbordist[ti], kb);
    }

    atoms[py::str("entropy")] = entropy;
}

void calculate_average_entropy(py::dict& atoms){
    double entsum;
    vector<double> entropy = atoms[py::str("entropy")].cast<vector<double>>();
    vector<vector<int>> neighbors = atoms[py::str("neighbors")].cast<vector<vector<int>>>();
	int nop = neighbors.size();

	vector<double> avg_entropy(nop);
    
    for(int ti=0; ti<nop; ti++){
        entsum = entropy[ti];
        for(int tj=0; tj<neighbors[ti].size(); tj++){
            entsum += entropy[neighbors[ti][tj]];
        }
        avg_entropy[ti] = entsum/(double(neighbors[ti].size() + 1));
    }

    atoms[py::str("average_entropy")] = avg_entropy;
}