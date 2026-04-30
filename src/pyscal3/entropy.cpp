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

/*
 * Optimised entropy calculation.
 *
 * Key optimisation: Gaussian windowing.
 * Each neighbour's Gaussian contribution is only non-negligible within
 * ±4σ of its distance.  By accumulating only within that window
 * we reduce the inner-loop work from  nsteps × n_neighbors  to
 * n_neighbors × (8σ/h).  With σ=0.2 and h=0.001 this is a ~3× reduction
 * in exp() calls, plus improved cache locality.
 *
 * The original trapezoidal rule quirk (skipping j=nsteps-1) and the
 * local-density bug (using atom 0's density for all atoms when rho==0)
 * are preserved for backward compatibility.
 */

void calculate_entropy(py::dict& atoms,
    double sigma,
    double rho,
    double rstart,
    double rstop,
    double h,
    double kb){

    vector<vector<int>> neighbors =
        atoms[py::str("neighbors")].cast<vector<vector<int>>>();
    vector<double> cutoff_vec =
        atoms[py::str("cutoff")].cast<vector<double>>();
    vector<vector<double>> neighbordist =
        atoms[py::str("neighbordist")].cast<vector<vector<double>>>();
    int nop = (int)neighbors.size();

    /* Match original: when rho==0 (local mode), set rho from atom 0
       and use for all atoms (original code modifies the rho parameter
       in first iteration and never resets it). */
    if (rho == 0) {
        rho = neighbors[0].size() /
            (4.1887902047863905 * pow(cutoff_vec[0], 3));
    }

    int nsteps = (int)((rstop - rstart) / h);

    /* Pre-compute constants */
    double sigma2      = sigma * sigma;
    double inv_2sigma2 = -1.0 / (2.0 * sigma2);          // negative for exp arg
    double inv_fsigma  =  1.0 / sqrt(2.0 * PI * sigma2);
    double gauss_window = 5.0 * sigma;                    // ≈ 1.0 Å for σ=0.2

    /* Pre-compute the r-grid (shared across atoms) */
    vector<double> r_grid(nsteps + 1);
    for (int j = 0; j <= nsteps; j++)
        r_grid[j] = rstart + j * h;

    vector<double> entropy_out(nop);
    vector<double> raw_g(nsteps + 1);   // per-atom, reused

    for (int ti = 0; ti < nop; ti++) {

        int nn = (int)neighbordist[ti].size();

        /* ---- Accumulate raw Gaussian sum with windowing ---- */
        fill(raw_g.begin(), raw_g.end(), 0.0);

        for (int i = 0; i < nn; i++) {
            double rij = neighbordist[ti][i];
            int jmin = max(0,      (int)((rij - gauss_window - rstart) / h));
            int jmax = min(nsteps, (int)((rij + gauss_window - rstart) / h) + 1);
            for (int j = jmin; j <= jmax; j++) {
                double dr = r_grid[j] - rij;
                raw_g[j] += exp(dr * dr * inv_2sigma2);
            }
        }

        /* ---- Trapezoidal integration (matches original loop structure) ----
         *
         *  Original code: xstart at j=0,
         *                 summ   for j = 1 .. nsteps-2,
         *                 xend   at j = nsteps.
         *  We replicate that exactly using half-weight endpoints
         *  and skipping j = nsteps-1.
         */
        double summ = 0.0;

        for (int j = 0; j <= nsteps; j++) {
            double r   = r_grid[j];
            double r2  = r * r;
            double frho = 4.0 * PI * rho * r2;

            /* g(r) = raw_g / (4π·ρ·r²·√(2πσ²)) */
            double g = (frho > 0.0) ? inv_fsigma * raw_g[j] / frho : 0.0;

            /* Integrand: (g·ln(g) - g + 1)·r²
               When g is negligibly small, limit → r². */
            double integrand;
            if (g > 1e-30)
                integrand = (g * log(g) - g + 1.0) * r2;
            else
                integrand = r2;

            /* Weight: half for endpoints, full for interior,
               but skip j = nsteps-1 to match original code. */
            if (j == 0 || j == nsteps)
                summ += 0.5 * integrand;
            else if (j != nsteps - 1)
                summ += integrand;
        }

        entropy_out[ti] = -rho * kb * h * summ;
    }

    atoms[py::str("entropy")] = entropy_out;
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