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
 * rho == 0 selects the "local" mode in which every atom uses its own
 * density n_i / (4/3 pi r_c,i^3), with r_c,i the neighbour cutoff of that
 * atom; otherwise the global density rho is used for all atoms.
 */

py::array_t<double> calculate_entropy(const nl_index& offsets,
    const nl_values& neighbordist,
    const nl_values& cutoff,
    double sigma,
    double rho,
    double rstart,
    double rstop,
    double h,
    double kb){

    const std::int64_t* off = offsets.data();
    const double* dist = neighbordist.data();
    const double* cutoff_vec = cutoff.data();
    const int nop = (int)(offsets.shape(0) - 1);

    const bool local = (rho == 0.0);

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

    py::array_t<double> entropy_array(nop);
    double* entropy_out = entropy_array.mutable_data();
    vector<double> raw_g(nsteps + 1);   // per-atom, reused

    for (int ti = 0; ti < nop; ti++) {

        int nn = (int)(off[ti + 1] - off[ti]);

        double rho_i = rho;
        if (local) {
            double rc = cutoff_vec[ti];
            rho_i = (rc > 0.0 && nn > 0)
                ? nn / (4.1887902047863905 * rc * rc * rc) : 0.0;
        }
        if (rho_i <= 0.0) {
            entropy_out[ti] = 0.0;
            continue;
        }

        /* ---- Accumulate raw Gaussian sum with windowing ---- */
        fill(raw_g.begin(), raw_g.end(), 0.0);

        for (int i = 0; i < nn; i++) {
            double rij = dist[off[ti] + i];
            int jmin = max(0,      (int)((rij - gauss_window - rstart) / h));
            int jmax = min(nsteps, (int)((rij + gauss_window - rstart) / h) + 1);
            for (int j = jmin; j <= jmax; j++) {
                double dr = r_grid[j] - rij;
                raw_g[j] += exp(dr * dr * inv_2sigma2);
            }
        }

        /* ---- Trapezoidal integration over [rstart, rstart + nsteps*h] ---- */
        double summ = 0.0;

        for (int j = 0; j <= nsteps; j++) {
            double r   = r_grid[j];
            double r2  = r * r;
            double frho = 4.0 * PI * rho_i * r2;

            /* g(r) = raw_g / (4π·ρ·r²·√(2πσ²)) */
            double g = (frho > 0.0) ? inv_fsigma * raw_g[j] / frho : 0.0;

            /* Integrand: (g·ln(g) - g + 1)·r²
               When g is negligibly small, limit → r². */
            double integrand;
            if (g > 1e-30)
                integrand = (g * log(g) - g + 1.0) * r2;
            else
                integrand = r2;

            /* Trapezoid weights: half at the end points. */
            if (j == 0 || j == nsteps)
                summ += 0.5 * integrand;
            else
                summ += integrand;
        }

        entropy_out[ti] = -rho_i * kb * h * summ;
    }

    return entropy_array;
}

py::array_t<double> calculate_average_entropy(const nl_index& offsets,
    const nl_index& neighbors,
    const nl_values& entropy){
    const std::int64_t* off = offsets.data();
    const std::int64_t* nb = neighbors.data();
    const double* ent = entropy.data();
    const py::ssize_t nop = offsets.shape(0) - 1;
    double entsum;

    py::array_t<double> avg_entropy(nop);
    double* out = avg_entropy.mutable_data();
    for (py::ssize_t ti=0; ti<nop; ti++){
        entsum = ent[ti];
        for (std::int64_t ci=off[ti]; ci<off[ti+1]; ci++){
            entsum += ent[nb[ci]];
        }
        out[ti] = entsum/(double(off[ti+1] - off[ti] + 1));
    }
    return avg_entropy;
}