/*
  This file is part of IonMD.

  IonMD is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  IonMD is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU General Public License
  along with IonMD.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef IONMD_HPP
#define IONMD_HPP

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_histogram2d.h>
#include <armadillo>
#include "params.hpp"

using arma::vec;
using arma::dot;
using arma::mat;

// Constants
// =========

const double pi = 4*std::atan(1);
const double amu = 1.660538782e-27;
const double q_e = 1.602176487e-19;
const double c = 2.99792458e8;
const double OOFPEN = 8.9875518e+09;
const double HBAR = 1.0545716e-34;
const double kB = 1.3806503e-23;

// Ion struct
// ==========

struct Ion {
    vec x, v, a;
    double m,			// mass
	Z;			// charge
    int index;
    int lc;			// 1 for laser cooled
};

// Functions
// =========

// Vector functions
// ----------------

inline void copy_vector(vec a, double *b) {
    a[0] = b[0];
    a[1] = b[1];
    a[2] = b[2];
}

inline void normalize(vec a) {
    double mag = sqrt(dot(a, a));
    a[0] /= mag;
    a[1] /= mag;
    a[2] /= mag;
}

inline double dot(double *a, vec b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double dot(vec a, double *b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Utility functions
// -----------------

void write_xyz_file(char *fname, Ion **ions, Params *p);

extern "C" {
    // Parameter printing functions
    // ----------------------------

    void print_ion_statistics(Params *p);
    void print_params(Params *p);
    void print_positions(Ion *ion);

    // Ion functions
    // -------------

    Ion *initIon(double *x0, double *v0, int index, Params *p);
    void simCCDPoint(Ion *ion, gsl_histogram2d **ccd, Params *p);
    void updateIon(Ion *ion, Ion **ions, double t, mat Fcoullist, Params *p);
    vec FTrap(Ion *ion, double t, Params *p);
    vec FLaser(Ion *ion, Params *p);
    vec FCoulomb(Ion *ion, Ion **ions, Params *p);
    vec FSecular(Ion *ion, double t, Params *p);
    vec FStochastic(Ion *ion, Params *p);
    mat allCoulomb(Ion **ions, Params *p);

    // Main simulation function
    // ------------------------

    int simulate(double *x0, double *v0, Params *p);
}

#endif
