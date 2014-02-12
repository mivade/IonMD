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

//-------------//
//--CONSTANTS--//
//-------------//

const double pi = 4*std::atan(1);
const double amu = 1.660538782e-27;
const double q_e = 1.602176487e-19;
const double c = 2.99792458e8;
const double OOFPEN = 8.9875518e+09;
const double HBAR = 1.0545716e-34;
const double kB = 1.3806503e-23;
const double g_elastic = 0.017;	// elastic collision rate with 

typedef struct Ion {
    vec x(3), v(3), a(3);
    double m,			// mass
	Z;			// charge
    int index;
    int lc;			// 1 for laser cooled
} Ion;

//------------------------//
//--FUNCTION DEFINITIONS--//
//------------------------//

extern "C" {
    // Vector functions
    // inline void zeroVector(double *vec) {
    // 	vec[0] = vec[1] = vec[2] = 0.0;
    // }

    inline void copyVector(double *a, double *b) {
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
    }

    // inline double dot(double *a, double *b) {
    // 	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    // }

    inline void normalize(vec a) {
	double mag = sqrt(dot(a, a));
	a[0] /= mag;
	a[1] /= mag;
	a[2] /= mag;
    }

    // Utility functions
    void printIonStatistics(Params *p);
    void printParams(Params *p);
    void printPositions(Ion *ion);

    // Ion functions
    Ion *initIon(double *x0, double *v0, int index, Params *p);
    //void minimize(Ion **ions, Params *p);
    void simCCDPoint(Ion *ion, gsl_histogram2d **ccd, Params *p);
    void updateIon(Ion *ion, Ion **ions, double t, double *Fcoullist, Params *p);
    void FTrap(Ion *ion, double t, Params *p, vec F);
    void FLaser(Ion *ion, Params *p, vec F);
    void FCoulomb(Ion *ion, Ion **ions, Params *p, vec F);
    void FSecular(Ion *ion, double t, Params *p, vec F);
    void FStochastic(Ion *ion, Params *p, vec F);
    void allCoulomb(Ion **ions, Params *p, double *Flist);

    // Main simulation function
    int simulate(double *x0, double *v0, Params *p);
}

#endif
