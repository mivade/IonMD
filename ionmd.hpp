#ifndef IONMD_HPP
#define IONMD_HPP

#include <cmath>
#include <cstdlib>
#include <gsl/gsl_histogram2d.h>
#include "params.hpp"

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
    double x[3], v[3], a[3];	// position, velocity, and acceleration
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
    inline void zeroVector(double *vec) {
	vec[0] = vec[1] = vec[2] = 0.0;
    }

    inline void copyVector(double *a, double *b) {
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
    }

    inline double dot(double *a, double *b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    }

    inline void normalize(double *a) {
	double mag = sqrt(dot(a,a));
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
    void FTrap(Ion *ion, double t, Params *p, double *F);
    void FLaser(Ion *ion, Params *p, double *F);
    void FCoulomb(Ion *ion, Ion **ions, Params *p, double *F);
    void FSecular(Ion *ion, double t, Params *p, double *F);
    void FStochastic(Ion *ion, Params *p, double *F);
    void allCoulomb(Ion **ions, Params *p, double *Flist);
    void swapIons(Ion **ions, int i, int j);	// swap ions i and j

    // Main simulation function
    int simulate(double *x0, double *v0, Params *p);
}

#endif
