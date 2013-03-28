#ifndef IONMD_HPP
#define IONMD_HPP

#include <cmath>
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
double xhat[] = {-sqrt(2)/2, sqrt(2)/2, 0};
double yhat[] = {0,0,1};

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
    // Utility functions
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

    void printIonStatistics(Params *p);
    void printParams(Params *p);
    void simCCDPoint(FILE *fout, Ion *ion);

    // Ion functions
    Ion *initIon(double *x0, double *v0, int index, Params *p);
    void updateIon(Ion *ion, Ion **ions, double t, double *Fcoullist, Params *p);
    void FTrap(Ion *ion, double t, Params *p, double *F);
    void FLaser(Ion *ion, Params *p, double *F);
    void FCoulomb(Ion *ion, Ion **ions, Params *p, double *F);
    void FSecular(Ion *ion, double t, Params *p, double *F);
    void FStochastic(Ion *ion, Params *p, double *F);
    void allCoulomb(Ion **ions, Params *p, double *Flist);

    // Main simulation function
    int simulate(double *x0, double *v0, Params *p);
}

#endif
