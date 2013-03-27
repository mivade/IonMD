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
    double *x, *v;	// position and velocity
    double m,		// mass
	Z;		// charge
    int index;
    int lc;		// 1 for laser cooled
} Ion;

//------------------------//
//--FUNCTION DEFINITIONS--//
//------------------------//

extern "C" {
    // Utility functions
    double *zeroVector(int length);
    double dot(double *a, double *b);
    void printIonStatistics(Params *p);
    void printParams(Params *p);
    void simCCDPoint(FILE *fout, Ion *ion);

    // Ion functions
    Ion *initIon(double *x0, double *v0, int index, Params *p);
    void updateIon(Ion *ion, Ion **ions, double t, double **Fcoul, Params *p);
    double *FTrap(Ion *ion, double t, Params *p);
    double *FLaser(Ion *ion, Params *p);
    double *FCoulomb(Ion *ion, Ion **ions, Params *p);
    double *FSecular(Ion *ion, double t, Params *p);
    double *FStochastic(Ion *ion, Params *p);
    double **allCoulomb(Ion **ions, Params *p);

    // Main simulation function
    int simulate(double *x0, double *v0, Params *p);
}

#endif
