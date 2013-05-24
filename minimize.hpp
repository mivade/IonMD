#ifndef MINIMIZE_HPP
#define MINIMIZE_HPP

#include <vector>
#include "ionmd.hpp"

using std::vector;

// Struct for passing data to minimization function
typedef struct MinData {
    int i, fcalls;
    Params *p;
    Ion **ions;
} MinData;

extern "C" {
    // Minimization main controlling function
    void minimize(double *x0, Ion **ions, Params *p);
    
    // Minimization function for nlopt
    double minfunc(const vector<double> &x, vector<double> &grad, void *_data);

    // Constraint function
    // need to make this vector valued
    // the idea is, for x,y, return abs(x-r0); for z: return abs(z-z0/2.0)
    double cfunc(const vector<double> &x, vector<double> &grad, void *_data);

    // Potential energy functions called by minfunc
    double UTrap(int i, Ion **ions, Params *p);
    double ULaser(int i, Ion **ions, Params *p);
    double UCoulomb(int i, Ion **ions, Params *p);
}

#endif
