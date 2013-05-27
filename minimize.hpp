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

// Old-style minimization by treating all ions as laser cooled
void minimize(Ion **ions, Params *p);

// Minimization main controlling function for nlopt
void minimize(double *x0, Ion **ions, Params *p);

// Writes initial positions to a file
void writeInitPos(Ion **ions, Params *p);
    
// Minimization function for nlopt
double minfunc(const vector<double> &x, vector<double> &grad, void *_data);

// Potential energy functions called by minfunc
double UTrap(int i, Ion **ions, Params *p);
double ULaser(int i, Ion **ions, Params *p);
double UCoulomb(int i, Ion **ions, Params *p);

#endif
