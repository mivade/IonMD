#ifndef MINIMIZE_HPP
#define MINIMIZE_HPP

#include <vector>

using std::vector;

extern "C" {
    // Minimization function for nlopt
    void minfunc(const vector<double> &x, vector<double> &grad, void *data);

    // Potential energy functions called by minfunc
    double UTrap(int i, Ion **ions, Params *p);
    double ULaser(int i, Ion **ions, Params *p);
    double UCoulomb(int i, Ion **ions, Params *p);
}

#endif
