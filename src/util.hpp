#ifndef UTIL_HPP
#define UTIL_HPP

#include <cmath>
#include <armadillo>

using std::sqrt;
using arma::vec;
using arma::dot;
using arma::mat;

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

#endif
