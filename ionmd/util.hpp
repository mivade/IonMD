#ifndef UTIL_HPP
#define UTIL_HPP

#include <ctime>
#include <string>
#include <cmath>
#include <armadillo>

using std::sqrt;
using arma::vec;
using arma::dot;
using arma::mat;

namespace ionmd {

inline auto timestamp_str() -> std::string
{
    const auto now = std::time(nullptr);
    char buff[256];
    std::strftime(buff, sizeof(buff), "%FT%T%z", std::localtime(&now));
    return std::string(buff);
}

}  // namespace ionmd

//FIXME: figure out if below are still required

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
