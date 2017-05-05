#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <cmath>


namespace ionmd { namespace constants {

constexpr double pi = 4*std::atan(1);
constexpr double amu = 1.660538782e-27;
constexpr double q_e = 1.602176487e-19;
constexpr double c = 2.99792458e8;

constexpr double OOFPEN = 8.9875518e+27;  // microns...
//constexpr double OOFPEN = 8.9875518e9;  // meters...

constexpr double HBAR = 1.0545716e-34;
constexpr double kB = 1.3806503e-23;
constexpr double g_elastic = 0.017;

} }  // namespace ionmd::constants

#endif
