#ifndef LASER_HPP
#define LASER_HPP

#include <armadillo>

using arma::vec;

namespace ionmd {
    class Laser {
    public:
	double detuning;
	vec wave_vector;
    };
}

#endif
