#ifndef LASER_HPP
#define LASER_HPP

#include <vector>
#include <armadillo>


namespace ionmd {

using arma::vec;

class Laser {
public:
    double detuning;
    vec wave_vector;
};

typedef std::vector<Laser *> lasers_t;

}

#endif
