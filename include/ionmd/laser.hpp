#ifndef LASER_HPP
#define LASER_HPP

#include <vector>
#include <memory>
#include <armadillo>
#include <ionmd/constants.hpp>


namespace ionmd {

using arma::vec;

/**
 * Doppler cooling laser class. This currently uses a damping model of Doppler
 * cooling.
 *
 * Note that it is up to the user to keep track of which lasers should affect
 * which ions (for simulations involving multiple ion species). This class is
 * completely unaware of which ions it should affect.
 *
 */
class Laser
{
public:
    /// Angular frequency detuning from transition (not used in damping model)
    double detuning;

    /// The laser's k-vector
    vec wave_vector;

    /// Damping coefficient
    double beta;

    /// Constant force term
    double F0;

    /**
     * @param beta Damping coefficient
     * @param F0 Constant damping term
     * @param khat Unit vector indicating pointing direction
     */
    Laser(double beta, double F0, vec khat)
    {
        this->beta = beta;
        this->F0 = F0;
        this->wave_vector = khat / arma::norm(khat);  // Ensure khat is actually normalized
    }

    /**
     * @param detuning Detuning from transition (angular frequency)
     * @param beta Damping coefficient
     * @param F0 Constant damping term
     * @param wavelength Laser wavelength
     * @param khat Unit vector indicating pointing direction
     */
    Laser(double detuning, double beta, double F0, double wavelength, vec khat)
        : Laser(beta, F0, khat)
    {
        this->detuning = detuning;
        khat = khat / arma::norm(khat);  // Ensure khat is actually normalized
        this->wave_vector = khat * 2*constants::pi / wavelength;
    }
};

typedef std::shared_ptr<Laser> laser_ptr;
typedef std::vector<laser_ptr> lasers_ptr;

}

#endif
