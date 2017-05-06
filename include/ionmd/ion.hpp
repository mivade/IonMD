#ifndef ION_HPP
#define ION_HPP

#include <vector>
#include <memory>
#include <armadillo>
#include <ionmd/laser.hpp>
#include <ionmd/trap.hpp>
#include <ionmd/params.hpp>

namespace ionmd {

using arma::vec;
using arma::mat;


class Ion {
private:
    /// Common simulation parameters
    params_ptr p;

    /// Trap pointer for computing trap-related forces.
    trap_ptr trap;

    /// Doppler cooling lasers affecting this ion.
    lasers_ptr lasers;

    /// Charge in Coulombs
    const double charge;

    /// Compute the force from Doppler cooling lasers.
    vec doppler_force();

    /// Compute the secular motion (ponderomotive) force of the trap.
    vec secular_force();

    /**
     * Compute micromotion (rf) force due to the trap.
     * @param t
     */
    vec micromotion_force(double t);

    /**
     * Sum all forces due to other ions.
     * @param forces
     * @param index
     */
    vec coulomb_force(const mat &forces, const unsigned int &index);

    /// Compute stochastic forces (background gas collisions, etc.).
    vec stochastic_force();

public:
    // FIXME: consider making v, a, m private

    vec x;  /// Ion position
    vec v;  /// Ion velocity
    vec a;  /// Ion acceleration

    const double m;  /// Ion mass
    const double Z;  /// Ion charge in units of [e]

    /**
     * @param params
     * @param trap
     * @param m Ion mass
     * @param Z Ion charge
     */
    Ion(params_ptr params, trap_ptr trap, double m, double Z);

    /**
     * @param params
     * @param trap
     * @param m Ion mass
     * @param Z Ion charge
     * @param x0 Initial position vector
     */
    Ion(params_ptr params, trap_ptr trap, double m, double Z, vec x0);

    /**
     * @param params
     * @param trap
     * @param lasers Shared pointer to Doppler cooling lasers that affect this ion
     * @param m Ion mass
     * @param Z Ion charge in units of e
     * @param x0 Initial position
     */
    Ion(params_ptr params, trap_ptr trap, lasers_ptr lasers,
        double m, double Z, vec x0);

    /**
     * Apply a single time step of integration.
     * @param t Current time
     * @param forces Pre-computed Coulomb forces due to all other ions
     * @param index
     * @returns The new ion position
     */
    const vec update(const double &t, const mat &forces,
                     const unsigned int &index);

    /**
     * Pre-compute the Coluomb forces due to all other ions in the trap.
     * @param ions Vector of all ions in the trap.
     */
    const vec coulomb(const std::vector<Ion> &ions);
};

}  // namespace ionmd

#endif
