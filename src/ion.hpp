#ifndef ION_HPP
#define ION_HPP

#include <vector>
#include <memory>
#include <armadillo>
#include "ionmd.hpp"
#include "laser.hpp"
#include "trap.hpp"
#include "params.hpp"


namespace ionmd {

using arma::vec;
using arma::mat;

class Ion {
private:
  /// Unique ID
  unsigned int id;

  /// Common simulation parameters
  std::shared_ptr<SimParams> p;

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
   */
  vec coulomb_force(mat forces);

  /// Compute stochastic forces (background gas collisions, etc.).
  vec stochastic_force();

public:
  vec x;  /// Ion position
  vec v;  /// Ion velocity
  vec a;  /// Ion acceleration

  double m;  /// Ion mass.
  double Z;  /// Ion charg in units of e.

  /// Trap pointer for computing trap-related forces.
  std::shared_ptr<Trap> trap;

  /**
   */
  Ion();

  /**
   * @param p
   * @param trap
   * @param lasers Doppler cooling lasers that affect this ion
   * @param m Ion mass
   * @param Z Ion charge in units of e
   * @param x0 Initial position
   */
  Ion(std::shared_ptr<SimParams> params, std::shared_ptr<Trap> trap,
      const lasers_t lasers,
      const double m, const double Z, const vec x0);

  /**
   * Apply a single time step of integration.
   * @param t Current time
   * @param forces Pre-computed Coulomb forces due to all other ions
   */
  void update(double t, mat forces);

  /**
   * Pre-compute the Coluomb forces due to all other ions in the trap.
   * @param ions Vector of all ions in the trap.
   */
  vec coulomb(const std::vector<Ion> ions);

  /// Convenience function to pretty-print the ion position.
  void print_position();
};

}  // namespace ionmd

#endif
