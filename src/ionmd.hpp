#ifndef IONMD_HPP
#define IONMD_HPP

#include <armadillo>
#include "ion.hpp"
#include "trap.hpp"
#include "params.hpp"


namespace ionmd {

using arma::mat;

/**
 * Class that controls the overall simulation.
 */
class Simulation {
private:
    /// General simulation parameters.
    params_ptr p;

    /// Trap parameters.
    trap_ptr trap;

    /// All ions to simulate.
    std::vector<Ion> ions;

    /// Simulation is in progress.
    bool in_progress = false;

    /**
     * Precomputes all Coulomb interactions between ions that way they can be
     * applied all at once when advancing a time step.
     *
     * TODO: convert to OpenCL or similar.
     */
    mat precompute_coulomb();

public:
    Simulation();
    Simulation(SimParams p, Trap trap);
    Simulation(SimParams p, Trap trap, std::vector<Ion> ions);

    /**
     * Set or update simulation parameters. This method will only set parameters
     * when the simulation is not in progress.
     */
    void set_params(params_ptr p);

    /**
     * Set or update trap parameters. This method will only set parameters
     * when the simulation is not in progress.
     */
    void set_trap(trap_ptr trap);

    /**
     * Set ions. This method will only set parameters when the simulation is not
     * in progress.
     */
    void set_ions(std::vector<Ion> ions);

    /** Run the simulation. */
    void run();
};

}

#endif
