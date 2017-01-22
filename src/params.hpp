#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>


namespace ionmd {

/**
 * Container structure for all parameters of a simulation. Note that all units
 * are SI unless otherwise specified.
 */
struct SimParams {
    /// Total number of ions
    unsigned int num_ions;

    /// Time step
    double dt = 50e-9;

    /// Max time
    double t_max = 1e-3;

    /// Enable micromotion calculations (requires smaller time steps)
    bool micromotion_enabled = false;

    /// Enable Coulomb repulsion calculation
    bool coulomb_enabled = true;

    /// Enable stochastic force calculation
    bool stochastic_enabled = false;

    /// Enable Doppler cooling simulation
    bool doppler_enabled = false;

    // Output
    std::string output_filename;
    bool record_trajectories;
};

}

#endif
