#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>
#include <memory>
#include "util.hpp"


namespace ionmd {

/**
 * Container structure for all parameters of a simulation.
 */
struct SimParams {
    /// Time step [&mu;s]
    double dt = 0.01;

    /// Total number of time steps.
    unsigned int num_steps = 20000;

    /// Verbosity of output. Higher numbers increases verbosity.
    unsigned int verbosity = 0;

    /// Enable micromotion calculations (requires smaller time steps)
    bool micromotion_enabled = false;

    /// Enable Coulomb repulsion calculation
    bool coulomb_enabled = true;

    /// Enable stochastic force calculation
    bool stochastic_enabled = false;

    /// Enable Doppler cooling simulation
    bool doppler_enabled = false;

    /// Filename for writing trajectory data to.
    std::string filename = "out.bin";

    /// How many points in time to store before writing to disk.
    size_t buffer_size = 10000;
};

typedef std::shared_ptr<SimParams> params_ptr;

}

#endif
