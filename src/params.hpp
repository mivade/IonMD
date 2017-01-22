#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>


namespace ionmd {

struct SimParams {
  // Total number of ions
  unsigned int num_ions;

  // Simulation settings
  double dt, t_max;  // time step, max time
  unsigned int verbosity;

  // Enable/disable some types of forces
  bool micromotion_enabled;
  bool coulomb_enabled;
  bool stochastic_enabled;
  bool doppler_enabled;

  // Output
  std::string output_filename;
  bool record_trajectories;
};

}

#endif
