#ifndef PARAMS_HPP
#define PARAMS_HPP

//-------------------------//
//--SIMULATION PARAMETERS--//
//-------------------------//

// Parameters to pass from the simulation control scripts
typedef struct Params {
    // Ion parameters
    int N;		// total number of ions
    double *m, *Z;	// masses and charges
    int *lc;		// laser cooled or not

    // Laser parameters
    double *khat;	// direction
    double lmbda,	// wavelength
	r_l,		// beam radius
	delta,		// detuning
	s0,		// saturation parameter
	Gamma;		// laser cooled ion's linewidth

    // Trap parameters
    double r0, z0, kappa,	// radius, length, geometric paremter
	Omega,			// RF frequency
	V, U, UEC,		// RF amplitude, offset, end cap voltages
	Vsec, w;		// secular amplitude, frequency

    // Background gas parameters
    double gamma_col;	// collision rate

    // Simulation settings
    double dt, t_max;	// time step and max time
    int t_steps,	// number of time steps (easy to give with Python)
	use_rfmm,	// include RF micromotion
	use_coulomb,	// include the Coulomb interaction
	use_laser,	// include laser cooling
	use_secular,	// include secular excitation
	use_stochastic,	// include stochastic processes
	num_threads,	// number of threads to use for multiprocessing
	quiet;		// 1 to turn off progress reports, etc.

    // Data recording
    char *traj_fname,	// trajectory file name
	*fpos_fname,	// final positions file name
	*ccd_fname;	// file name for generating simulated CCD images
    int record_traj;	// record trajectories?
    double traj_start;	// time to start recording trajectory data
} Params;

#endif
