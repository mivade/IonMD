#ifndef PARAMS_HPP
#define PARAMS_HPP

//-------------------------//
//--SIMULATION PARAMETERS--//
//-------------------------//

// Parameters to pass from the simulation control scripts
typedef struct Params {
    // Ion parameters
    int N,		// total number of ions
	N_masses;	// number of unique masses
    double *m, *Z,	// masses and charges
	*masses;	// list of unique masses
    int *lc;		// laser cooled or not

    // Laser parameters
    double *khat;	// direction
    double lmbda,	// wavelength
	r_l,		// beam radius
	delta,		// detuning
	s,		// saturation parameter
	Gamma;		// laser cooled ion's linewidth

    // Trap parameters
    double r0, z0, kappa,	// radius, length, geometric paremter
	Omega,			// RF frequency
	V, U, UEC,		// RF amplitude, offset, end cap voltages
	Vsec, w;		// secular amplitude, frequency

    // Background gas parameters
    double gamma_col,	// collision rate
	m_gas,		// background gas mass
	T_gas;		// background gas temperature

    // CCD settings
    int sim_ccd;	// set to 1 to simulate CCD
    double ccd_bins,	// number of bins ("pixels") of the simulated CCD
	ccd_extent;	// physical extent of CCD in um

    // Simulation settings
    double dt, t_max,	// time step and max time
	min_time,	// amount of time to run minimization routine for
	abort_bounds;	// radial postion at which ions are considered
			// out of bounds
    int minimizing;	// internal variable; 1 when performing
			// minimization routine
    int t_steps,	// number of time steps (easy to give with Python)
	use_rfmm,	// include RF micromotion
	use_coulomb,	// include the Coulomb interaction
	use_laser,	// include laser cooling
	use_secular,	// include secular excitation
	use_stochastic,	// include stochastic processes
	use_abort,	// abort if ions are out of bounds
	num_threads,	// number of threads to use for multiprocessing
	quiet;		// 1 to turn off progress reports, etc.

    // Data recording
    char *traj_fname,	// trajectory file name
	*fpos_fname,	// final positions file name
	*ccd_fname,	// file name prefix for generating simulated CCD images
	*temp_fname;	// file name for recording ensemble temperature
    int record_traj;	// record trajectories?
    double traj_start;	// time to start recording trajectory data
    int T_steps;	// number of steps for averaging velocities
} Params;

#endif
