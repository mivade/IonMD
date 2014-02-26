/*
  This file is part of IonMD.

  IonMD is free software: you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  IonMD is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
  License for more details.

  You should have received a copy of the GNU General Public License
  along with IonMD.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef PARAMS_HPP
#define PARAMS_HPP

//-------------------------//
//--SIMULATION PARAMETERS--//
//-------------------------//

// Parameters to pass from the simulation control scripts.
struct Params {
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
	Gamma,		// laser cooled ion's linewidth
	beta,		// cooling parameter (for constant cooling rate)
	F0;		// laser pressure force (~3.2E-20 kg m/s^2 for Ba+ at s = 1)

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
    int quit_after_minimizing;	// For debugging purposes.
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
	*com_fname,	// COM coordinate file name
	*fpos_fname,	// final positions file name
	*ccd_fname,	// file name prefix for generating simulated CCD images
	*temp_fname;	// file name for recording ensemble temperature
    int record_traj;	// record trajectories?
    double traj_start;	// time to start recording trajectory data
    int T_steps;	// number of steps for averaging velocities
};

#endif
