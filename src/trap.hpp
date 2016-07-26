#ifndef TRAP_HPP
#define TRAP_HPP

namespace ionmd {
    struct Trap {
	// Geometric parameters
	double r0, z0, kappa;

	// rf frequency
	double omega_rf;

	// rf amplitude, offset, and endcap voltage
	double V_rf, U_dc, U_ec;
    };
}

#endif
