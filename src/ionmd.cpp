#include <iostream>
#include <cmath>
#include <vector>
#include <thread>

#include <omp.h>

#include <armadillo>

#include "ion.hpp"
#include "trap.hpp"

using arma::vec;
using arma::mat;
using std::cout;
using std::cerr;
using std::endl;
using namespace ionmd;


/**
   Precomputes all Coulomb interactions (that way they can be applied
   all at once instead of having one ion updated before the Coulomb
   interaction is computed for the rest).

   TODO: convert to OpenCL
*/
mat precompute_coulomb(std::vector<Ion> ions) {
    int i = 0;
    mat Flist(3, ions.size());
    #pragma omp parallel for
    for (auto ion: ions) {
	Flist.col(i) = ion.coulomb(ions);
	i++;
    }
    return Flist;
}


/** Main entry point to run simulations. */
int simulate(SimParams params, Trap trap) {
    auto params_p = std::make_shared<SimParams>(params);
    auto trap_p = std::make_shared<Trap>(trap);

    // Initialize RNG
    // std::mersenne_twister_engine<double> rng;

    // Set number of threads for multiprocessing
    // TOOD: convert to OpenCL and use GPU when available
    const int num_threads = std::thread::hardware_concurrency();
    omp_set_num_threads(num_threads);

    // Storage of pre-computed Coulomb force data
    mat coulomb_forces(3, params_p->num_ions);
    coulomb_forces.zeros();

    // TODO: Initialize CCD

    // TODO: Initialize lasers

    // Initialize ions
    std::vector<Ion> ions;
    const vec x0 = arma::zeros<vec>(3);

    for (unsigned int i = 0; i < params_p->num_ions; i++) {
	// TODO: place on grid
	// TODO: figure out how to specify mass and charge in params
        ions.push_back(Ion(params_p, trap_p, 40, 1, x0));
    }

    // TODO: recording initialization

    // Run simulation
    int index = 0;
    int t_10 = (int)(params_p->t_max/params_p->dt) / 10;
    cout << "Simulating..." << endl;

    for(double t = 0; t < params_p->t_max; t += params_p->dt) {
        // Calculate Coulomb forces
        if (params_p->coulomb_enabled)
            coulomb_forces = precompute_coulomb(ions);

	// Progress update
	if (index % t_10 == 0 && params_p->verbosity != 0) {
	    cout << int(10*index/t_10) << "% complete; "
		 << "t = " << t << "\n";
	}

        // Update each ion
	for (auto ion: ions) {
            // TODO: Record data

	    ion.update(t, coulomb_forces);

	    // TODO: Check bounds
        }

        index++;
    }

    return 0;
}
