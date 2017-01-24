#include <iostream>
#include <cmath>
#include <vector>
#include <thread>

#include <omp.h>

#include "ionmd.hpp"

namespace ionmd {

using arma::vec;
using arma::mat;
using std::cout;
using std::cerr;
using std::endl;


Simulation::Simulation()
{
}


Simulation::Simulation(SimParams p, Trap trap) {
    this->p = std::make_shared<SimParams>(p);
    this->trap = std::make_shared<Trap>(trap);
}


Simulation::Simulation(SimParams p, Trap trap, std::vector<Ion> ions)
    : Simulation(p, trap)
{
    this->ions = ions;
}


mat Simulation::precompute_coulomb() {
    int i = 0;
    mat Flist(3, ions.size());
    #pragma omp parallel for
    for (auto ion: ions) {
	Flist.col(i) = ion.coulomb(ions);
	i++;
    }
    return Flist;
}


void Simulation::run() {
    // Initialize RNG
    // std::mersenne_twister_engine<double> rng;

    // Set number of threads for multiprocessing
    // TOOD: convert to OpenCL and use GPU when available
    const int num_threads = std::thread::hardware_concurrency();
    omp_set_num_threads(num_threads);

    // Storage of pre-computed Coulomb force data
    mat coulomb_forces(3, p->num_ions);
    coulomb_forces.zeros();

    // TODO: Initialize CCD

    // TODO: Initialize lasers

    // Initialize ions
    std::vector<Ion> ions;
    const vec x0 = arma::zeros<vec>(3);

    for (unsigned int i = 0; i < p->num_ions; i++) {
	// TODO: place on grid
	// TODO: figure out how to specify mass and charge in params
        ions.push_back(Ion(p, trap, 40, 1, x0));
    }

    // TODO: recording initialization

    // Run simulation
    int index = 0;
    int t_10 = (int)(p->t_max/p->dt) / 10;
    cout << "Simulating..." << endl;

    for(double t = 0; t < p->t_max; t += p->dt) {
        // Calculate Coulomb forces
        if (p->coulomb_enabled)
            coulomb_forces = precompute_coulomb();

	// Progress update
	if (index % t_10 == 0 && p->verbosity != 0) {
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
}

}  // namespace ionmd
