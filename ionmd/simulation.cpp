#include <iostream>
#include <cmath>
#include <vector>
#include <array>

#include "simulation.hpp"
#include "util.hpp"

namespace ionmd {

using arma::vec;
using arma::mat;


Simulation::Simulation() {
    status = SimStatus::IDLE;
}


Simulation::Simulation(SimParams p, Trap trap)
    : Simulation()
{
    this->p = std::make_shared<SimParams>(p);
    this->trap = std::make_shared<Trap>(trap);
}


Simulation::Simulation(SimParams p, Trap trap, std::vector<Ion> ions)
    : Simulation(p, trap)
{
    for (auto ion: ions) {
        this->ions.push_back(ion);
    }
    // BOOST_LOG_TRIVIAL(debug) << "Number of ions: " << this->ions.size();
}


mat Simulation::precompute_coulomb() {
    mat Flist(3, ions.size());

    #pragma omp parallel for
    for (unsigned int i = 0; i < ions.size(); i++) {
        // FIXME
        // Flist.col(i) = ion.coulomb(ions);
        i++;
    }
    return Flist;
}


void Simulation::set_params(SimParams new_params) {
    if (status != SimStatus::RUNNING) {
        p = std::make_shared<SimParams>(new_params);
    }
    else {
        // BOOST_LOG_TRIVIAL(error) << "Can't change parameters while simulation is running!";
    }
}



void Simulation::set_trap(Trap new_trap) {
    if (status != SimStatus::RUNNING) {
        trap = std::make_shared<Trap>(new_trap);
    }
    else {
        // BOOST_LOG_TRIVIAL(error) << "Can't set a new trap while simulation is running!";
    }
}


Ion Simulation::make_ion(const double &m, const double &Z,
                         const std::vector<double> &x0)
{
    return Ion(p, trap, m, Z, x0);
}


void Simulation::add_ion(const double &m, const double &z,
                         const std::vector<double> &x0)
{
    if (status != SimStatus::RUNNING) {
        ions.push_back(make_ion(m, z, x0));
    }
    else {
        // BOOST_LOG_TRIVIAL(error) << "Can't add an ion while simulation in progress!";
    }
}


void Simulation::set_ions(std::vector<Ion> ions) {
    if (status != SimStatus::RUNNING) {
        ions.clear();

        for (auto ion: ions) {
            ions.push_back(ion);
        }
    }
    else {
        // BOOST_LOG_TRIVIAL(error) << "Can't set new ions while simulation is running!";
    }
}


void Simulation::run() {
    // Initialize RNG
    // std::mersenne_twister_engine<double> rng;

    // Storage of pre-computed Coulomb force data
    mat coulomb_forces(3, ions.size());
    coulomb_forces.zeros();

    // FIXME: Initialization should happen before here
    // TODO: Initialize CCD
    // TODO: Initialize lasers

    // Stores every ion's position in one iteration
    auto current_positions = std::vector<vec>();
    for (unsigned int i = 0; i < p->buffer_size; i++) {
        current_positions.push_back(vec(3));
    }

    // Stores all ion positions
    auto trajectories = mat(p->num_steps, ions.size() * 3);

    // Run simulation
    // BOOST_LOG_TRIVIAL(info) << "Start simulation: " << timestamp_str() << "\n";
    status = SimStatus::RUNNING;

    auto t = double(0);
    for (unsigned int step = 0; step < p->num_steps; step++) {
        // Calculate Coulomb forces
        if (p->coulomb_enabled) {
            coulomb_forces = precompute_coulomb();
        }

        // Progress update
        // FIXME

        // Update each ion
        // TODO: update to use Boost.compute
        #pragma omp parallel for
        for (unsigned int i = 0; i < ions.size(); i++) {
            const auto x = ions[i].update(t, coulomb_forces);

            // Record trajectory position
            for (unsigned int j = 0; j < 3; j++) {
                trajectories(step, 3*i + j) = x[j];
            }
            current_positions[i] = x;

            // TODO: Check bounds
        }

        t += p->dt;
    }

    trajectories.save(p->filename, arma::raw_binary);
    status = SimStatus::FINISHED;
}

}  // namespace ionmd
