#include <iostream>
#include <cmath>
#include <random>
#include <boost/format.hpp>
#include <armadillo>
#include "ion.hpp"
#include "util.hpp"
#include "constants.hpp"

using std::cout;
using arma::vec;
using arma::mat;
using namespace ionmd;

std::default_random_engine rng;
std::uniform_real_distribution<double> uniform(0, 1);

Ion::Ion(SimParams *p, Trap *trap, lasers_t lasers, double m, double Z, vec x0) {
    this->p = p;
    this->trap = trap;
    this->doppler_lasers = lasers;
    this->m = m;
    this->Z = Z;
    this->x = x0;
}

void Ion::print_position() {
    cout << boost::format("%e %e %e\n") % this->x[0] % this->x[1] % this->x[2];
}

void Ion::update(double t, mat forces) {
    vec F, a;

    this->x += this->v*this->p->dt + 0.5*this->a*pow(this->p->dt, 2);

    F = this->secular_force();
    if (this->p->micromotion_enabled)
	F += this->micromotion_force(t);
    if (this->p->coulomb_enabled)
	F += this->coulomb_force(forces);
    if (this->p->stochastic_enabled)
	F += this->stochastic_force();
    if (this->p->doppler_enabled)
	F += this->doppler_force();

    a = F/this->m;
    this->v += 0.5*(this->a + a)*this->p->dt;
    this->a = a;
}

vec Ion::doppler_force() {
    vec F(3);
    F.zeros();
    // double beta, F0;
    // beta = this->p->beta;
    // F0 = this->p->F0;

    // if(this->p->minimizing) {
    // 	beta = 1e-20; // unrealistically large damping for minimizing
    // 	if(this->doppler_coolable == 0) // don't use constant pressure term on non-lc'ed ions
    // 	    F0 = 0;
    // }

    for (auto laser: this->doppler_lasers) {
	// use all lasers here
    }

    // for (int i = 0; i < 3; i++) {
    //     F[i] = F0*p->khat[i] - beta*this->v[i];
    // }
    return F;
}

vec Ion::coulomb_force(mat forces) {
}

vec Ion::secular_force() {
    double A, B;
    vec F(3);

    //A = p->kappa*p->UEC/pow(p->z0,2);
    A = this->Z*pow(this->trap->V_rf, 2)/(this->m*pow(this->trap->omega_rf, 2)*pow(this->trap->r0, 4));
    B = this->trap->kappa*this->trap->U_ec/(2*pow(this->trap->z0,2));
    F[0] = -2*this->Z*(A-B)*this->x[0];
    F[1] = -2*this->Z*(A-B)*this->x[1];
    F[2] = -4*this->Z*B*this->x[2];
    return F;
}

vec Ion::micromotion_force(double t) {
    // TODO
    vec F(3);
    F.zeros();
    return F;
}

vec Ion::stochastic_force() {
    double v;
    vec F(3);
    F.zeros();
    vec hat = {uniform(rng), uniform(rng), uniform(rng)};
    normalize(hat);
    v = sqrt(2*kB*p->gamma_col*p->dt/this->m);
    F = this->m*v*hat/p->dt;
    return F;
}
