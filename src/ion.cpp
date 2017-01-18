#include <iostream>
#include <cmath>
#include <random>
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


Ion::Ion(std::shared_ptr<SimParams> params, std::shared_ptr<Trap> trap,
    std::shared_ptr<lasers_t> lasers,
    double m, double Z, vec x0)
{
    this->p = params;
    this->trap = trap;

    // FIXME
    // this->doppler_lasers = lasers;

    this->m = m;
    this->Z = Z;
    this->x = x0;
}


void Ion::print_position() {
    cout << "(" << std::scientific << x[0] << ", "
         << std::scientific << x[1] << ", "
         << std::scientific << x[2] << ")\n";
}


void Ion::update(double t, mat forces) {
    vec F, accel;

    x += v*p->dt + 0.5*a*pow(p->dt, 2);

    F = secular_force();
    if (p->micromotion_enabled)
	F += micromotion_force(t);
    if (p->coulomb_enabled)
	F += coulomb_force(forces);
    if (p->stochastic_enabled)
	F += stochastic_force();
    if (p->doppler_enabled)
	F += doppler_force();

    accel = F/m;
    v += 0.5*(a + accel)*p->dt;
    a = accel;
}


vec Ion::doppler_force() {
    vec F(3);
    F.zeros();

    // FIXME
    // double beta, F0;
    // beta = this->p->beta;
    // F0 = this->p->F0;

    // if(this->p->minimizing) {
    // 	beta = 1e-20; // unrealistically large damping for minimizing
    // 	if(this->doppler_coolable == 0) // don't use constant pressure term on non-lc'ed ions
    // 	    F0 = 0;
    // }

    // for (auto laser: this->doppler_lasers) {
	// use all lasers here
    //}

    // for (int i = 0; i < 3; i++) {
    //     F[i] = F0*p->khat[i] - beta*this->v[i];
    // }
    return F;
}


vec Ion::coulomb_force(mat forces) {
    // FIXME
    vec F(3);
    F.zeros();
    return F;
}


vec Ion::secular_force() {
    double A, B;
    vec F(3);

    //A = p->kappa*p->UEC/pow(p->z0,2);
    A = Z*pow(trap->V_rf, 2)/(m*pow(trap->omega_rf, 2)*pow(trap->r0, 4));
    B = trap->kappa*trap->U_ec/(2*pow(trap->z0,2));
    F[0] = -2*Z*(A-B)*x[0];
    F[1] = -2*Z*(A-B)*x[1];
    F[2] = -4*Z*B*x[2];
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
    // FIXME
    v = 0; // sqrt(2*kB*p->gamma_col*p->dt/this->m);
    F = this->m*v*hat/p->dt;
    return F;
}
