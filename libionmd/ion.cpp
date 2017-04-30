#include <iostream>
#include <cmath>
#include <random>

#include <armadillo>

#include <ionmd/ion.hpp>
#include <ionmd/util.hpp>
#include <ionmd/constants.hpp>

using std::cout;
using arma::vec;
using arma::mat;
using arma::normalise;

using namespace ionmd;

std::default_random_engine rng;
std::uniform_real_distribution<double> uniform(0, 1);


Ion::Ion(params_ptr params, trap_ptr trap, const double m, const double Z)
    : v(3), a(3), m(m), Z(Z)
{
    this->p = params;
    this->trap = trap;
    this->v.zeros();
    this->a.zeros();
}


Ion::Ion(params_ptr params, trap_ptr trap, const double m, const double Z,
         const vec x0)
    : Ion(params, trap, m, Z)
{
    this->x = x0;
    // BOOST_LOG_TRIVIAL(debug) << "m=" << m << ", Z=" << Z << ", x0="
    //                          << x0[0] << "," << x0[1] << "," << x0[2];
}


Ion::Ion(params_ptr params, trap_ptr trap, lasers_t lasers,
         const double m, const double Z, const vec x0)
    : Ion(params, trap, m, Z)
{
    // FIXME
    // this->doppler_lasers = lasers;
    this->x = x0;
}


const vec Ion::update(double t, mat forces)
{
    vec dx = v*p->dt + 0.5*a*pow(p->dt, 2);
    x += dx;

    vec F = secular_force()
        + micromotion_force(t)
        + coulomb_force(forces)
        + stochastic_force()
        + doppler_force();

    auto accel = F/m;
    v += 0.5*(a + accel)*p->dt;
    a = accel;

    return x;
}


vec Ion::doppler_force() {
    vec F = arma::zeros<vec>(3);

    if (p->doppler_enabled) {
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
    else {
        return F;
    }
}


vec Ion::coulomb_force(mat forces)
{
    vec F = arma::zeros<vec>(3);

    if (p->coulomb_enabled) {
        // FIXME
        return F;
    }
    else {
        return F;
    }
}


vec Ion::secular_force()
{
    vec F(3);

    //A = p->kappa*p->UEC/pow(p->z0,2);
    const auto mass = m*amu;
    const double A = Z*pow(trap->V_rf, 2)/(mass*pow(trap->omega_rf, 2)*pow(trap->r0, 4));
    const double B = trap->kappa*trap->U_ec/(2*pow(trap->z0,2));
    F[0] = -2*Z*(A-B)*x[0];
    F[1] = -2*Z*(A-B)*x[1];
    F[2] = -4*Z*B*x[2];
    return F;
}


vec Ion::micromotion_force(double t)
{
    vec F = arma::zeros<vec>(3);

    if (p->micromotion_enabled) {
        // FIXME
        return F;
    }
    else {
        return F;
    }
}


vec Ion::stochastic_force()
{
    vec F = arma::zeros<vec>(3);

    if (p->stochastic_enabled) {
        // FIXME
        return F;
    }
    else {
        return F;
    }

    // FIXME
    // double v;
    // vec direction = {uniform(rng), uniform(rng), uniform(rng)};
    // vec hat = normalise(direction);
    // v = 0; // sqrt(2*kB*p->gamma_col*p->dt/this->m);
    // F = this->m*v*hat/p->dt;
    // return F;
}
