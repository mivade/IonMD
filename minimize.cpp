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

#include <cstdio>
#include <iostream>
#include <cstdlib>
//#include <nlopt.hpp>
#include "minimize.hpp"
#include "ionmd.hpp"

using std::copy;
using std::cout;
using std::cerr;
using std::endl;

void minimize(Ion **ions, Params *p) {
    double t = 0;
    double *Fclist = new double[p->N*3];

    // Store settings
    int use_rfmm = p->use_rfmm,
	use_coulomb = p->use_coulomb,
	use_laser = p->use_laser,
	use_stochastic = p->use_stochastic,
	i, j;

    // Turn off everything but Coulomb, trap, and laser
    p->use_rfmm = 0;
    p->use_coulomb = 1;
    p->use_laser = 1;
    p->use_stochastic = 0;

    // Begin minimization
    p->minimizing = 1;
    while(true) {	// TODO: better end decision
	allCoulomb(ions, p, Fclist);
	for(i=0; i<p->N; i++) {
	    updateIon(ions[i], ions, t, Fclist, p);
	}
	if(t >= p->min_time)
	    break;
	t += p->dt;
    }
    for(i=0; i<p->N; i++) {
	for(j=0; j<3; j++) {
	    ions[i]->v[j] = 0;
	    ions[i]->a[j] = 0;
	}
    }
    p->minimizing = 0;

    // Restore settings
    p->use_rfmm = use_rfmm;
    p->use_coulomb = use_coulomb;
    p->use_laser = use_laser;
    p->use_stochastic = use_stochastic;
    delete[] Fclist;

    // Write initial positions
    writeInitPos(ions, p);
}

/*void minimize(double *x0, Ion **ions, Params *p) {
    int i;
    double Umin;

    // Vectors for passing stuff to NLopt
    vector<double> ubounds(3*p->N), lbounds(3*p->N), x(3*p->N);

    // Data for passing to minimization routine
    MinData *min_data = new MinData;
    min_data->i = 0;
    min_data->fcalls = 0;
    min_data->p = p;
    min_data->ions = ions;

    // Setup optimization
    nlopt::opt opt(nlopt::LN_PRAXIS, 3*p->N);
    for(i=0; i<3*p->N; i+=3) {
	ubounds[i] = p->r0/2.;
	ubounds[i+1] = ubounds[i];
	ubounds[i+2] = p->z0;
	lbounds[i] = -ubounds[i];
	lbounds[i+1] = -ubounds[i+1];
	lbounds[i+2] = -ubounds[i+2];
	x[i] = x[i+1] = x[i+2] = 0;
    }
    opt.set_upper_bounds(ubounds);
    opt.set_lower_bounds(lbounds);
    opt.set_min_objective(minfunc, min_data);
    opt.set_initial_step(1e-6);
    opt.set_xtol_rel(2e-8);
    opt.set_ftol_rel(0);
    //copy(x0, x0+3*p->N, x.begin());
    opt.optimize(x, Umin);
    printf("Finished minimization in %d function calls.\n", min_data->fcalls);
    cout << "Minimum energy: " << Umin << endl;
    delete min_data;
    writeInitPos(ions, p);
    }*/

void writeInitPos(Ion **ions, Params *p) {
    FILE *init_pos_file = fopen("ipos.xyz", "w");
    fprintf(init_pos_file, "%d\nInitial positions in microns\n", p->N);
    for(int i=0; i<p->N; i++) {
	//printf("%e %e %e\n", ions[i]->x[0], ions[i]->x[1], ions[i]->x[2]);
	fprintf(init_pos_file, "%d %e %e %e\n", int(ions[i]->m/amu),
		ions[i]->x[0]/1e-6, ions[i]->x[1]/1e-6, ions[i]->x[2]/1e-6);
    }
    fclose(init_pos_file);
}

double minfunc(const vector<double> &x, vector<double> &grad, void *_data) {
    int i,j;
    double U = 0;
    MinData *data = (MinData*)_data;
    #pragma omp parallel for
    for(i=0; i<data->p->N; i++) {
	for(j=0; j<3; j++)
	    data->ions[i]->x[j] = x[i+j];
    }
    #pragma omp parallel for
    for(i=0; i<data->p->N; i++) {
	U += UTrap(i, data->ions, data->p) + \
	    ULaser(i, data->ions, data->p) + \
	    UCoulomb(i, data->ions, data->p);
    }
    data->fcalls++;
    return U;
}

double UTrap(int i, Ion **ions, Params *p) {
    double C, D, U;
    C = ions[i]->Z*pow(p->V,2)/(ions[i]->m*pow(p->Omega,2)*pow(p->r0,4));
    D = p->kappa*p->UEC/(2*pow(p->z0,2));
    U = ions[i]->Z*((C-D)*pow(ions[i]->x[0],2) + (C-D)*pow(ions[i]->x[1],2) + \
		    2*D*pow(ions[i]->x[2],2));
    return U;
}

double ULaser(int i, Ion **ions, Params *p) {
    double F0, k, s, s0, Gamma, beta, delta, U;
    k = 2*pi/p->lmbda;
    Gamma = p->Gamma;
    s = p->s;
    delta = p->delta + k*dot(p->khat, ions[i]->v);
    s0 = s*(1 + pow(2*delta/Gamma,2));
    beta = -HBAR*k*k*4*s0*delta/Gamma/pow(1+s0+pow(2*delta/Gamma,2),2);
    if(ions[i]->lc)
	F0 = abs(HBAR*k*s*Gamma/(2*(1+s)));
    else
	F0 = 0;
    // abs?
    U = -F0*dot(ions[i]->x, p->khat) + 0.5*beta*dot(ions[i]->x, ions[i]->x);
    return U;
}

double UCoulomb(int i, Ion **ions, Params *p) {
    double U = 0;
    double r[3];
    #pragma omp parallel for
    for(int j=0; j<p->N; j++) {
	if(i == j)
	    continue;
	for(int k=0; k<3; k++)
	    r[k] = ions[i]->x[k] - ions[j]->x[k];
	U += OOFPEN*ions[i]->Z/sqrt(dot(r,r));
    }
    return U;
}
