#include "minimize.hpp"
#include "ionmd.hpp"

void minfunc(const vector<double> &x, vector<double> &grad, void *data) {
    double U = 0;
    for(int i=0; i<p->N; i++) {
	U += UTrap(i, x) + ULaser(i, x) + UCoulomb(i, x, p);
    }
    return U;
}

void UTrap(int i, Ion **ions, Params *p) {
    double C, D, U;
    C = ions[i]->Z*pow(p->V,2)/(ions[i]->m*pow(p->Omega,2)*pow(p->r0,4));
    D = p->kappa*p->UEC/(2*pow(p->z0,2));
    U = ions[i]->Z*((C-D)*pow(ions[i]->x[0],2) + (C-D)*pow(ions[i]->x[1],2) + \
		    2*D*pow(ions[i]->x[2],2));
    return U;
}

void ULaser(int i, Ion **ions, Params *p) {
    double F0, k, s, s0, Gamma, beta, delta, U;
    k = 2*pi/p->lmbda;
    Gamma = p->Gamma;
    s = p->s;
    delta = p->delta + k*dot(p->khat, ion->v);
    s0 = s*(1 + pow(2*delta/Gamma,2));
    beta = -HBAR*pow(k,2)*4*s0*delta/Gamma/pow(1+s0+pow(2*delta/Gamma,2),2);
    if(ions[i]->lc != 0)
	F0 = abs(HBAR*k*s*Gamma/(2*(1+s)));
    else
	F0 = 0;
    U = -F0*dot(ions[i]->x, p->khat) + 0.5*beta*dot(ions[i]->x, ions[i]->x);
    return U;
}

void UCoulomb(int i, Ion **ions, Params *p) {
    double U = 0;
    double r[3];
    #pragma parallel for
    for(int j=0; j<p->N; j++) {
	if(i == j)
	    continue;
	for(int k=0; k<3; k++)
	    r[k] = ions[i]->x[k] - ions[j]->x[k];
	U += OOFPEN*ions[i]->Z/sqrt(dot(r,r));
    }
    return U;
}
