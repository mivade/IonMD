#include <cstdio>
#include <cmath>
#include <omp.h>
#include "boost/multi_array.hpp"
#include "ionmd.hpp"

using namespace std;

//---------------------//
//--UTILITY FUNCTIONS--//
//---------------------//

double *zeroVector() {
    double *v = new double[3];
    for(int i=0; i<3; i++)
	v[i] = 0;
    return v;
}

double dot(double *a, double *b) {
    double c = 0;
    for(int i=0; i<3; i++)
	c += a[i]*b[i];
    return c;
}

void printIonStatistics(Params *p) {
    for(int i=0; i<p->N; i++) {
	printf("Ion %d: m = %.1f, Z = %.1f\n", i+1, p->m[i]/amu, p->Z[i]/q_e);
    }
}

void printParams(Params *p) {
    printf("Total number of ions: %d\n", p->N);
    printf("Time parameters:\n");
    printf("   dt = %.1e, t_max = %.1e, t_steps = %d\n\n",
	   p->dt, p->t_max, p->t_steps);
    printf("Laser parameters:\n");
    printf("   lambda = %.1f, delta = %.1e*Gamma, Gamma = 2*pi*%.1e, s0 = %.1f\n\n",
	   p->lmbda/1e-9, p->delta/p->Gamma, p->Gamma/(2*pi), p->s0);
    printf("Trap parameters:\n");
    printf("   V = %.1f, U = %.1f, UEC = %.1f\n   r0 = %.2e, z0 = %2e, kappa = %.1e\n",
	   p->V, p->U, p->UEC, p->r0, p->z0, p->kappa);
    printf("   Vsec = %.1f, w = 2*pi*%.1f\n\n", p->Vsec, p->w/(2*pi));
    printf("Simulation options:\n");
    printf("   RF micromotion: %s\n", (p->use_rfmm ? "on" : "off"));
    printf("   Coulomb interaction: %s\n", (p->use_coulomb ? "on" : "off"));
    printf("   Laser interaction: %s\n", (p->use_laser ? "on" : "off"));
    printf("   Secular excitation: %s\n", (p->use_secular ? "on" : "off"));
    printf("   Stochastic forces: %s\n", (p->use_stochastic ? "on" : "off"));
    return;
}

// Write point to the CCD file for creating a simulated CCD image
// later. Mass given in amu, x and y points in microns.
void simCCDPoint(FILE *fout, Ion *ion) {
    float m = (float)(ion->m/amu),
	x = (float)(dot(ion->x, xhat)/1e-6),
	y = (float)(dot(ion->x, yhat)/1e-6);
    fwrite(&m, sizeof(float), 1, fout);
    fwrite(&x, sizeof(float), 1, fout);
    fwrite(&y, sizeof(float), 1, fout);
    return;
}

//-----------------//
//--ION FUNCTIONS--//
//-----------------//

Ion *initIon(double *x0, double*v0, int i, Params *p) {
    Ion *ion = new Ion;
    ion->x = new double[3];
    ion->v = new double[3];
    for(int j=0; j<3; j++) {
    	ion->x[j] = x0[j];
    	ion->v[j] = v0[j];
    }
    ion->m = p->m[i];
    ion->Z = p->Z[i];
    ion->index = i;
    ion->lc = p->lc[i];
    return ion;
}

void updateIon(Ion *ion, Ion **ions, double t, double **Fcoul, Params *p) {
    double *F, *Ft, *Fl, *Fc, *Fsec, *Fs;
    F = zeroVector();
    Ft = FTrap(ion, t, p);
    if(p->use_laser && ion->lc)
	Fl = FLaser(ion, p);
    else
	Fl = zeroVector();
    if(p->use_coulomb)
	Fc = Fcoul[ion->index];
    else
	Fc = zeroVector();
    if(p->use_secular)
	Fsec = FSecular(ion, t, p);
    else
	Fsec = zeroVector();
    if(p->use_stochastic)
	Fs = FStochastic(ion, p);
    else
	Fs = zeroVector();
    for(int i=0; i<3; i++) {
        F[i] += Ft[i] + Fl[i] + Fc[i] + Fsec[i] + Fs[i];
        ion->v[i] += F[i]*p->dt/ion->m;
        ion->x[i] += ion->v[i]*p->dt;
    }
    delete[] F;
    delete[] Ft;
    delete[] Fl;
    delete[] Fc;
    delete[] Fs;
}

//----------//
//--FORCES--//
//----------//

// Trap force
double *FTrap(Ion *ion, double t, Params *p) {
    double *F = new double[3];
    double A, B;
    A = p->kappa*p->UEC/pow(p->z0,2);
    if(p->use_rfmm) {
	B = 4*p->V*cos(p->Omega*t)/pow(p->r0,2);
	F[0] = ion->Z*(A - B)*ion->x[0];
	F[1] = ion->Z*(A + B)*ion->x[1];
    }
    else {
	B = -2*ion->Z*pow(p->V,2)/(ion->m*pow(p->Omega,2)*pow(p->r0,4));
	F[0] = ion->Z*(A + B)*ion->x[0];
	F[1] = ion->Z*(A + B)*ion->x[1];
    }
    F[2] = -2*ion->Z*A*ion->x[2];
    return F;
}

// Laser cooling
double *FLaser(Ion *ion, Params *p) {
    double *F = zeroVector();
    double F0, k, s0, Gamma, beta, delta;
    k = 2*pi/p->lmbda;
    s0 = p->s0;
    Gamma = p->Gamma;
    delta = p->delta - k*dot(p->khat, ion->v);
    beta = -HBAR*pow(k,2)*4*s0*delta/Gamma/pow(1+s0+pow(2*delta/Gamma,2),2);
    for(int i=0; i<3; i++) {
	//F0 = HBAR*k*p->khat[i]*Gamma/(s0/(s0 + 1));
	//F0 = 2*pi*HBAR*c*s0*Gamma/(3*pow(p->lmbda,3)*pow(p->r_l,2))*p->khat[i];
	//F0 = HBAR*k*p->khat[i]*Gamma/2/(1 + s0 + pow(2*delta/Gamma,2));
	F0 = 0;
	F[i] = F0 - beta*ion->v[i];
    }
    return F;
}

// Coulomb interaction
double *FCoulomb(Ion *ion, Ion **ions, Params *p) {
    double *F = zeroVector(), *r = zeroVector();
    double r_mag = 0;
    int i = ion->index;
    int j,k;
    for(j=0; j<p->N; j++) {
        if(i == j)
            continue;
	for(k=0; k<3; k++)
            r[k] = (ion->x[k] - ions[j]->x[k]);
        r_mag = sqrt(dot(r,r));
        for(k=0; k<3; k++)
	    F[k] += OOFPEN*ion->Z*ions[j]->Z*r[k]/pow(r_mag,3);
    }
    delete[] r;
    return F;
}

// Secular excitations
double *FSecular(Ion *ion, double t, Params *p) {
    double *F = zeroVector();
    F[0] = ion->Z*p->Vsec*ion->x[0]*cos(p->w*t);
    return F;
}

// Stochastic processes, e.g. collisions with background gases
double *FStochastic(Ion *ion, Params *p) {
    double *F = zeroVector();
    return F;
}

// Precomputes all Coulomb interactions (that way they can be applied
// all at once instead of having one ion updated before the Coulomb
// interaction is computed for the rest).
double **allCoulomb(Ion **ions, Params *p) {
    int i;
    double **F = new double*[p->N];
    #pragma omp parallel for
    for(i=0; i<p->N; i++)
	F[i] = FCoulomb(ions[i], ions, p);
    return F;
}

//--------------//
//--SIMULATION--//
//--------------//

int simulate(double *x0, double *v0, Params *p) {
    double **Fc;
    int i,j;
    float tmp;

    // Initialize ions
    double *x = new double[3],
	*v = new double[3];
    Ion **ions = new Ion*[p->N];
    for(i=0; i<p->N; i++) {
	for(j=0; j<3; j++) {
	    x[j] = x0[3*i+j];
	    v[j] = v0[3*i+j];
	}
	ions[i] = initIon(x, v, i, p);
    }
    delete[] x;
    delete[] v;

    // Data recording initialization
    FILE *traj_file = fopen(p->traj_fname, "wb");
    FILE *ccd_file = fopen(p->ccd_fname, "wb");
    
    // Run simulation
    omp_set_num_threads(p->num_threads);
    int t_i = 0;
    for(double t=0; t<p->t_max; t+=p->dt) {
	// Calculate Coulomb forces
	if(p->use_coulomb)
	    Fc = allCoulomb(ions, p);

	// Update each ion
	for(i=0; i<p->N; i++) {
	    // Record data
	    if(p->record_traj) {
		if(t > p->traj_start) {
		    for(j=0; j<3; j++) {
			tmp = (float)(ions[i]->x[j]/1e-3);
			fwrite(&tmp, sizeof(float), 1, traj_file);
		    }
		    simCCDPoint(ccd_file, ions[i]);
		}
	    }

	    // Update
	    updateIon(ions[i], ions, t, Fc, p);
	}
	
	// Cleanup
	if(p->use_coulomb) {
	    // for(i=0; i<p->N; i++)
	    // 	delete[] Fc[i];
	    delete[] Fc;
	}
	t_i++;
    }

    // Save data
    FILE *fpos_file = fopen(p->fpos_fname, "w");
    FILE *fvel_file = fopen("fvel.txt", "w");
    fprintf(fpos_file, "%d\n", p->N);
    fprintf(fpos_file, "Simulated ion crystal, positions in microns\n");
    for(i=0; i<p->N; i++) {
	fprintf(fpos_file, "%d ", (int)(ions[i]->m/amu));
	//fprintf(fpos_file, "BA ");
	for(j=0; j<3; j++) {
	    fprintf(fpos_file, "%e ", ions[i]->x[j]/1e-6);
	    fprintf(fvel_file, "%e ", ions[i]->v[j]);
	}
	fprintf(fpos_file, "\n");
	fprintf(fvel_file, "\n");
    }

    // Cleanup
    fclose(traj_file);
    fclose(fpos_file);
    fclose(ccd_file);
    for(i=0; i<p->N; i++)
	delete ions[i];
    delete[] ions;
    // SEGFAULT!
    // if(p->use_coulomb) {
    // 	for(i=0; i<p->N; i++)
    // 	    delete[] Fc[i];
    // 	delete[] Fc;
    // }// FIX!!!!!!!!
    return 1;
}
