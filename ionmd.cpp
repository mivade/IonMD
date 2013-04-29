#include <cstdio>
#include <cmath>
#include <omp.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include "ionmd.hpp"

using namespace std;

#define dbg(wat) (std::cout << wat << "\n")

// Random number generator
const gsl_rng_type *rng_T = gsl_rng_mt19937;
gsl_rng *rng = gsl_rng_alloc(rng_T);

//---------------------//
//--UTILITY FUNCTIONS--//
//---------------------//

void printIonStatistics(Params *p) {
    for (int i = 0; i < p->N; i++) {
        printf("Ion %d: m = %.1f, Z = %.1f\n", i+1, p->m[i]/amu, p->Z[i]/q_e);
    }
}

void printParams(Params *p) {
    printf("Total number of ions: %d\n", p->N);
    printf("Time parameters:\n");
    printf("   dt = %.1e, t_max = %.1e, t_steps = %d\n\n", p->dt, p->t_max, p->t_steps);
    printf("Laser parameters:\n");
    printf("   lambda = %.1f, delta = %.1e*Gamma, Gamma = 2*pi*%.1e, s0 = %.1f\n\n", p->lmbda/1e-9, p->delta/p->Gamma, p->Gamma/(2*pi), p->s0);
    printf("Trap parameters:\n");
    printf("   V = %.1f, U = %.1f, UEC = %.1f\n   r0 = %.2e, z0 = %2e, kappa = %.1e\n", p->V, p->U, p->UEC, p->r0, p->z0, p->kappa);
    printf("   Vsec = %.1f, w = 2*pi*%.1f\n\n", p->Vsec, p->w/(2*pi));
    printf("Background gas parameters:\n");
    printf("   gamma_col = %.5f\n\n", p->gamma_col);
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

Ion *initIon(double *x0, double *v0, int i, Params *p) {
    Ion *ion = new Ion;
    copyVector(ion->x, x0);
    copyVector(ion->v, v0);
    ion->m = p->m[i];
    ion->Z = p->Z[i];
    zeroVector(ion->a);
    ion->index = i;
    ion->lc = p->lc[i];
    return ion;
}

void updateIon(Ion *ion, Ion **ions, double t, double *Fcoullist, Params *p) {
    int i;
    double F[3], Ft[3], Fl[3], Fc[3], Fsec[3], Fs[3], a[3];
    zeroVector(F);
    for(i=0; i<3; i++)
	ion->x[i] = ion->x[i] + ion->v[i]*p->dt + 0.5*ion->a[i]*pow(p->dt,2);
    FTrap(ion, t, p,Ft);
    if(p->use_laser && ion->lc) {
	//if(ion->lc || (!ion->lc && t < p->traj_start))
	FLaser(ion, p, Fl);
    }
    else
        zeroVector(Fl);
    if(p->use_coulomb)
        copyVector(Fc, &Fcoullist[ion->index*3]); //alternately: Fcoullist+ion->index*3
    else
        zeroVector(Fc);
    if(p->use_secular)
        FSecular(ion, t, p, Fsec);
    else
        zeroVector(Fsec);
    if(p->use_stochastic)
        FStochastic(ion, p, Fs);
    else
        zeroVector(Fs);
    for(int i = 0; i < 3; i++) {
        F[i] += Ft[i] + Fl[i] + Fc[i] + Fsec[i] + Fs[i];
	a[i] = F[i]/ion->m;
	//ion->x[i] = ion->x[i] + ion->v[i]*p->dt + 0.5*ion->a[i]*pow(p->dt,2);
	ion->v[i] = ion->v[i] + 0.5*(ion->a[i] + a[i])*p->dt;
	//ion->x[i] += ion->v[i]*p->dt;
        //ion->v[i] += F[i]*p->dt/ion->m;
    }
    copyVector(ion->a, a);
}

//----------//
//--FORCES--//
//----------//

// Trap force
void FTrap(Ion *ion, double t, Params *p, double *F) { //result stored in F
    double A, B;
    A = p->kappa*p->UEC/pow(p->z0,2);
    if (p->use_rfmm) {
        B = 4*p->V*cos(p->Omega*t)/pow(p->r0,2);
        F[0] = ion->Z*(A - B)*ion->x[0];
        F[1] = ion->Z*(A + B)*ion->x[1];
    } else {
        B = 2*ion->Z*pow(p->V,2)/(ion->m*pow(p->Omega,2)*pow(p->r0,4));
        F[0] = ion->Z*(A - B)*ion->x[0];
        F[1] = ion->Z*(A - B)*ion->x[1];
    }
    F[2] = -2*ion->Z*A*ion->x[2];
}

// Laser cooling
// (result stored in F)
void FLaser(Ion *ion, Params *p, double *F) {
    zeroVector(F);
    double F0, k, s0, Gamma, beta, delta;
    k = 2*pi/p->lmbda;
    s0 = p->s0;
    Gamma = p->Gamma;
    delta = p->delta - k*dot(p->khat, ion->v);
    beta = -HBAR*pow(k,2)*4*s0*delta/Gamma/pow(1+s0+pow(2*delta/Gamma,2),2);
    for (int i = 0; i < 3; i++) {
        F0 = HBAR*k*p->khat[i]*Gamma/2/(s0/(s0 + 1));
        //F0 = 0;
        F[i] = F0 - beta*ion->v[i];
    }
}

// Coulomb interaction
// (result stored in F)
void FCoulomb(Ion *ion, Ion **ions, Params *p, double *F) { 
    double r[3];
    zeroVector(r);
    zeroVector(F);
    double r_mag = 0;
    int i = ion->index;
    int j,k;
    for(j = 0; j < p->N; j++) {
        if(i == j)
            continue;
        for(k = 0; k < 3; k++)
            r[k] = (ion->x[k] - ions[j]->x[k]);
        r_mag = sqrt(dot(r,r));
        for(k = 0; k < 3; k++)
            F[k] += OOFPEN*ion->Z*ions[j]->Z*r[k]/pow(r_mag,3);
    }
}

// Secular excitations
// (result stored in F)
void FSecular(Ion *ion, double t, Params *p, double *F) {
    zeroVector(F);
    F[0] = ion->Z*p->Vsec*ion->x[0]*cos(p->w*t);
}

// Stochastic processes, e.g. collisions with background gases
// (result stored in F)
void FStochastic(Ion *ion, Params *p, double *F) {
    zeroVector(F);
    if(gsl_rng_uniform(rng) <= exp(-p->gamma_col*p->dt))
       return;	// no collision
    double hat[] = {gsl_rng_uniform(rng), gsl_rng_uniform(rng), gsl_rng_uniform(rng)};
    normalize(hat);
    double hbark = HBAR*2*pi/p->lmbda;
    for(int i=0; i<3; i++)
	F[i] = hbark*hat[i]/p->dt;
}

// Precomputes all Coulomb interactions (that way they can be applied
// all at once instead of having one ion updated before the Coulomb
// interaction is computed for the rest).
// (result stored in Flist)
void allCoulomb(Ion **ions, Params *p, double *Flist) { 
    int i;
    #pragma omp parallel for
    for (i = 0; i < p->N; i++)
        FCoulomb(ions[i], ions, p,&Flist[i*3]); //Alternately: Flist+i*3
}

//--------------//
//--SIMULATION--//
//--------------//

int simulate(double *x0, double *v0, Params *p) {
    // RNG seeding
    gsl_rng_set(rng, time(0));

    //every %3 element is the start of a new vector 
    //keeps from having to reallocate every time
    //don't have to manually manage the memory for each vector
    double *Fclist = new double[p->N*3]; 
    int i,j;
    float tmp;
    
    // Initialize ions
    Ion **ions = new Ion*[p->N];
    for(i = 0; i < p->N; i++) {
        ions[i] = initIon(&x0[i*3], &v0[i*3], i, p);
    }
    
    // Data recording initialization
    FILE *traj_file = fopen(p->traj_fname, "wb");
    FILE *ccd_file = fopen(p->ccd_fname, "wb");
    
    // Run simulation
    omp_set_num_threads(p->num_threads);
    int t_i = 0;
    int t_10 = (int)(p->t_max/p->dt)/10;
    for(double t=0; t<p->t_max; t+=p->dt) {
        // Calculate Coulomb forces
        if (p->use_coulomb)
            allCoulomb(ions, p, Fclist);

	// Progress update
	if (t_i % t_10 == 0 && !p->quiet)
	    printf("%d%% complete; t = %f\n", (int)10*t_i/t_10, t);

        // Update each ion
        for (i = 0; i < p->N; i++) {
            // Record data
            if (p->record_traj) {
                if(t > p->traj_start) {
		    if(i == 0) {
			for (j=0; j<3; j++) {
			    tmp = (float)(ions[i]->x[j]/1e-3);
			    fwrite(&tmp, sizeof(float), 1, traj_file);
			}
                    }
                    simCCDPoint(ccd_file, ions[i]);
                }
            }
            // Update
            updateIon(ions[i], ions, t, Fclist, p);
        }
        t_i++;
    }

    // Save data
    FILE *fpos_file = fopen(p->fpos_fname, "w");
    FILE *fvel_file = fopen("fvel.txt", "w");
    fprintf(fpos_file, "%d\n", p->N);
    fprintf(fpos_file, "Simulated ion crystal, positions in microns\n");
    for (i = 0; i < p->N; i++) {
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
    fclose(fvel_file);
    fclose(ccd_file);
    for (i = 0; i < p->N; i++)
        delete ions[i];
    delete[] ions;
    delete[] Fclist;
    gsl_rng_free(rng);
    return 1;
}
