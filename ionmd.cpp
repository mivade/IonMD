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
#include <cmath>
#include <vector>
#include <omp.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "ionmd.hpp"
#include "minimize.hpp"

using namespace std;
using arma::vec;
using arma::dot;
using arma::mat;

#define dbg(wat) (std::cout << wat << "\n")
#define err(e) (cerr << "ERROR: " << e << endl)
const vec xhat = {-sqrt(2)/2, sqrt(2)/2, 0}, yhat = {0,0,1};

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
    printf("   lambda = %.1f, delta = %.1e*Gamma, Gamma = 2*pi*%.1e, s = %.1f\n", p->lmbda/1e-9, p->delta/p->Gamma, p->Gamma/(2*pi), p->s);
    printf("   khat = [%.1f, %.1f, %.1f]\n\n", p->khat[0], p->khat[1], p->khat[2]);
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
    printf("   Abort on ion out of bounds: %s\n", (p->use_abort ? "on" : "off"));
    if(p->use_abort)
	printf("   Abort bounds: %.3e\n", p->abort_bounds);
    return;
}

void printPositions(Ion *ion) {
    printf("%e %e %e\n", ion->x[0], ion->x[1], ion->x[2]);
}

//-----------------//
//--ION FUNCTIONS--//
//-----------------//

Ion *initIon(double *x0, double *v0, int i, Params *p) {
    Ion *ion = new Ion;
    ion->x = vec(3);
    ion->v = vec(3);
    ion->a = arma::zeros<vec>(3);
    copyVector(ion->x, x0);
    copyVector(ion->v, v0);
    ion->m = p->m[i];
    ion->Z = p->Z[i];
    ion->index = i;
    ion->lc = p->lc[i];
    return ion;
}

void simCCDPoint(Ion *ion, gsl_histogram2d **ccd, Params *p) {
    if(!p->sim_ccd)
	return;
    int j = -1;
    if(p->N_masses == 1)
	j = 0;
    else {
	for(int i=0; i<p->N_masses; i++) {
	    if(ion->m == p->masses[i])
		j = i;
	}
    }
    if(j == -1) {
	err("Mass list not configured correctly. Not simulating CCD images.");
	p->sim_ccd = 0;
	return;
    }
    double x = dot(ion->x, xhat)/1e-6,
	y = dot(ion->x, yhat)/1e-6;
    gsl_histogram2d_increment(ccd[j], x, y);
    return;
}

void updateIon(Ion *ion, Ion **ions, double t, double *Fcoullist, Params *p) {
    // TODO: These all need to be length 3
    vec F, Ft, Fl, Fc, Fsec, Fs, a;
    F.zeros();
    ion->x += ion->v*p->dt + 0.5*ion->a*pow(p->dt, 2);
    FTrap(ion, t, p, &Ft);
    if((p->use_laser && ion->lc) || p->minimizing)
	FLaser(ion, p, &Fl);
    else
        Fl.zeros();
    if(p->use_coulomb)
        copyVector(Fc, &Fcoullist[ion->index*3]); //alternately: Fcoullist+ion->index*3
    else
        Fc.zeros();
    if(p->use_secular)
        FSecular(ion, t, p, &Fsec);
    else
        Fsec.zeros();
    if(p->use_stochastic)
        FStochastic(ion, p, &Fs);
    else
        Fs.zeros();
    F = Ft + Fl + Fc + Fsec + Fs;
    a = F/ion->m;
    ion->v += 0.5*(ion->a + a)*p->dt;
    ion->a = a;
}

//----------//
//--FORCES--//
//----------//

// Trap force
// (result stored in F)
void FTrap(Ion *ion, double t, Params *p, vec *F) {
    double A, B;
    //A = p->kappa*p->UEC/pow(p->z0,2);
    if (p->use_rfmm) {
        //B = 4*p->V*cos(p->Omega*t)/pow(p->r0,2);
	// TODO: Fix
	A = 0; B = 0;
        F[0] = ion->Z*(A - B)*ion->x[0];
        F[1] = ion->Z*(A + B)*ion->x[1];
    }
    else {
	A = ion->Z*pow(p->V,2)/(ion->m*pow(p->Omega,2)*pow(p->r0,4));
        B = p->kappa*p->UEC/(2*pow(p->z0,2));
        F[0] = -2*ion->Z*(A-B)*ion->x[0];
        F[1] = -2*ion->Z*(A-B)*ion->x[1];
    }
    F[2] = -4*ion->Z*B*ion->x[2];
}

// Laser cooling
// (result stored in F)
void FLaser(Ion *ion, Params *p, vec *F) {
    F->zeros();
    double beta, F0;
    beta = p->beta;
    F0 = p->F0;
    if(p->minimizing) {
	beta = 1e-20; // unrealistically large damping for minimizing
	if(ion->lc == 0) // don't use constant pressure term on non-lc'ed ions
	    F0 = 0;
    }
    for (int i = 0; i < 3; i++) {
        F[i] = F0*p->khat[i] - beta*ion->v[i];
    }
}

// Coulomb interaction
// (result stored in F)
void FCoulomb(Ion *ion, Ion **ions, Params *p, vec *F) { 
    vec r = vec(3);
    r.zeros();
    F->zeros();
    double r_mag = 0;
    int i = ion->index, j, k;
    //#pragma omp parallel for
    for(j = 0; j < p->N; j++) {
        if(i == j)
            continue;
	r = ion->x - ions[j]->x;
        r_mag = sqrt(dot(r, r));
        for(k = 0; k < 3; k++)
	    F[k] += OOFPEN*ion->Z*ions[j]->Z*r[k]/pow(r_mag, 3);
	// F += OOFPEN*ion->Z*ions[j]->Z*r/pow(r_mag, 3);
    }
}

// Secular excitations
// (result stored in F)
void FSecular(Ion *ion, double t, Params *p, vec *F) {
    F->zeros();
    F[0] = ion->Z*p->Vsec*ion->x[0]*cos(p->w*t);
}

// Stochastic processes, e.g. collisions with background gases
// (result stored in F)
void FStochastic(Ion *ion, Params *p, vec *F) {
    double v;
    F->zeros();
    //if(gsl_rng_uniform(rng) <= exp(-p->gamma_col*p->dt))
    //return;	// no collision
    vec hat = {gsl_rng_uniform(rng),
	       gsl_rng_uniform(rng),
	       gsl_rng_uniform(rng)};
    normalize(hat);
    v = sqrt(2*kB*p->gamma_col*p->dt/ion->m);
    for(int i=0; i<3; i++)
	F[i] = ion->m*v*hat[i]/p->dt;
}

// Precomputes all Coulomb interactions (that way they can be applied
// all at once instead of having one ion updated before the Coulomb
// interaction is computed for the rest).
// (result stored in Flist)
void allCoulomb(Ion **ions, Params *p, mat *Flist) { 
    int i, j;
    #pragma omp parallel for
    for(i=0; i<p->N; i++) {
        //FCoulomb(ions[i], ions, p, &Flist[i*3]); //Alternately: Flist+i*3
	FCoulomb(ions[i], ions, p, &Flist->col(i));
    }
}

//--------------//
//--SIMULATION--//
//--------------//

int simulate(double *x0, double *v0, Params *p) {
    // RNG seeding
    gsl_rng_set(rng, time(0));
    
    // Set number of threads for multiprocessing
    omp_set_num_threads(p->num_threads);

    // Every %3 element is the start of a new vector 
    // Keeps from having to reallocate every time
    // Don't have to manually manage the memory for each vector
    //double *Fclist = new double[p->N*3];
    mat *Fclist = new mat(3, p->N);
    int i, j,
	abort = 0,
	T_ctr = 0;
    double vT = 0;
    float tmp;
    vec xcom(3);

    // Initialize CCD
    gsl_histogram2d *ccd[p->N_masses];
    for(i=0; i<p->N_masses; i++) {
	ccd[i] = gsl_histogram2d_alloc(p->ccd_bins, p->ccd_bins);
	gsl_histogram2d_set_ranges_uniform(ccd[i],
					   -p->ccd_extent, p->ccd_extent,
					   -p->ccd_extent, p->ccd_extent);
    }
    
    // Initialize ions
    float M = 0;
    Ion **ions = new Ion*[p->N];
    for(i=0; i<p->N; i++) {
        ions[i] = initIon(&x0[i*3], &v0[i*3], i, p);
	M += (float)ions[i]->m;
    }

    // Do minimization routine
    printf("Minimizing...\n");
    minimize(ions, p);
    //p->minimizing = 0;
    //minimize(x0, ions, p);

    // Reset initial velocities
    for(i=0; i<p->N; i++) {
	ions[i]->v.zeros();
	//copyVector(ions[i]->v, &v0[i*3]);
    }
    
    // Data recording initialization
    FILE *traj_file = fopen(p->traj_fname, "wb");
    FILE *com_file = fopen(p->com_fname, "wb");
    FILE *temp_file = fopen(p->temp_fname, "w");
    
    // Run simulation
    int t_i = 0;
    int t_10 = (int)(p->t_max/p->dt)/10;
    printf("Simulating...\n");
    for(double t=0; t<p->t_max; t+=p->dt) {
        // Calculate Coulomb forces
        if (p->use_coulomb)
            allCoulomb(ions, p, Fclist);

	// Progress update
	if (t_i % t_10 == 0 && !p->quiet)
	    printf("%d%% complete; t = %f\n", (int)10*t_i/t_10, t);

        // Update each ion
	xcom.zeros();
        for(i=0; i<p->N; i++) {
            // Record data
            if(p->record_traj) {
                if(t > p->traj_start) {
		    for(j=0; j<3; j++) {
			tmp = (float)(ions[i]->x[j]/1e-3);
			xcom[j] += (float)ions[i]->m*tmp/M;
			if(i == 0)
			    fwrite(&tmp, sizeof(float), 1, traj_file);
                    }
		    simCCDPoint(ions[i], ccd, p);
                }
            }

	    // Update "temperature"
	    // TODO: DTRT
	    // TODO: arma: fix this
	    for(j=0; j<3; j++)
		vT += pow(ions[i]->v[j], 2);
	    if(T_ctr == p->T_steps) {
		T_ctr = 0;
		vT = sqrt(vT)/p->T_steps/p->N;
		fprintf(temp_file, "%e %e\n", t, vT);
	    }
	    else T_ctr++;

            // Update
            updateIon(ions[i], ions, t, Fclist, p);

	    // Check bounds
	    if(p->use_abort) {
		for(j=0; j<2; j++) {
		    if(abs(ions[i]->x[j]) >= p->abort_bounds || ions[i]->x[j] != ions[i]->x[j]) {
			err("Ion out of bounds! Aborting...");
			fprintf(stderr, "x[%i] = %.3e\n",
				j, ions[i]->x[j]);
			abort = 1;
		    }
		}
		if(abort) break;
	    }
        }
	if(p->record_traj) {
	    for(j=0; j<3; j++) {
		float tmp = float(xcom[j]);
		fwrite(&tmp, sizeof(float), 1, com_file);
	    }
	}
        t_i++;
	if(abort) break;
    }

    // Save data
    FILE *fpos_file = fopen(p->fpos_fname, "w");
    FILE *fvel_file = fopen("fvel.txt", "w");
    fprintf(fpos_file, "%d\n", p->N);
    fprintf(fpos_file, "Simulated ion crystal, positions in microns\n");
    for (i=0; i<p->N; i++) {
        fprintf(fpos_file, "%d ", (int)(ions[i]->m/amu));
        for(j=0; j<3; j++) {
            fprintf(fpos_file, "%e ", ions[i]->x[j]/1e-6);
            fprintf(fvel_file, "%e ", ions[i]->v[j]);
        }
        fprintf(fpos_file, "\n");
        fprintf(fvel_file, "\n");
    }
    if(p->sim_ccd) {
	for(i=0; i<p->N_masses; i++) {
	    char s[256];
	    sprintf(s, "%s_%i.dat", p->ccd_fname, i+1);
	    FILE *ccd_file = fopen(s, "wb");
	    gsl_histogram2d_fwrite(ccd_file, ccd[i]);
	    fclose(ccd_file);
	}
    }

    // Cleanup
    fclose(com_file);
    fclose(traj_file);
    fclose(fpos_file);
    fclose(fvel_file);
    fclose(temp_file);
    for (i = 0; i < p->N; i++)
        delete ions[i];
    delete[] ions;
    delete[] Fclist;
    for(i=0; i<p->N_masses; i++)
	gsl_histogram2d_free(ccd[i]);
    if(p->use_abort && abort == 1)
	return 0;
    return 1;
}
