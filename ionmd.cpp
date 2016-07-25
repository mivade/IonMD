#include <cstdio>
#include <cmath>
#include <vector>
#include <omp.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <boost/format.hpp>
#include "ionmd.hpp"

using namespace std;

#define dbg(wat) (std::cout << wat << "\n")
#define err(e) (cerr << "ERROR: " << e << endl)
const vec xhat = {-sqrt(2)/2, sqrt(2)/2, 0}, yhat = {0,0,1};

// Random number generator
const gsl_rng_type *rng_T = gsl_rng_mt19937;
gsl_rng *rng = gsl_rng_alloc(rng_T);

//-----------------//
//--ION FUNCTIONS--//
//-----------------//

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

// Precomputes all Coulomb interactions (that way they can be applied
// all at once instead of having one ion updated before the Coulomb
// interaction is computed for the rest).
// (result stored in Flist)
//void allCoulomb(Ion **ions, Params *p, mat *Flist) {
mat allCoulomb(Ion **ions, Params *p) {
    int i;
    mat Flist(3, p->N);
    #pragma omp parallel for
    for(i=0; i<p->N; i++) {
        //FCoulomb(ions[i], ions, p, &Flist[i*3]); //Alternately: Flist+i*3
	Flist.col(i) = FCoulomb(ions[i], ions, p);
    }
    return Flist;
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
    mat Fclist(3, p->N);
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
    for(double t = 0; t < p->t_max; t += p->dt) {
        // Calculate Coulomb forces
        if (p->use_coulomb)
            Fclist = allCoulomb(ions, p);

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
    for(i=0; i<p->N_masses; i++)
	gsl_histogram2d_free(ccd[i]);
    if(p->use_abort && abort == 1)
	return 0;
    return 1;
}
