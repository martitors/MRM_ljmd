/* compute forces */
#include "../include/force_compute.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void force(mdsys_t *sys)
{
    double ffac,c12,c6,rcsq,rsq;
    double rx,ry,rz;
    int i,j, offset;

    double *fx, *fy, *fz; // local pointers for openmp
    double epot=0.0; // needed for reduction with openmp

    /* zero energy and forces */
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
    rcsq = sys->rcut * sys->rcut;

    #ifdef _OPENMP
    #pragma omp parallel private(fx,fy,fz,i,j,rx,ry,rz,rsq,rsq_inv,r6,ffac,offset) reduction(+:epot)
    #endif
    {

        // who is who
        int tid; // thread id
        int nthreads; // number of threads

        #ifdef _OPENMP
        	tid = omp_get_thread_num();
        	nthreads = omp_get_num_threads();
		int rest = sys->natoms % nthreads;
                if (nthreads < rest) {
                	i_omp = sys->natoms/nthreads + 1;
       		        fromidx = tid * i_omp;
         		toidx = (tid + 1) * i_omp;
    		} else {
        		i_omp = sys->natoms/nthreas;
        		fromidx = tid * i_omp + rest;
        		toidx = (tid + 1) * i_omp + rest;
    		}
        #else
        	tid = 0;
        	nthreads = 1;
		i_omp = sys->natoms;
    		fromidx = 0;
    		toidx = sys->natoms;
        #endif

        // assign slices of force array to threads
        fx = sys->fx + (tid * sys->natoms);
        fy = sys->fy + (tid * sys->natoms);
        fz = sys->fz + (tid * sys->natoms);

        // zero energy and forces
        azzero(fx,sys->natoms);
        azzero(fy,sys->natoms);
        azzero(fz,sys->natoms);

	#ifdef _OPENMP
	#ifndef CHUNKSIZE
	#define CHUNKSIZE 1
	#endif
	#pragma omp for schedule(guided,CHUNKSIZE)
	#endif

    	for(i=fromidx; i < toidx ; ++i) {
        	for(j=i+1; j < (sys->natoms); ++j) {

            /* particles have no interactions with themselves */
            //if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double r6,rinv; rinv=1.0/rsq; r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                sys->epot += r6*(c12*r6 - c6);


                sys->fx[i] += rx*ffac;          sys->fx[j] -= rx*ffac;
                sys->fy[i] += ry*ffac;          sys->fy[j] -= ry*ffac;
                sys->fz[i] += rz*ffac;          sys->fz[j] -= rz*ffac;
            }
        }
    }
	#ifdef _OPENMP
	#pragma omp for
	#endif
        for (i=0; i<sys->natoms; ++i) {
            for (offset=sys->natoms; offset<nthreads*sys->natoms; offset+=sys->natoms) {
                sys->fx[i] += sys->fx[offset + i];
                sys->fy[i] += sys->fy[offset + i];
                sys->fz[i] += sys->fz[offset + i];
            }
        }
    }//
     sys->epot = epot;
}

void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    for (i=0; i<sys->natoms; ++i) {
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
