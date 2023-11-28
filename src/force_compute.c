/* compute forces */
#include "../include/force_compute.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void force(mdsys_t *sys)
{
    double ffac,rsq;
    int i,j, offs, fromidx ,toidx;

    double epot=0.0; // needed for reduction with openmp

    // zero energy and force
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);
    

    double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
    double rcsq = sys->rcut * sys->rcut;
    int tid ;

    #if defined(_OPENMP)
    #pragma omp parallel private(i,j,rcsq,offs) reduction(+:epot) 
	#endif
      { 
        int i;
        double *fx, *fy, *fz; // local pointers for openmp
        #ifdef _OPENMP
            epot=0.0;
        	int nthreads = omp_get_max_threads();
            tid = omp_get_thread_num();
            
        #else
        	tid = 0;
        	int nthreads = 1;
        #endif

        // assign slices of force array to threads
        fx = sys->fx + (tid * sys->natoms);
        fy = sys->fy + (tid * sys->natoms);
        fz = sys->fz + (tid * sys->natoms);

        // zero forces
        azzero(fx,sys->natoms);
        azzero(fy,sys->natoms);
        azzero(fz,sys->natoms);

	double rx1, ry1, rz1;
	double rx, ry, rz;

    	//inside the box defined based on thread ID
	//for(int ii=fromidx; ii < toidx ; ++ii) {,
    for(int i=0; i < (sys->natoms -1); i += nthreads) {
        int ii = i + tid;
        if (ii >= (sys->natoms -1)) break;
		        rx1=sys->rx[ii];
                ry1=sys->ry[ii];
                rz1=sys->rz[ii];

        	for(int jj=ii+1; jj < sys->natoms ; ++jj) {

            	/* particles have no interactions with themselves */
            	//if (i==j) continue;

        		/* get distance between particle i and j */
            		rx=pbc(rx1 - sys->rx[jj], 0.5*sys->box);
            		ry=pbc(ry1 - sys->ry[jj], 0.5*sys->box);
            		rz=pbc(rz1 - sys->rz[jj], 0.5*sys->box);
            		rsq = rx*rx + ry*ry + rz*rz;

            		/* compute force and energy if within cutoff */
            		if (rsq < rcsq) {
                		double r6,rinv; rinv=1.0/rsq; r6=rinv*rinv*rinv;
                		ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                		epot += r6*(c12*r6 - c6);


                		fx[ii] += rx*ffac;          fx[jj] -= rx*ffac;
                		fy[ii] += ry*ffac;          fy[jj] -= ry*ffac;
                		fz[ii] += rz*ffac;          fz[jj] -= rz*ffac;
            		}
        	}
    	}

	#if defined (_OPENMP)
	#pragma omp barrier
	#endif
           printf("tid: %d  numthreads:%d \n",tid,nthreads);
            i = (sys->natoms/nthreads) + 1;
            fromidx = tid * i;
            toidx = (tid + 1) * i;
            if (toidx > sys->natoms) toidx = sys->natoms;
           for (int i=1; i< nthreads; ++i) {
             int offs = i*sys->natoms;
			 for (int j=fromidx; j < toidx; ++j) {
				        sys->fx[i] += sys->fx[offs + i];
                		sys->fy[i] += sys->fy[offs + i];
                		sys->fz[i] += sys->fz[offs + i];
            		}
        	}

    }//end of parallel
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
