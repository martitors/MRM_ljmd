/* compute forces */
#include "../include/force_compute.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void force(mdsys_t *sys)
{
    double ffac,rsq;
    int i,j, offset;

    double epot=0.0; // needed for reduction with openmp

    // zero energy and force
    sys->epot=0.0;
    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    #ifdef _OPENMP
    #pragma omp parallel private(i,j,rsq,rsq_inv,r6,ffac,offset) reduction(+:epot)
    {   
	#endif
    	double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    	double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
    	double rcsq = sys->rcut * sys->rcut;
	   
        double *fx, *fy, *fz; // local pointers for openmp
        
        int i, fromidx, toidx;  //interval variables
        int tid;        // thread id
        int nthreads;   // number of threads

        #ifdef _OPENMP
        	tid = omp_get_thread_num();
        	nthreads = omp_get_num_threads();
		int rest = sys->natoms % nthreads;
                if (tid < rest) {
                	i = sys->natoms/nthreads + 1;
       		        fromidx = tid * i;
         		toidx = (tid + 1) * i;
    		} else {
        		i = sys->natoms/nthreads;
        		fromidx = tid * i + rest;
        		toidx = (tid + 1) * i + rest;
    		}
        #else
        	tid = 0;
        	nthreads = 1;
		i = sys->natoms;
    		fromidx = 0;
    		toidx = sys->natoms;
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
	for(int ii=fromidx; ii < toidx ; ++ii) {
		rx1=sys->rx[ii];
                ry1=sys->ry[ii];
                rz1=sys->rz[ii];

        	for(int jj=ii+1; jj < toidx ; ++jj) {

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


                		sys->fx[ii] += rx*ffac;          sys->fx[jj] -= rx*ffac;
                		sys->fy[ii] += ry*ffac;          sys->fy[jj] -= ry*ffac;
                		sys->fz[ii] += rz*ffac;          sys->fz[jj] -= rz*ffac;
            		}
        	}
    	}

	//with outside the box defined based on the thread ID
	for(int ii=fromidx; ii<toidx ; ++ii) {
		int kk = ii + toidx;
		rx1=sys->rx[ii];
                ry1=sys->ry[ii];
                rz1=sys->rz[ii];
		
		for(int jj=kk ; jj<sys->natoms ; ++jj) {
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


                                sys->fx[ii] += rx*ffac;          sys->fx[jj] -= rx*ffac;
                                sys->fy[ii] += ry*ffac;          sys->fy[jj] -= ry*ffac;
                                sys->fz[ii] += rz*ffac;          sys->fz[jj] -= rz*ffac;
                        }
                }
        }



	#if defined (_OPENMP)
	#pragma omp barrier
	#endif
                for (i=1; i < nthreads; ++i) {
			int offs = i*sys->natoms;
			for (int j=fromidx; j < toidx; ++j) {
				sys->fx[i] += sys->fx[offset + i];
                		sys->fy[i] += sys->fy[offset + i];
                		sys->fz[i] += sys->fz[offset + i];
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
