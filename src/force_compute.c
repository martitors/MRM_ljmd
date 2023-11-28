/* compute forces */
#include "../include/force_compute.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void force(mdsys_t *sys)
{

    double epot; 
    epot =0.0 ;
    //int fromidx, toidx;

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:epot) 
	#endif
    {    


        int tid ;
        double rx1, ry1, rz1;
        double rx, ry, rz;
        double *fx, *fy, *fz;
        //double *f;
        double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
        double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
        double rcsq = sys->rcut * sys->rcut;
        double ffac,rsq;
        int j, fromidx ,toidx;

        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #else
        	tid = 0;
        #endif


        fx = sys->fx + (tid * sys->natoms);
        fy = sys->fy + (tid * sys->natoms);
        fz = sys->fz + (tid * sys->natoms);

        // zero forces
        azzero(fx,sys->natoms);
        azzero(fy,sys->natoms);
        azzero(fz,sys->natoms);

    	//inside the box defined based on thread ID
	    //for(int ii=fromidx; ii < toidx ; ++ii) {,
        for (int i = 0; i < sys->natoms - 1; i += sys->nthreads)
        {   int ii= i + tid;
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


            int i = (sys->natoms/sys->nthreads) + 1;
            fromidx = tid * i;
            toidx = fromidx  + i;
            if (toidx > sys->natoms) toidx = sys->natoms;
                for (int i=1; i< sys->nthreads; ++i) {
                int offs = i*sys->natoms;
                for (int j=fromidx; j < toidx; ++j) {
                    sys->fx[j] += sys->fx[offs + j];
                    sys->fy[j] += sys->fy[offs + j];
                    sys->fz[j] += sys->fz[offs + j];
            		}
        	}

    }//end of parallel

    //printf("%f\n",epot);
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
