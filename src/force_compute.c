/* compute forces */
#include "../include/force_compute.h"
#include <omp.h>

void force(mdsys_t *sys)
{
    double ffac,c12,c6,rcsq,rsq;
    double rx,ry,rz;
    int i,j, nthreads;
    
    /* zero energy and forces */
    sys->epot=0.0;
    double epot_tmp = 0.0;

    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
    rcsq = sys->rcut * sys->rcut;

    #if defined(_OPENMP)
    //double *omp_fx = sys->fx;
    //double *omp_fy = sys->fy;
    //double *omp_fz = sys->fz;

    #pragma omp parallel default(shared) private(i,j,rx,ry,rz,ffac)
    #endif
    {

    double *omp_fx = sys->fx;
    double *omp_fy = sys->fy;
    double *omp_fz = sys->fz;

    /*#if defined(_OPENMP)
    int tid=omp_get_thread_num();
    #else
    int tid=0;
    #endif*/

    //omp_fx=sys->fx + (tid*sys->natoms); azzero(omp_fx,sys->natoms);
    //omp_fy=sys->fy + (tid*sys->natoms); azzero(omp_fy,sys->natoms);
    //omp_fz=sys->fz + (tid*sys->natoms); azzero(omp_fz,sys->natoms);

    #pragma omp for reduction(+:epot_tmp, omp_fx[:sys->natoms],omp_fy[:sys->natoms],omp_fz[:sys->natoms])
    	
     for(i=0; i < (sys->natoms) -1 ; ++i) {
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
                epot_tmp += r6*(c12*r6 - c6);


                omp_fx[i] += rx*ffac;          omp_fx[j] -= rx*ffac;
                omp_fy[i] += ry*ffac;          omp_fy[j] -= ry*ffac;
                omp_fz[i] += rz*ffac;          omp_fz[j] -= rz*ffac;
            }
        }
    }

    /*for(i=0 ; i < sys->natoms ; ++i) {
        	for(j=0; j < (sys->natoms); ++j) {

        	    	/* particles have no interactions with themselves */
            		//if (i==j) continue;

            		/* get distance between particle i and j */
            		/*rx=pbc(sys->rx[i] - sys->rx[j], 0.5*sys->box);
            		ry=pbc(sys->ry[i] - sys->ry[j], 0.5*sys->box);
            		rz=pbc(sys->rz[i] - sys->rz[j], 0.5*sys->box);
            		r = sqrt(rx*rx + ry*ry + rz*rz);

            		/* compute force and energy if within cutoff */
            		/*if (r < sys->rcut) {
                		ffac = -4.0*sys->epsilon*(-12.0*pow(sys->sigma/r,12.0)/r
                                         +6*pow(sys->sigma/r,6.0)/r);

                		epot_tmp += 0.5*4.0*sys->epsilon*(pow(sys->sigma/r,12.0)
                                               -pow(sys->sigma/r,6.0));

				omp_fx[i] += rx / r * ffac;
                		omp_fy[i] += ry / r * ffac;
                		omp_fz[i] += rz / r * ffac;
            		}
        	}
    	}*/

	}//parallel region
    sys->epot = epot_tmp;
    
    /*#pragma omp critical
        {
            for (i = 0; i < sys->natoms; ++i) {
                sys->fx[i] += omp_fx[i];
                sys->fy[i] += omp_fy[i];
                sys->fz[i] += omp_fz[i];
            }
        }*/
  
}

void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    double ekin_tmp = 0.0;

#if defined(_OPENMP)
#pragma omp parallel for reduction(+ : ekin_tmp)
#endif
    for (i=0; i<sys->natoms; ++i) {
        ekin_tmp += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
    }
    sys->ekin = ekin_tmp;
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
