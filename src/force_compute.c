/* compute forces */
#include "../include/force_compute.h"
#ifdef _OPENMP
#include <omp.h>
#endif

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

    #ifdef _OPENMP
    #pragma omp parallel default(shared) private(i,j,rx,ry,rz,ffac)
    {
    double *omp_fx = sys->fx;
    double *omp_fy = sys->fy;
    double *omp_fz = sys->fz;
    #pragma omp for reduction(+:epot_tmp, omp_fx[:sys->natoms],omp_fy[:sys->natoms],omp_fz[:sys->natoms])
    #endif	
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
                #ifdef _OPENMP
                epot_tmp += r6*(c12*r6 - c6);
                #else
                sys->epot += r6*(c12*r6 - c6);
                #endif

                #ifdef _OPENMP
                omp_fx[i] += rx*ffac;          omp_fx[j] -= rx*ffac;
                omp_fy[i] += ry*ffac;          omp_fy[j] -= ry*ffac;
                omp_fz[i] += rz*ffac;          omp_fz[j] -= rz*ffac;
                #else
                sys->fx[i] += rx*ffac;          sys->fx[j] -= rx*ffac;
                sys->fy[i] += ry*ffac;          sys->fy[j] -= ry*ffac;
                sys->fz[i] += rz*ffac;          sys->fz[j] -= rz*ffac;
                #endif
            }
        }
    }
    #ifdef _OPENMP
    }
    sys->epot = epot_tmp;
    #endif
  
}

void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    double ekin_tmp = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : ekin_tmp)
#endif
    for (i=0; i<sys->natoms; ++i) {
        #ifdef _OPENMP
        ekin_tmp += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
        #else
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
        #endif
    }
    #ifdef _OPENMP
    sys->ekin = ekin_tmp;
    #endif
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}
