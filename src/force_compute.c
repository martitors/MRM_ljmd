/* compute forces */
#include "../include/force_compute.h"
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _MPI
#include "mpi.h"
#endif //_MPI

void force(mdsys_t *sys)
{
    double ffac,c12,c6,rcsq,rsq;
    double rx,ry,rz;
    int ii;
    int i,j, nthreads;
    
    /* zero energy and forces */
    double epot = 0.0;

    #ifdef _MPI
    azzero( sys->cx, sys->natoms );
    azzero( sys->cy, sys->natoms );
    azzero( sys->cz, sys->natoms );
    
    MPI_Bcast( sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    #else
    sys->epot=0.0;

    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

    #endif

    c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);
    rcsq = sys->rcut * sys->rcut;
    
#ifdef _MPI
    for(i=0; i < (sys->natoms - 1); i += sys->npes) {
            ii = i + sys->rank;
            if (ii >= (sys->natoms - 1)) break;

#elif defined(_OPENMP)
    #pragma omp parallel default(shared) private(ii,j,rx,ry,rz,ffac)
    {
    double *omp_fx = sys->fx;
    double *omp_fy = sys->fy;
    double *omp_fz = sys->fz;
    #pragma omp for reduction(+:epot, omp_fx[:sys->natoms],omp_fy[:sys->natoms],omp_fz[:sys->natoms])
    for(ii=0; ii < (sys->natoms); ++ii) {
#else

    for(ii=0; ii < (sys->natoms); ++ii) {
#endif
        for(j=ii+1; j < (sys->natoms); ++j) {
            /* particles have no interactions with themselves */
            //if (i==j) continue;

            /* get distance between particle i and j */
            rx=pbc(sys->rx[ii] - sys->rx[j], 0.5*sys->box);
            ry=pbc(sys->ry[ii] - sys->ry[j], 0.5*sys->box);
            rz=pbc(sys->rz[ii] - sys->rz[j], 0.5*sys->box);
            rsq = rx*rx + ry*ry + rz*rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double r6,rinv; rinv=1.0/rsq; r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                epot += r6*(c12*r6 - c6);

            #ifdef _MPI

                sys->cx[ii] += rx*ffac;          sys->cx[j] -= rx*ffac;
                sys->cy[ii] += ry*ffac;          sys->cy[j] -= ry*ffac;
                sys->cz[ii] += rz*ffac;          sys->cz[j] -= rz*ffac;

            #elif defined(_OPENMP)

                omp_fx[ii] += rx*ffac;          omp_fx[j] -= rx*ffac;
                omp_fy[ii] += ry*ffac;          omp_fy[j] -= ry*ffac;
                omp_fz[ii] += rz*ffac;          omp_fz[j] -= rz*ffac;

            #else

                sys->fx[ii] += rx*ffac;          sys->fx[j] -= rx*ffac;
                sys->fy[ii] += ry*ffac;          sys->fy[j] -= ry*ffac;
                sys->fz[ii] += rz*ffac;          sys->fz[j] -= rz*ffac;
            
            #endif
            }
        }
    }

    #ifdef _MPI
    MPI_Reduce( sys->cx, sys->fx, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( sys->cy, sys->fy, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( sys->cz, sys->fz, sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    #elif defined (_OPENMP)
    }
    sys->epot = epot;
    #else
    sys->epot = epot;
    #endif
  
}



