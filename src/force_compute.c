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
    double epot = 0.0 ;

    /* zero energy and forces */
    sys->epot=0.0;

    #ifdef _MPI

    MPI_Bcast( sys->rx, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( sys->ry, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    MPI_Bcast( sys->rz, sys->natoms, MPI_DOUBLE, 0, MPI_COMM_WORLD );

    #endif

    #if defined(_OPENMP)
    #pragma omp parallel reduction(+:epot) 
	#endif
    {    

        int tid ;
        double rx1, ry1, rz1;
        double rx, ry, rz;
        double *fx, *fy, *fz;
        double *cx,*cy,*cz;

        // calculate pow() outside the loop
        double c12=4.0*sys->epsilon*pow(sys->sigma,12.0);
        double c6 =4.0*sys->epsilon*pow(sys->sigma, 6.0);

        double rcsq = sys->rcut * sys->rcut;   // to avoid expensive sqrt

        double ffac,rsq;
        int j, fromidx ,toidx;

        #ifdef _OPENMP
            tid = omp_get_thread_num();
        #else
        	tid = 0;
        #endif

        #ifdef _MPI

        cx = sys->cx + (tid * sys->natoms);
        cy = sys->cy + (tid * sys->natoms);
        cz = sys->cz + (tid * sys->natoms);

        azzero( cx, sys->natoms );
        azzero( cy, sys->natoms );
        azzero( cz, sys->natoms );    

        #else

        fx = sys->fx + (tid * sys->natoms);
        fy = sys->fy + (tid * sys->natoms);
        fz = sys->fz + (tid * sys->natoms);

        // zero forces
        azzero(fx,sys->natoms);
        azzero(fy,sys->natoms);
        azzero(fz,sys->natoms);

        #endif

        for(int k=0; k < (sys->natoms)-1; k += sys->npes*sys->nthreads) {
            int ii = k + sys->rank + (sys->npes*tid);
            if (ii >= (sys->natoms - 1)) break;
            
                    rx1=sys->rx[ii];
                    ry1=sys->ry[ii];
                    rz1=sys->rz[ii];

        for(j=ii+1; j < (sys->natoms); ++j) {

            rx=pbc(rx1 - sys->rx[j], 0.5*sys->box);
            ry=pbc(ry1 - sys->ry[j], 0.5*sys->box);
            rz=pbc(rz1 - sys->rz[j], 0.5*sys->box);

            rsq = rx*rx + ry*ry + rz*rz;

            /* compute force and energy if within cutoff */
            if (rsq < rcsq) {
                double r6,rinv; rinv=1.0/rsq; r6=rinv*rinv*rinv;
                ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                epot += r6*(c12*r6 - c6);

            #ifdef _MPI

                cx[ii] += rx*ffac;          cx[j] -= rx*ffac;
                cy[ii] += ry*ffac;          cy[j] -= ry*ffac;
                cz[ii] += rz*ffac;          cz[j] -= rz*ffac;

            #else 
                fx[ii] += rx*ffac;          fx[j] -= rx*ffac;
                fy[ii] += ry*ffac;          fy[j] -= ry*ffac;
                fz[ii] += rz*ffac;          fz[j] -= rz*ffac;
            
            #endif
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
                    #ifdef _MPI
                    sys->cx[j] += sys->cx[offs + j];
                    sys->cy[j] += sys->cy[offs + j];
                    sys->cz[j] += sys->cz[offs + j];
                    #else
                    sys->fx[j] += sys->fx[offs + j];
                    sys->fy[j] += sys->fy[offs + j];
                    sys->fz[j] += sys->fz[offs + j];
                    #endif
            		}
        	}

    }


    #ifdef _MPI
    MPI_Reduce( sys->cx, sys->fx, sys->nthreads*sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( sys->cy, sys->fy, sys->nthreads*sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( sys->cz, sys->fz, sys->nthreads*sys->natoms, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &epot, &sys->epot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    #else
    sys->epot = epot;
    #endif

}

