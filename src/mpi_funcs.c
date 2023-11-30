#include "types.h"
#include "../include/mpi_funcs.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef _MPI
#include "mpi.h"
#endif //_MPI

#ifdef _OPENMP
#include <omp.h>
#endif

void init_mpi_omp(int argc, char **argv,mdsys_t *sys){

    #ifdef _MPI
    
    MPI_Init( &argc, &argv ); 
    MPI_Comm_size( MPI_COMM_WORLD, &sys->npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &sys->rank );    

    if (sys->rank==0) printf("MPI correctly initialized\n");
    #else
    sys->npes=1;
    sys->rank=0;
    #endif

    #ifdef _OPENMP
        sys->nthreads = omp_get_max_threads();
    #else
        sys->nthreads = 1;
    #endif
}


void mpi_bcasts(mdsys_t *sys){

    
    #ifdef _MPI

    MPI_Bcast(&sys->natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys->epsilon, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys->sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys->rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys->box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys->nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    sys->cx = (double*) malloc( sys->nthreads * sys->natoms * sizeof(double) );
    sys->cy = (double*) malloc( sys->nthreads * sys->natoms * sizeof(double) );
    sys->cz = (double*) malloc( sys->nthreads * sys->natoms * sizeof(double) );
    
    #endif
}

void mpi_fin(mdsys_t *sys){
 #ifdef _MPI
    free( sys->cx ); 
    free( sys->cy ); 
    free( sys->cz ); 

  MPI_Finalize();
 #endif
}