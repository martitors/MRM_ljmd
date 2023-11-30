#ifndef MPI_FUNCS_H
#ifdef __cplusplus
extern "C"{
    #endif
    #define MPI_FUNCS_H
        #include <sys/time.h>
        #include <stdio.h>  
        #include "types.h"

           void init_mpi_omp(int argc, char **argv,mdsys_t *sys);
           void mpi_bcasts(mdsys_t *sys);
           void mpi_fin(mdsys_t *sys);

    #ifdef __cplusplus
   } 
#endif
#endif