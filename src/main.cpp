#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "types.h"
#include "constants.h"
#include "utilities.h"
#include "verlet.h"
#include "output.h"
#include "input.h"
#ifdef _CELL
#include "cell.h"
#else
#include "force_compute.h"
#endif
#include "cleanup.h"
#include "allocate.h"

#ifdef _MPI
#include "mpi.h"
#endif //_MPI

#ifdef _OPENMP
#include <omp.h>
#endif

#define LJMD_VERSION 1.0

int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    #ifdef _MPI
    
    MPI_Init( &argc, &argv ); 

    MPI_Comm_size( MPI_COMM_WORLD, &sys.npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &sys.rank );    

    if (sys.rank==0) printf("MPI correctly initialized\n");
    #else
    sys.npes=1;
    sys.rank=0;
    #endif

    #ifdef _OPENMP
        sys.nthreads = omp_get_max_threads();
    #else
        sys.nthreads =1;
    #endif


    if (sys.rank==0){
    printf("LJMD version %3.1f\n", LJMD_VERSION);
    
    printf("\nnthreads: %d\n", sys.nthreads);
    t_start = wallclock();

    read_input(line, restfile, trajfile, ergfile, &sys, &nprint);


    }

    #ifdef _MPI

    MPI_Bcast(&sys.natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.epsilon, 1, MPI_DOUBLE, 0,MPI_COMM_WORLD);
    MPI_Bcast(&sys.sigma, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.rcut, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.box, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sys.nsteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
    sys.cx = (double*) malloc( sys.nthreads * sys.natoms * sizeof(double) );
    sys.cy = (double*) malloc( sys.nthreads * sys.natoms * sizeof(double) );
    sys.cz = (double*) malloc( sys.nthreads * sys.natoms * sizeof(double) );
    
    #endif

    allocate(&sys);

    if(sys.rank == 0){

        /* read restart */
        char restPath[256]; // Adjust the size based on your needs
        snprintf(restPath, sizeof(restPath), "../examples/%s", restfile);
    

        fp=fopen(restPath,"r"); 
        if(fp) {
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.rx+i, sys.ry+i, sys.rz+i);
            }
            for (i=0; i<sys.natoms; ++i) {
                fscanf(fp,"%lf%lf%lf",sys.vx+i, sys.vy+i, sys.vz+i);
            }
            fclose(fp);
        } else {
            perror("cannot read restart file");
            return 3;
        }
    }
    /* initialize forces and energies.*/
    sys.nfi=0;
    #ifdef _CELL
    sys.clist = NULL;
    sys.plist = NULL;
    build_cells(&sys);
    update_cells(&sys);
    #endif
    force(&sys);

    if(sys.rank == 0){
        ekin(&sys);

        #ifdef _CELL
        update_cells(&sys);
        #endif
        char ergPath[256]; // Adjust the size based on your needs
        snprintf(ergPath, sizeof(ergPath), "../examples/%s", ergfile);

        char trajPath[256]; // Adjust the size based on your needs
        snprintf(trajPath, sizeof(trajPath), "../examples/%s", trajfile);

        erg=fopen(ergPath,"w");
        traj=fopen(trajPath,"w");

        printf("Startup time: %10.3fs\n", wallclock()-t_start);
        printf("Starting simulation with %d atoms for %d steps.\n",sys.natoms, sys.nsteps);
        printf("     NFI            TEMP            EKIN                 EPOT              ETOT\n");
        output(&sys, erg, traj);
        
        t_start = wallclock();

    }
    /* reset timer */
    /**************************************************/
    /* main MD loop */
    for(sys.nfi=1; sys.nfi <= sys.nsteps; ++sys.nfi) {
        /* write output, if requested */
        if ((sys.nfi % nprint) == 0 && sys.rank == 0)
            output(&sys, erg, traj);

        /* propagate system and recompute energies */
         if(sys.rank == 0)velverlet_1(&sys);
        force(&sys);   
         if(sys.rank == 0)velverlet_2(&sys);
         if(sys.rank == 0)ekin(&sys);
        #ifdef _CELL
        update_cells(&sys);
        #endif
        }
    /**************************************************/
    if(sys.rank == 0){
    /* clean up: close files, free memory */

    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
    #ifdef _CELL
    free_cells(&sys);
    #endif
    }
    cleanup(sys, erg, traj);
    #ifdef _MPI
    free( sys.cx ); 
    free( sys.cy ); 
    free( sys.cz ); 

    MPI_Finalize();
    #endif
    return 0;
}
