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
#include "force_compute.h"
#include "cleanup.h"
#include "allocate.h"
#include "mpi_funcs.h"

#define LJMD_VERSION 1.0

int main(int argc, char **argv)
{
    int nprint, i;
    char restfile[BLEN], trajfile[BLEN], ergfile[BLEN], line[BLEN];
    FILE *fp,*traj,*erg;
    mdsys_t sys;
    double t_start;

    init_mpi_omp(argc,argv,&sys);


    if (sys.rank==0){
    printf("LJMD version %3.1f\n", LJMD_VERSION);
    
    printf("\nnthreads: %d\n", sys.nthreads);
    t_start = wallclock();

    read_input(line, restfile, trajfile, ergfile, &sys, &nprint);
    }

    mpi_bcasts(&sys);

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
    force(&sys);

    if(sys.rank == 0){
        ekin(&sys);

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
        velverlet_1(&sys);
        force(&sys);   
        velverlet_2(&sys);
        ekin(&sys);
        }
    /**************************************************/
    if(sys.rank == 0){
    /* clean up: close files, free memory */

    printf("Simulation Done. Run time: %10.3fs\n", wallclock()-t_start);
    }
    cleanup(sys, erg, traj);
    mpi_fin(&sys);
    return 0;
}
