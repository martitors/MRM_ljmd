// unit test example with test fixture
#include "gtest/gtest.h"
#include "types.h"
#include "utilities.h"
#include "allocate.h"
#include <random>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef _CELL
#include "cell.h"
#else
#include "force_compute.h"
#endif

#ifdef _MPI
#include "mpi.h"
#endif  

// Test case for the force function
TEST(ForceCalculation, ComputesForces) {
    // Create a molecular system with 3 particles
    
    mdsys_t sys1;
    sys1.natoms = 4;
    sys1.epsilon = 2.0;
    sys1.sigma = 1.0;
    sys1.rcut = 4.0;
    sys1.box = 20.0;


    mdsys_t sys2;
    sys2.natoms = 3;
    sys2.mass = 2.0;
    sys2.epsilon = 1.0;
    sys2.sigma = 1.0;
    sys2.rcut = 3.0;
    sys2.box = 15.0;
    

    #ifdef _OPENMP
        sys1.nthreads = omp_get_max_threads();
        sys2.nthreads = omp_get_max_threads();
    #else
        sys1.nthreads =1;
        sys2.nthreads =1;
    #endif

    // Allocate memory
    #ifdef _MPI
    MPI_Init( NULL, NULL ); 

    MPI_Comm_size( MPI_COMM_WORLD, &sys1.npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &sys1.rank ); 
    MPI_Comm_size( MPI_COMM_WORLD, &sys2.npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &sys2.rank ); 
    #else
        sys1.rank = 0;
        sys1.npes = 1;
        sys2.rank = 0;
        sys2.npes = 1;
    #endif

    // Allocate memory for the system
    allocate(&sys1);
    allocate(&sys2);

    // Initialize positions
    if (sys1.rank == 0){

    sys1.rx[0] = 0.5; sys1.ry[0] = 0.5; sys1.rz[0] = 0.5;
    sys1.rx[1] = 2.5; sys1.ry[1] = 0.0; sys1.rz[1] = 3.0;
    sys1.rx[2] = 10.0; sys1.ry[2] = 0.0; sys1.rz[2] = 0.0;
    sys1.rx[3] = 5.0; sys1.ry[3] = 0.0; sys1.rz[3] = 0.0;
    }
    if (sys2.rank == 0){
    sys2.rx[0] = 1.0; sys2.ry[0] = 0.0; sys2.rz[0] = 2.0;
    sys2.rx[1] = 0.0; sys2.ry[1] = 0.0; sys2.rz[1] = 5.0;
    sys2.rx[2] = 0.0; sys2.ry[2] = 2.0; sys2.rz[2] = 1.0;
    }


    #ifdef _MPI
    sys1.cx = (double*) malloc( sys1.nthreads*sys1.natoms * sizeof(double) );
    sys1.cy = (double*) malloc( sys1.nthreads*sys1.natoms * sizeof(double) );
    sys1.cz = (double*) malloc( sys1.nthreads*sys1.natoms * sizeof(double) );
    sys2.cx = (double*) malloc( sys2.nthreads*sys2.natoms * sizeof(double) );
    sys2.cy = (double*) malloc( sys2.nthreads*sys2.natoms * sizeof(double) );
    sys2.cz = (double*) malloc( sys2.nthreads*sys2.natoms * sizeof(double) );
    #endif
    // Call the force function

    #ifdef _CELL
    sys1.clist = NULL;
    sys1.plist = NULL;
    sys2.clist = NULL;
    sys2.plist = NULL;
    build_cells(&sys1);
    build_cells(&sys2);
    update_cells(&sys1);
    update_cells(&sys2);
    #endif

    force(&sys1);
    force(&sys2);

    // Check forces against expected values
    if (sys1.rank==0){

    ASSERT_DOUBLE_EQ(sys1.fx[2],0.0);
    ASSERT_DOUBLE_EQ(sys1.fy[2],0.0);
    ASSERT_DOUBLE_EQ(sys1.fz[2],0.0);
    ASSERT_DOUBLE_EQ(sys2.fx[1],0.0);
    ASSERT_DOUBLE_EQ(sys2.fy[1],0.0);
    ASSERT_DOUBLE_EQ(sys2.fz[1],0.0);
    
    // Check the computed forces against expected values
    ASSERT_NEAR(sys1.fx[0], 0.0078842986764635671,1.e-5);
    ASSERT_NEAR(sys2.fy[2],-0.036694101508916367,1.e-5);
    ASSERT_NEAR(sys2.fz[0],-sys2.fz[2],1.e-5);
    }


    // Clean up: free memory
    free(sys1.rx);
    free(sys1.ry);
    free(sys1.rz);
    free(sys1.vx);
    free(sys1.vy);
    free(sys1.vz);
    free(sys1.fx);
    free(sys1.fy);
    free(sys1.fz);

    free(sys2.rx);
    free(sys2.ry);
    free(sys2.rz);
    free(sys2.vx);
    free(sys2.vy);
    free(sys2.vz);
    free(sys2.fx);
    free(sys2.fy);
    free(sys2.fz);

    #ifdef _MPI
    free(sys1.cx);
    free(sys1.cy);
    free(sys1.cz);
    free(sys2.cx);
    free(sys2.cy);
    free(sys2.cz);
    #endif

    #ifdef _CELL
    free_cells(&sys1);
    free_cells(&sys2);
    #endif
}
