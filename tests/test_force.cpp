// unit test example with test fixture
#include "gtest/gtest.h"
#include "types.h"
#include "utilities.h"
#include "allocate.h"
#include "force_compute.h"

#ifdef _MPI
#include "mpi.h"
#endif  

TEST(ForceCalculation, ComputesForces) {
    // Create two small molecular system with 3 and 4 particles
    
    mdsys_t sys1;
    sys1.natoms = 3;
    sys1.mass = 1.0;
    sys1.epsilon = 2.0;
    sys1.sigma = 1.0;
    sys1.rcut = 4.0;
    sys1.box = 10.0;

    // mdsys_t sys2;
    // sys2.natoms = 3;
    // sys2.mass = 2.0;
    // sys2.epsilon = 1.0;
    // sys2.sigma = 1.0;
    // sys2.rcut = 3.0;
    // sys2.box = 8.0;

    int mpirank = 0;

    #ifdef _MPI
    
    MPI_Init( NULL, NULL ); 

    MPI_Comm_size( MPI_COMM_WORLD, &sys1.npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &sys1.rank ); 

    mpirank = sys1.rank;
    #endif

    // Allocate memory

    allocate(&sys1);
    //allocate(&sys2);

    #if defined(_MPI)
    sys1.cx = (double*) malloc( sys1.natoms * sizeof(double) );
    sys1.cy = (double*) malloc( sys1.natoms * sizeof(double) );
    sys1.cz = (double*) malloc( sys1.natoms * sizeof(double) );
    // sys2.cx = (double*) malloc( sys2.natoms * sizeof(double) );
    // sys2.cy = (double*) malloc( sys2.natoms * sizeof(double) );
    // sys2.cz = (double*) malloc( sys2.natoms * sizeof(double) );
    #endif

    if (mpirank==0){
    // Initialize positions, velocities, and forces
    sys1.rx[0] = 0.5; sys1.ry[0] = 11.0; sys1.rz[0] = 2.0;
    sys1.rx[1] = 2.5; sys1.ry[1] = 0.0; sys1.rz[1] = 3.0;
    sys1.rx[2] = 3.0; sys1.ry[2] = 0.0; sys1.rz[2] = 2.0;
    //sys1.rx[3] = 0.0; sys1.ry[3] = 5.0; sys1.rz[3] = 0.0;

    // sys2.rx[0] = 1.0; sys2.ry[0] = 0.0; sys2.rz[0] = 2.0;
    // sys2.rx[1] = 0.0; sys2.ry[1] = 0.0; sys2.rz[1] = 5.0;
    // sys2.rx[2] = 0.0; sys2.ry[2] = 2.0; sys2.rz[2] = 1.0;

    azzero(sys1.vx, sys1.natoms);
    azzero(sys1.vy, sys1.natoms);
    azzero(sys1.vz, sys1.natoms);
    // azzero(sys2.vx, sys2.natoms);
    // azzero(sys2.vy, sys2.natoms);
    // azzero(sys2.vz, sys2.natoms);
    
    azzero(sys1.fx, sys1.natoms);
    azzero(sys1.fy, sys1.natoms);
    azzero(sys1.fz, sys1.natoms);
    // azzero(sys2.fx, sys2.natoms);
    // azzero(sys2.fy, sys2.natoms);
    // azzero(sys2.fz, sys2.natoms);
    }

    //initialise forces to zero
    #ifdef _MPI
    azzero(sys1.cx, sys1.natoms);
    azzero(sys1.cy, sys1.natoms);
    azzero(sys1.cz, sys1.natoms);
    // azzero(sys2.cx, sys2.natoms);
    // azzero(sys2.cy, sys2.natoms);
    // azzero(sys2.cz, sys2.natoms);
    #endif

    // Call the force function
    force(&sys1);
    // force(&sys2);

    //Check that f of particles extern to the cutoff is zero
    if (mpirank==0){

    //ASSERT_DOUBLE_EQ(sys1.fx[3],0.0);
    //ASSERT_DOUBLE_EQ(sys1.fy[3],0.0);
    //ASSERT_DOUBLE_EQ(sys1.fz[3],0.0);
    //ASSERT_DOUBLE_EQ(sys2.fx[1],0.0);
    //ASSERT_DOUBLE_EQ(sys2.fy[1],0.0);
    //ASSERT_DOUBLE_EQ(sys2.fz[1],0.0);
    
    // Check the computed forces against expected values
    //ASSERT_NEAR(sys1.fx[0],0.11659418934778376,1e-5);
    //ASSERT_NEAR(sys2.fy[2],-0.036694101508916326,1e-5);
    //ASSERT_NEAR(sys2.fz[0],-sys2.fz[2],1e-5);
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

    // free(sys2.rx);
    // free(sys2.ry);
    // free(sys2.rz);
    // free(sys2.vx);
    // free(sys2.vy);
    // free(sys2.vz);
    // free(sys2.fx);
    // free(sys2.fy);
    // free(sys2.fz);

    #if defined(_MPI)
    
    // free(sys2.cx);
    // free(sys2.cy);
    // free(sys2.cz);
    free(sys1.cx);
    free(sys1.cy);
    free(sys1.cz);

    MPI_Finalize();
    #endif
}

