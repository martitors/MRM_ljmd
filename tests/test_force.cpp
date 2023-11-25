// unit test example with test fixture
#include "gtest/gtest.h"
#include "verlet.h"
#include "types.h"
#include "allocate.h"
#include <random>


TEST(ForceCalculation, ComputesForces) {
    // Create two small molecular system with 3 and 4 particles
    mdsys_t sys1;
    sys1.natoms = 4;
    sys1.mass = 1.0;
    sys1.epsilon = 2.0;
    sys1.sigma = 1.0;
    sys1.rcut = 4.0;
    sys1.box = 10.0;

    mdsys_t sys2;
    sys2.natoms = 3;
    sys2.mass = 2.0;
    sys2.epsilon = 1.0;
    sys2.sigma = 1.0;
    sys2.rcut = 3.0;
    sys2.box = 8.0;

    // Allocate memory
    allocate(&sys1);
    allocate(&sys2);


    // Initialize positions, velocities, and forces
    sys1.rx[0] = 0.5; sys1.ry[0] = 11.0; sys1.rz[0] = 2.0;
    sys1.rx[1] = 2.5; sys1.ry[1] = 0.0; sys1.rz[1] = 3.0;
    sys1.rx[2] = 3.0; sys1.ry[2] = 0.0; sys1.rz[2] = 2.0;
    sys1.rx[3] = 0.0; sys1.ry[3] = 5.0; sys1.rz[3] = 0.0;

    sys1.vx[0] = 0.0; sys1.vy[0] = 0.0; sys1.vz[0] = 0.0;
    sys1.vx[1] = 0.0; sys1.vy[0] = 1.0; sys1.vz[1] = 0.0;
    sys1.vx[2] = 0.0; sys1.vy[2] = 0.0; sys1.vz[2] = 0.0;
    sys1.vx[3] = 0.0; sys1.vy[3] = 1.0; sys1.vz[3] = 0.0;
  

    sys2.rx[0] = 1.0; sys2.ry[0] = 0.0; sys2.rz[0] = 2.0;
    sys2.rx[1] = 0.0; sys2.ry[1] = 0.0; sys2.rz[1] = 5.0;
    sys2.rx[2] = 0.0; sys2.ry[2] = 2.0; sys2.rz[2] = 1.0;

    sys2.vx[0] = 0.0; sys2.vy[0] = 0.0; sys2.vz[0] = 0.0;
    sys2.vx[1] = 0.0; sys2.vy[1] = 0.0; sys2.vz[1] = 0.0;
    sys2.vx[2] = 0.0; sys2.vy[2] = 0.0; sys2.vz[2] = 0.0;

    //initialise forces to zero
    azzero(sys1.fx, sys1.natoms);
    azzero(sys1.fy, sys1.natoms);
    azzero(sys1.fz, sys1.natoms);
    azzero(sys2.fx, sys2.natoms);
    azzero(sys2.fy, sys2.natoms);
    azzero(sys2.fz, sys2.natoms);

    // Call the force function
    force(&sys1);
    force(&sys2);

    //Check that f of particles extern to the cutoff is zero
    ASSERT_DOUBLE_EQ(sys1.fx[3],0.0);
    ASSERT_DOUBLE_EQ(sys1.fy[3],0.0);
    ASSERT_DOUBLE_EQ(sys1.fz[3],0.0);
    ASSERT_DOUBLE_EQ(sys2.fx[1],0.0);
    ASSERT_DOUBLE_EQ(sys2.fy[1],0.0);
    ASSERT_DOUBLE_EQ(sys2.fz[1],0.0);
    
    // Check the computed forces against expected values
    ASSERT_DOUBLE_EQ(sys1.fx[0],0.11659418934778376);
    ASSERT_DOUBLE_EQ(sys2.fy[2],-0.036694101508916367);
    ASSERT_DOUBLE_EQ(sys2.fz[0],-sys2.fz[2]);



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
}

