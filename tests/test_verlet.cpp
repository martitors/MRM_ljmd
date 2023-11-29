// unit test example with test fixture
#include "gtest/gtest.h"
#include "verlet.h"
#include "allocate.h"


TEST(Verlet, TimeIntegration){

    // Create a small molecular system with 2 particles 
    mdsys_t sys1;
    sys1.natoms = 2;
    sys1.mass = 1.0;
    sys1.epsilon = 2.0;
    sys1.sigma = 1.0;
    sys1.rcut = 4.0;
    sys1.box = 10.0;
    sys1.dt = 2.0;
    sys1.nthreads = 1;


    // Allocate memory
    allocate(&sys1);

    // Initialize positions, velocities, and forces
    sys1.rx[0] = 0.5; sys1.ry[0] = 1.5; sys1.rz[0] = 2.0;
    sys1.rx[1] = 2.5; sys1.ry[1] = 0.0; sys1.rz[1] = 3.0;

    sys1.vx[0] = 1.0; sys1.vy[0] = 3.0; sys1.vz[0] = 0.0;
    sys1.vx[1] = 0.0; sys1.vy[1] = 1.0; sys1.vz[1] = 2.0;
  
    sys1.fx[0] = 0.5; sys1.fy[0] = -1.0; sys1.fz[0] = 0.0;
    sys1.fx[1] = -0.5; sys1.fy[1] = 1.0; sys1.fz[1] = 0.0;

    // Call the verlet function

    velverlet_1(&sys1);

    // Check the computed forces against expected values
    ASSERT_DOUBLE_EQ(sys1.rx[0],2.5004183999999725);
    ASSERT_DOUBLE_EQ(sys1.ry[0],7.499163200000055);
    ASSERT_DOUBLE_EQ(sys1.rz[0],2.0);
    ASSERT_DOUBLE_EQ(sys1.vx[0],1.0002091999999863);
    ASSERT_DOUBLE_EQ(sys1.vy[0],2.9995816000000275);
    ASSERT_DOUBLE_EQ(sys1.vz[0],0.0);

    velverlet_1(&sys1);

    // Check the computed forces against expected values
    ASSERT_DOUBLE_EQ(sys1.rx[1],2.4987448000000825);
    ASSERT_DOUBLE_EQ(sys1.ry[1],4.002510399999835);
    ASSERT_DOUBLE_EQ(sys1.rz[1],11.0);
    ASSERT_DOUBLE_EQ(sys1.vx[1],-0.00041839999997254782);
    ASSERT_DOUBLE_EQ(sys1.vy[1],1.000836799999945);
    ASSERT_DOUBLE_EQ(sys1.vz[1],2.0);


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

}