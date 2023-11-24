#include "gtest/gtest.h"
#include "../include/force_compute.h" // Include the header file with function declarations

class ForceTest : public ::testing::Test {
protected:
    mdsys_t sys;

    void SetUp() override {
        sys.natoms = 2;
        sys.epsilon = 1.0;
        sys.sigma = 1.0;
        sys.box = 10.0;
        sys.rcut = 2.0;
        sys.rx = new double[2]{0.0, 1.5};
        sys.ry = new double[2]{0.0, 0.0};
        sys.rz = new double[2]{0.0, 0.0};
        sys.fx = new double[2]; // set to zero in the beginning
        sys.fy = new double[2]; // set to zero in the beginning
        sys.fz = new double[2]; // set to zero in the beginning

        // Set the forces to zero
        azzero(sys.fx, sys.natoms);
        azzero(sys.fy, sys.natoms);
        azzero(sys.fz, sys.natoms);
    }

    void TearDown() override {
        delete[] sys.rx;
        delete[] sys.ry;
        delete[] sys.rz;
        delete[] sys.fx;
        delete[] sys.fy;
        delete[] sys.fz;
    }
};

TEST_F(ForceTest, TwoParticles) {
    force(&sys);

    // Expected Values
    double expected_epot = -0.3203;
    double expected_fx[2] = {2.2627/2, 0.0};
    double expected_fy[2] = {0.0, 0.0};
    double expected_fz[2] = {0.0, 0.0};

    // Check computed values against expected values
    ASSERT_NEAR(sys.epot, expected_epot,1e-4);

    
        ASSERT_NEAR(sys.fx[0], expected_fx[0],1e-1);
        ASSERT_NEAR(sys.fy[0], expected_fy[0],1e-1);
        ASSERT_NEAR(sys.fz[0], expected_fz[0],1e-1);
    
}