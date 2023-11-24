// unit test example with test fixture
#include "gtest/gtest.h"
#include "verlet.h"

class VerletTest: public ::testing::Test {

protected:

    mdsys_t *sys;

    void SetUp()
    {
        sys = new mdsys_t;
        sys->natoms = 2;
        sys->mass = 1.0;
        sys->rx = new double[2];
        sys->vx = new double[2];
        sys->fx = new double[2];
        sys->rx[0] = -1.0;
        sys->rx[1] = 1.0;
        sys->vx[0] = 0.0;
        sys->vx[1] = 0.0;
        sys->fx[0] = 1.0;
        sys->fx[1] = 0.2;
    }

    void TearDown()
        {
            delete[] sys->rx;
            delete[] sys->vx;
            delete[] sys->fx;

            delete sys;
        }
};

TEST_F(VerletTest, step1)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0],-1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0],0.0);

    // Save initial values for comparison
    double initial_rx_0 = sys->rx[0];
    double initial_vx_0 = sys->vx[0];

    // Perform the velocity Verlet step 1
    velverlet_1(sys);

    // Calculate the expected values after the step
    double expected_rx_0 = initial_rx_0 + 0.5 * sys->dt * initial_vx_0;
    double expected_vx_0 = initial_vx_0 + 0.5 * sys->dt / mvsq2e * sys->fx[0] / sys->mass;

    // Check if the values match the expectations
    ASSERT_DOUBLE_EQ(sys->rx[0], expected_rx_0);
    ASSERT_DOUBLE_EQ(sys->vx[0], expected_vx_0);
}

TEST_F(VerletTest, step2)
{
    ASSERT_NE(sys,nullptr);
    ASSERT_DOUBLE_EQ(sys->rx[0],-1.0);
    ASSERT_DOUBLE_EQ(sys->vx[0],0.0);

    // Save initial values for comparison
    double initial_vx_0 = sys->vx[0];

    // Perform the velocity Verlet step 2
    velverlet_2(sys);

    // Calculate the expected value of vx[0] after the step
    double expected_vx_0 = initial_vx_0 + 0.5 * sys->dt / mvsq2e * sys->fx[0] / sys->mass;

    // Check if the values match the expectations
    ASSERT_DOUBLE_EQ(sys->rx[0], -1.0); // rx[0] should remain unchanged in step 2
    ASSERT_DOUBLE_EQ(sys->vx[0], expected_vx_0);
}
