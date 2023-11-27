
#include "gtest/gtest.h"
#include "../include/utilities.h"
#include <iostream>
#include <chrono>
#include <thread>

#ifdef _OPENMP
 #include <omp.h>
#endif

void sleep(double seconds){
    std::chrono::duration<double> duration(seconds);
    std::this_thread::sleep_for(duration);
    };

TEST(TestWallClock, one) {

    double t_test = -wallclock();
    sleep(0.1); // Change to a longer sleep duration
    t_test += wallclock();

    const double expected_min_time = 0.09; // Adjust these thresholds accordingly
    const double expected_max_time = 1.01;

    ASSERT_GE(t_test, expected_min_time);
    ASSERT_LE(t_test, expected_max_time);
}

TEST(TestAzzero, doubles)
{
    double *buf = (double*)malloc(10*sizeof(double));
    for (int i=0; i<10; ++i) buf[i]=i+1;
    ASSERT_DOUBLE_EQ(buf[1],2.0);
    ASSERT_DOUBLE_EQ(buf[5],6.0);
    ASSERT_DOUBLE_EQ(buf[9],10.0);

    azzero(buf,10);
    ASSERT_DOUBLE_EQ(buf[1],0.0);
    ASSERT_DOUBLE_EQ(buf[5],0.0);
    ASSERT_DOUBLE_EQ(buf[9],0.0);
    free(buf);
}
