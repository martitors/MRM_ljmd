#include <gtest/gtest.h>
#include "../include/force_compute.h"

TEST(Kinetic_Energy, KE_2atoms) {
        mdsys_t sys;
        sys.natoms = 2;
        sys.mass = 39.948; //argon
        sys.vx = new double[2];
        sys.vy = new double[2];
        sys.vz = new double[2];
        sys.vx[0] = 1.5;
        sys.vy[0] = 2.5;
        sys.vz[0] = 3.5;
        sys.vx[1] = 10.5;
        sys.vy[1] = 15.6;
        sys.vz[1] = 20.7;
        ekin(&sys);
	free(sys.vx);
	free(sys.vy);
	free(sys.vz);
        ASSERT_DOUBLE_EQ(sys.ekin, 38327260.75777286);
}

TEST(kinetic_Energy, KE_3atoms) {
	mdsys_t sys;
	sys.natoms = 3;
	sys.mass = 10.0;
	sys.vx = new double[3];
	sys.vy = new double[3];
	sys.vz = new double[3];	
        sys.vx[0] = 1.5;
        sys.vy[0] = 2.5;
        sys.vz[0] = 3.5;	
        sys.vx[1] = 10.5;
        sys.vy[1] = 15.6;
        sys.vz[1] = 20.7;
        sys.vx[2] = 12.5;
        sys.vy[2] = 18.4;
        sys.vz[2] = 24.2;
	ekin(&sys);
	free(sys.vx);
	free(sys.vy);
	free(sys.vz);
        ASSERT_DOUBLE_EQ(sys.ekin, 22505975.144880109);
}
