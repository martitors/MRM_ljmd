#ifndef TYPES_H
#ifdef __cplusplus
extern "C" {
	#endif
	
	#define TYPES_H
	
    		/* structure to hold the complete information
   		 * about the MD system */

    		struct _mdsys {
        		int natoms,nfi,nsteps;
        		double dt, mass, epsilon, sigma, box, rcut;
        		double ekin, epot, temp;
        		double *rx, *ry, *rz;
        		double *vx, *vy, *vz;
        		double *fx, *fy, *fz;
				int nthreads;
    		};
    		typedef struct _mdsys mdsys_t;
	#ifdef __cplusplus
	}
	#endif
#endif
