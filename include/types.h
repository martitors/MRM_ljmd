#ifndef TYPES_H
#ifdef __cplusplus
extern "C" {
	#endif
	#define TYPES_H
	
    		/* structure to hold the complete information
   		 * about the MD system */
			struct _cell {
				int atoms;                     
				int *cell;  
				double center_x,center_y,center_z;
			};        
			typedef struct _cell cell_t;

    		struct _mdsys {
        		int natoms,nfi,nsteps;
        		double dt, mass, epsilon, sigma, box, rcut;
        		double ekin, epot, temp;
        		double *rx, *ry, *rz;
        		double *vx, *vy, *vz;
        		double *fx, *fy, *fz;
				int nthreads;
				int rank; 
				int npes; 
				double *cx, *cy, *cz; 

				cell_t *clist ;
				int *plist;
				double cell_size;
				int ncells,npair,cells_per_side;
    		};
    		typedef struct _mdsys mdsys_t;
	#ifdef __cplusplus
	}
	#endif
#endif
