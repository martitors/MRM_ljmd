#include "../include/cleanup.h"

#ifdef _MPI
#include "mpi.h"
#endif //_MPI
void cleanup(mdsys_t sys, FILE *erg, FILE *traj) {

	#ifdef _MPI
	if(sys.rank==0){
	fclose(erg);
    	fclose(traj);
	}
		free(sys.rx);
    	free(sys.ry);
    	free(sys.rz);
    	free(sys.vx);
    	free(sys.vy);
    	free(sys.vz);
    	free(sys.fx);
    	free(sys.fy);
    	free(sys.fz);
	#else
	fclose(erg);
    	fclose(traj);

    	free(sys.rx);
    	free(sys.ry);
    	free(sys.rz);
    	free(sys.vx);
    	free(sys.vy);
    	free(sys.vz);
    	free(sys.fx);
    	free(sys.fy);
    	free(sys.fz);
	#endif
}
