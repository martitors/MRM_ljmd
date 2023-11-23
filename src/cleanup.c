#include "cleanup.h"

void cleanup(mdsys_t *sys, FILE *erg, FILE *traj) {

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

}
