#ifndef FORCE_COMPUTE_H
#ifdef __cplusplus
extern "C" {	
    #endif
    #define FORCE_COMPUTE_H
    #include <math.h>

    #include "types.h"
    #include "constants.h"
    #include "utilities.h"
    #include <omp.h>

    void force(mdsys_t *sys);
    void ekin(mdsys_t *sys);

    #ifdef  __cplusplus
	}
	#endif
#endif
