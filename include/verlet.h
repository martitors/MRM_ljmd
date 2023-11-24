#ifndef VERLET_H
#ifdef __cplusplus
extern "C" {	
	#endif

    
    #include "types.h"
    #include "force_compute.h"
    #include "constants.h"

    void velverlet_1(mdsys_t *sys);
    void velverlet_2(mdsys_t *sys);


	#ifdef  __cplusplus
	}
	#endif
#endif
