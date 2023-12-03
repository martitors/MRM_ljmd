#ifndef CELL_H
#ifdef __cplusplus
extern "C" {	
    #endif
    #define CELL_H

    #include "types.h"
    #include "constants.h"
    #include "utilities.h"
    #include <math.h>


    void build_cells(mdsys_t *sys);
    void update_cells(mdsys_t *sys);
    void free_cells(mdsys_t *sys);
    void force(mdsys_t *sys);



    #ifdef  __cplusplus
	}
	#endif
#endif
