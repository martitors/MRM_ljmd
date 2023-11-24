#ifndef OUTPUT_H
#ifdef __cplusplus
extern "C"{
    #endif
    #define OUTPUT_H
        #include <stdlib.h>
        #include <stdio.h>
        #include "types.h"
        void output(mdsys_t *sys, FILE *erg, FILE *traj);
    #ifdef __cplusplus
   } 
#endif
#endif