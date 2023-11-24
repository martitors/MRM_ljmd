#ifndef INPUT_H
#ifdef __cplusplus
extern "C"{
    #endif
        #define INPUT_H
        #include <stdlib.h>
        #include <stdio.h>
        #include "utilities.h" // for the get_a_line function
        void read (char *line, char *restfile, char *trajfile, char *ergfile, mdsys_t *sys, int *nprint);
    #ifdef __cplusplus
   } 
#endif
#endif
