#ifndef UTILITIES_H
#ifdef __cplusplus
extern "C"{
    #endif
    #define UTILITIES_H
        #include <sys/time.h>
        #include <stdio.h>
        #include <stdlib.h>
        #include <string.h>
        #include <ctype.h>
        #include "constants.h"
        #include "types.h"

            double wallclock();
            void azzero(double *d, const int n);
            double pbc(double x, const double boxby2);
            int get_a_line(FILE *fp, char *buf);
            void ekin(mdsys_t *sys);

    #ifdef __cplusplus
   } 
#endif
#endif