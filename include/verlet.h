#ifndef VERLET_H
#define VERLET_H
    
    #include "types.h"
    #include "force_compute.h"
    #include "constants.h"

    void velverlet_1(mdsys_t *sys);
    void velverlet_2(mdsys_t *sys);

#endif