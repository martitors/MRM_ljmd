#ifndef UTILITIES_H
#define UTILITIES_H
#include "types.h"
#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "constants.h"

    double wallclock();
    void azzero(double *d, const int n);
    double pbc(double x, const double boxby2);
    int get_a_line(FILE *fp, char *buf);

#endif