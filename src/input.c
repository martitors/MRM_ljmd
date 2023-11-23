#include "input.h"

void read(char *line, char *restfile, char *trajfile, char *ergfile, mdsys_t *sys, int *nprint)
{
    if (get_a_line(stdin, line))
        exit(1);
    sys->natoms = atoi(line);

    if (get_a_line(stdin, line))
        exit(1);
    sys->mass = atof(line);

    if (get_a_line(stdin, line))
        exit(1);
    sys->epsilon = atof(line);

    if (get_a_line(stdin, line))
        exit(1);
    sys->sigma = atof(line);

    if (get_a_line(stdin, line))
        exit(1);
    sys->rcut = atof(line);

    if (get_a_line(stdin, line))
        exit(1);
    sys->box = atof(line);

    if (get_a_line(stdin, restfile))
        exit(1);
    if (get_a_line(stdin, trajfile))
        exit(1);
    if (get_a_line(stdin, ergfile))
        exit(1);

    if (get_a_line(stdin, line))
        exit(1);
    sys->nsteps = atoi(line);

    if (get_a_line(stdin, line))
        exit(1);
    sys->dt = atof(line);

    if (get_a_line(stdin, line))
        exit(1);
    *nprint = atoi(line);
}