
#include "../include/utilities.h"
/* helper function: get current time in seconds since epoch */
double wallclock()
{
        struct timeval t;
        gettimeofday(&t,0);
        return ((double) t.tv_sec) + 1.0e-6*((double) t.tv_usec);
}

/* helper function: zero out an array */
void azzero(double *d, const int n)
{
    int i;
    for (i=0; i<n; ++i) {
        d[i]=0.0;
    }
}

/* helper function: apply minimum image convention */
double pbc(double x, const double boxby2)
{
    while (x >  boxby2) x -= 2.0*boxby2;
    while (x < -boxby2) x += 2.0*boxby2;
    return x;
}

int get_a_line(FILE *fp, char *buf)
{
    char tmp[BLEN], *ptr;

    /* read a line and cut of comments and blanks */
    if (fgets(tmp,BLEN,fp)) {
        int i;

        ptr=strchr(tmp,'#');
        if (ptr) *ptr= '\0';
        i=strlen(tmp); --i;
        while(isspace(tmp[i])) {
            tmp[i]='\0';
            --i;
        }
        ptr=tmp;
        while(isspace(*ptr)) {++ptr;}
        i=strlen(ptr);
        strcpy(buf,tmp);
        return 0;
    } else {
        perror("problem reading input");
        return -1;
    }
    return 0;
}

void ekin(mdsys_t *sys)
{
    int i;

    sys->ekin=0.0;
    double ekin_tmp = 0.0;

//#ifdef _OPENMP
//#pragma omp parallel for reduction(+ : ekin_tmp)
//#endif
    for (i=0; i<sys->natoms; ++i) {
        //#ifdef _OPENMP
        //ekin_tmp += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
        //#else
        sys->ekin += 0.5*mvsq2e*sys->mass*(sys->vx[i]*sys->vx[i] + sys->vy[i]*sys->vy[i] + sys->vz[i]*sys->vz[i]);
        //#endif
    }
    //#ifdef _OPENMP
    //sys->ekin = ekin_tmp;
    //#endif
    sys->temp = 2.0*sys->ekin/(3.0*sys->natoms-3.0)/kboltz;
}


static void update_cells(mdsys_t *sys)
{
        
    if (sys->clist == NULL) {
        
        int nidx;
        
        int cell_per_dim = floor( sys->box / sys->rcut);
        sys->cell_size  = sys->box/ cell_per_dim;
        sys->ncells = cell_per_dim * cell_per_dim * cell_per_dim;
        int boxoffs = 0.5 * ( sys->box - sys->cell_size );

        sys->clist = (cell_t *) malloc(sys->ncells*sizeof(cell_t));
        sys->plist = (int *) malloc(2*sys->ncells*sys->ncells*sizeof(int));

    }

}