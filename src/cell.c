#include "../include/cell.h"

void build_cells(mdsys_t *sys){

    /*ratio between cellsize and rcut*/
    
    double cell_rcut_ratio = 2.0;
    double rsq;

    double cell_size = cell_rcut_ratio * sys->rcut;

    if (cell_size <= 0.0) {
        fprintf(stderr, "Error: cell_size is zero or negative.\n");
        exit(1);
    }

    double c = floor(sys->box / cell_size);

    if (c <= 0.0) {
        fprintf(stderr, "Error: 'number of cells per side' is zero or negative.\n");
        exit(1);
    }

    int cells_per_side = c;
    
    sys->cell_size = sys->box / cells_per_side;

    double rcsq = (sys->rcut + sqrt(3.0) * sys->cell_size) * (sys->rcut + sqrt(3.0) * sys->cell_size);

    int ncells = cells_per_side * cells_per_side * cells_per_side;

    sys->clist = (cell_t *)malloc(ncells * sizeof(cell_t));

    int avg_density = sys->natoms / ncells;

    int idx;
    double rx,ry,rz;
    for (int z=0; z < cells_per_side; ++z) {
         for (int y=0; y < cells_per_side; ++y) {
             for (int x=0; x < cells_per_side; ++x) {
        
        idx = z * (cells_per_side * cells_per_side) + y * (cells_per_side) + x;
        sys->clist[idx].cell = (int *) malloc( (3 * avg_density + 2)* sizeof(int));
        sys->clist[idx].atoms = 0 ;
        sys->clist[idx].center_x = (x * sys->cell_size + sys->cell_size/2.0) - 0.5 * (sys->box -  sys->cell_size/2.0) ;
        sys->clist[idx].center_y = (y * sys->cell_size + sys->cell_size/2.0) - 0.5 * (sys->box -  sys->cell_size/2.0) ;
        sys->clist[idx].center_z = (z * sys->cell_size + sys->cell_size/2.0) - 0.5 * (sys->box -  sys->cell_size/2.0) ;
             }
         }
    }

    /*Allocate and fill the pair list*/
    sys->plist = (int*) malloc( 2*ncells*ncells*sizeof(int));
    int npair = 0;
    for (int i=0; i < ncells-1; ++i) {
        for (int j=i+1; j < ncells; ++j) {
        
        rx=pbc(sys->clist[i].center_x - sys->clist[j].center_x, 0.5* sys->box);
        ry=pbc(sys->clist[i].center_y - sys->clist[j].center_y, 0.5* sys->box);
        rz=pbc(sys->clist[i].center_z - sys->clist[j].center_z, 0.5* sys->box);

        rsq = rx*rx + ry*ry + rz*rz;

        if (rsq >rcsq) continue;

        sys->plist[2 * npair] = i;
        sys->plist[2 * npair + 1] = j;
        ++npair;
        }
    }
    sys->npair = npair;
    sys->ncells = ncells;
    sys->cells_per_side= cells_per_side;

}

void update_cells(mdsys_t *sys){

    for (int j = 0;j < sys->ncells;j++) {
        sys->clist[j].atoms = 0;
        }

    for (int i=0; i < sys->natoms; ++i) {

        int x = 0;
        int y = 0;
        int z = 0;
        double rx,ry,rz;
        
        rx= pbc(sys->rx[i], 0.5* sys->box);
        ry= pbc(sys->ry[i], 0.5* sys->box);
        rz= pbc(sys->rz[i], 0.5* sys->box);

        while (rx > (sys->clist[x].center_x + sys->cell_size/2.0)) {
             ++x;
            };

        while (ry > (sys->clist[y].center_y + sys->cell_size/2.0)) {
            y += (sys->cells_per_side);
            };

        while (rz > (sys->clist[z].center_z + sys->cell_size/2.0)) {
            z += (sys->cells_per_side * sys->cells_per_side) ;
            };

        int idx = z + y + x;

        sys->clist[idx].cell[sys->clist[idx].atoms] = i;
        sys->clist[idx].atoms += 1;

        if (sys->clist[idx].atoms > (3 * sys->natoms/ sys->ncells + 2 )) {
            printf("Cell overflow!");
            exit(1);
        }
    }

}

void free_cells(mdsys_t *sys)
{
    int i;
    if (sys->clist == NULL) 
        return;
    
    for (i=0; i < sys->ncells; ++i) {
        free(sys->clist[i].cell);
    }
    free(sys->clist);
    sys->clist = NULL;
    sys->ncells = 0;
}

/* compute forces */
void force(mdsys_t *sys) 
{
    
    sys->epot=0.0;

    double rx,ry,rz,rsq;
    double rx1, ry1, rz1;
    double c12,c6,boxby2,rcsq;
    int i,natoms;

    /* precompute some constants */
    c12 = 4.0*sys->epsilon*pow(sys->sigma,12.0);
    c6  = 4.0*sys->epsilon*pow(sys->sigma, 6.0);
    rcsq= sys->rcut * sys->rcut;
    boxby2 = 0.5*sys->box;
    natoms = sys->natoms;

    azzero(sys->fx,sys->natoms);
    azzero(sys->fy,sys->natoms);
    azzero(sys->fz,sys->natoms);

        /* self interaction of atoms in cell */
        for(i=0; i < sys->ncells; i += 1) {
            cell_t *c1;

            if (i >= (sys->ncells)) break;
            c1=sys->clist + i;

            for (int j=0; j < c1->atoms-1; ++j) {

                int ii=c1->cell[j];
                rx1=sys->rx[ii];
                ry1=sys->ry[ii];
                rz1=sys->rz[ii];
        
                for(int k=j+1; k < c1->atoms; ++k) {

                    int jj=c1->cell[k];

                    /* get distance between particle i and j */
                    rx = pbc(rx1 - sys->rx[jj], 0.5 * sys->box);
                    ry = pbc(ry1 - sys->ry[jj], 0.5 * sys->box);
                    rz = pbc(rz1 - sys->rz[jj], 0.5 * sys->box);
                    
                    rsq = rx*rx + ry*ry + rz*rz;

                    /* compute force and energy if within cutoff */
                    if (rsq < rcsq) {
                        double r6,rinv,ffac;

                        rinv=1.0/rsq;
                        r6=rinv*rinv*rinv;
                    
                        ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                        sys->epot += r6*(c12*r6 - c6);

                        sys->fx[ii] += rx*ffac;         sys->fx[jj] -= rx*ffac;
                        sys->fy[ii] += ry*ffac;         sys->fy[jj] -= ry*ffac;
                        sys->fz[ii] += rz*ffac;         sys->fz[jj] -= rz*ffac;
                        
                    }
                }
            }
        }    

        /* interaction of atoms in different cells */
        for(i=0; i < sys->npair; i += sys->nthreads) {

            cell_t *c1, *c2;
            if (i >= (sys->npair)) break;

            c1=sys->clist + sys->plist[2*i];
            c2=sys->clist + sys->plist[2*i+1];
        
            for (int j=0; j < c1->atoms; ++j) {

                int ii=c1->cell[j];
                rx1=sys->rx[ii];
                ry1=sys->ry[ii];
                rz1=sys->rz[ii];
        
                for(int k=0; k < c2->atoms; ++k) {
                    int jj;
                
                    jj=c2->cell[k];
                
                    /* get distance between particle i and j */
                    rx=pbc(rx1 - sys->rx[jj], 0.5 * sys->box);
                    ry=pbc(ry1 - sys->ry[jj], 0.5 * sys->box);
                    rz=pbc(rz1 - sys->rz[jj], 0.5 * sys->box);
                    rsq = rx*rx + ry*ry + rz*rz;
                
                    /* compute force and energy if within cutoff */
                    if (rsq < rcsq) {
                        double r6,rinv,ffac;

                        rinv=1.0/rsq;
                        r6=rinv*rinv*rinv;
                    
                        ffac = (12.0*c12*r6 - 6.0*c6)*r6*rinv;
                        sys->epot += r6*(c12*r6 - c6);

                        sys->fx[ii] += rx*ffac;         sys->fx[jj] -= rx*ffac;
                        sys->fy[ii] += ry*ffac;         sys->fy[jj] -= ry*ffac;
                        sys->fz[ii] += rz*ffac;         sys->fz[jj] -= rz*ffac;
                    }
                }
            }
        }
}

