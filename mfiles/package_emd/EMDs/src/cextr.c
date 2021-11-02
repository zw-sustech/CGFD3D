/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF EXTREMA STRUCTURE                                  */
/*                                                                      */
/************************************************************************/

extrema_t init_extr(int n) {
    extrema_t ex;
    ex.ind_min=(int *)malloc(n*sizeof(int));
    ex.ind_max=(int *)malloc(n*sizeof(int));
    ex.x_min=(double *)malloc(n*sizeof(double));
    ex.x_max=(double *)malloc(n*sizeof(double));
    ex.ry_min=(double *)malloc(n*sizeof(double));
    ex.ry_max=(double *)malloc(n*sizeof(double));
    ex.iy_min=(double *)malloc(n*sizeof(double));
    ex.iy_max=(double *)malloc(n*sizeof(double));
    return ex;
}


/************************************************************************/
/*                                                                      */
/* DETECTION OF LOCAL EXTREMA                                           */
/*                                                                      */
/************************************************************************/

void extr(double x[], COMPLEX_T z[], double phi,int n,extrema_t *ex) {
    int cour;
    double val,valp,valn;
    #ifdef C99_OK
    COMPLEX_T eiphi;
    #endif
    ex->n_min=0;
    ex->n_max=0;
    
    #ifdef C99_OK
    eiphi = cexp(I*phi);
    #endif
    #ifdef C99_OK
    val = CREAL(eiphi*z[0]);
    #else
    val = crealeiphi(phi,z[0]);
    #endif

    #ifdef C99_OK
    valn = CREAL(eiphi*z[1]);
    #else
    valn = crealeiphi(phi,z[1]);
    #endif

    /* search for extrema in direction phi*/
    for(cour=1;cour<(n-1);cour++) {
        valp = val;
        val = valn;
        #ifdef C99_OK
        valn = CREAL(eiphi*z[cour+1]);
        #else
        valn = crealeiphi(phi,z[cour+1]);
        #endif

        if (val<valp && val<valn) /* local minimum */ {
            ex->x_min[ex->n_min+NBSYM]=x[cour];
            ex->ry_min[ex->n_min+NBSYM]=CREAL(z[cour]);
            ex->iy_min[ex->n_min+NBSYM]=CIMAG(z[cour]);
            ex->ind_min[ex->n_min+NBSYM]=cour;
            ex->n_min++;
        }
        if (val>valp && val>valn) /* local maximum */ {
            ex->x_max[ex->n_max+NBSYM]=x[cour];
            ex->ry_max[ex->n_max+NBSYM]=CREAL(z[cour]);
            ex->iy_max[ex->n_max+NBSYM]=CIMAG(z[cour]);
            ex->ind_max[ex->n_max+NBSYM]=cour;
            ex->n_max++;
        }
    }
}

/************************************************************************/
/*                                                                      */
/* EXTRAPOLATION OF EXTREMA TO LIMIT BORDER EFFECTS                     */
/*                                                                      */
/************************************************************************/

void boundary_conditions(double x[],COMPLEX_T z[],double phi,int n,extrema_t *ex) {
    int cour,nbsym;
    #ifdef C99_OK
    COMPLEX_T eiphi;
    #endif
    #ifdef C99_OK
    eiphi = cexp(I*phi);
    #endif
    
    nbsym = NBSYM;
    /* reduce the number of symmetrized points if there is not enough extrema */
    while(ex->n_min < nbsym+1 && ex->n_max < nbsym+1) nbsym--;
    if (nbsym < NBSYM) {
        for(cour=0;cour<ex->n_max;cour++) {
            ex->ind_max[nbsym+cour] = ex->ind_max[NBSYM+cour];
            ex->x_max[nbsym+cour] = ex->x_max[NBSYM+cour];
            ex->ry_max[nbsym+cour] = ex->ry_max[NBSYM+cour];
            ex->iy_max[nbsym+cour] = ex->iy_max[NBSYM+cour];
        }
        for(cour=0;cour<ex->n_min;cour++) {
            ex->ind_min[nbsym+cour] = ex->ind_min[NBSYM+cour];
            ex->x_min[nbsym+cour] = ex->x_min[NBSYM+cour];
            ex->ry_min[nbsym+cour] = ex->ry_min[NBSYM+cour];
            ex->iy_min[nbsym+cour] = ex->iy_min[NBSYM+cour];
        }
    }
    
    /* select the symmetrized points and the axis of symmetry at the beginning of the signal*/
    if (ex->x_max[nbsym] < ex->x_min[nbsym]) { /* first = max */
        #ifdef C99_OK
        if (CREAL(eiphi*z[0]) > CREAL(eiphi*z[ex->ind_min[nbsym]])) { /* the edge is not a min */
        #else
        if (crealeiphi(phi,z[0]) > crealeiphi(phi,z[ex->ind_min[nbsym]])) { /* the edge is not a min */
        #endif

            if (2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*ex->x_max[nbsym]-ex->x_max[2*nbsym-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-cour];
                    ex->x_min[cour] = 2*ex->x_max[nbsym]-ex->x_min[2*nbsym-1-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
                }
            }
        } else { /* edge is a min -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-2-cour];
                ex->ry_min[cour] = ex->ry_min[2*nbsym-2-cour];
                ex->iy_min[cour] = ex->iy_min[2*nbsym-2-cour];
            }
            ex->x_min[nbsym-1] = x[0];
            ex->ry_min[nbsym-1] = CREAL(z[0]);
            ex->iy_min[nbsym-1] = CIMAG(z[0]);
        }
    } else { /* first = min */
        #ifdef C99_OK
        if (CREAL(eiphi*z[0]) < CREAL(eiphi*z[ex->ind_max[nbsym]])) { /* the edge is not a max */
        #else
        if (crealeiphi(phi,z[0]) < crealeiphi(phi,z[ex->ind_max[nbsym]])) { /* the edge is not a max */
        #endif

            if (2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1] > x[0]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-1-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[cour] = 2*ex->x_min[nbsym]-ex->x_max[2*nbsym-1-cour];
                    ex->ry_max[cour] = ex->ry_max[2*nbsym-1-cour];
                    ex->iy_max[cour] = ex->iy_max[2*nbsym-1-cour];
                    ex->x_min[cour] = 2*ex->x_min[nbsym]-ex->x_min[2*nbsym-cour];
                    ex->ry_min[cour] = ex->ry_min[2*nbsym-cour];
                    ex->iy_min[cour] = ex->iy_min[2*nbsym-cour];
                }
            }
        } else { /* edge is a max -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_min[cour] = 2*x[0]-ex->x_min[2*nbsym-1-cour];
                ex->ry_min[cour] = ex->ry_min[2*nbsym-1-cour];
                ex->iy_min[cour] = ex->iy_min[2*nbsym-1-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_max[cour] = 2*x[0]-ex->x_max[2*nbsym-2-cour];
                ex->ry_max[cour] = ex->ry_max[2*nbsym-2-cour];
                ex->iy_max[cour] = ex->iy_max[2*nbsym-2-cour];
            }
            ex->x_max[nbsym-1] = x[0];
            ex->ry_max[nbsym-1] = CREAL(z[0]);
            ex->iy_max[nbsym-1] = CIMAG(z[0]);
        }
    }
    
    
    (ex->n_min) += nbsym-1;
    (ex->n_max) += nbsym-1;
    
    
    /* select the symmetrized points and the axis of symmetry at the end of the signal*/
    if (ex->x_max[ex->n_max] < ex->x_min[ex->n_min]) { /* last is a min */
        #ifdef C99_OK
        if (CREAL(eiphi*z[n-1]) < CREAL(eiphi*z[ex->ind_max[ex->n_max]])) { /* the edge is not a max */
        #else
        if (crealeiphi(phi,z[n-1]) < crealeiphi(phi,z[ex->ind_max[ex->n_max]])) { /* the edge is not a max */
        #endif

            if (2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_max[ex->n_max-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*ex->x_min[ex->n_min]-ex->x_min[ex->n_min-1-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-1-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-1-cour];
                }
            }
        } else { /* edge is a max -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_max[ex->n_max+2+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                ex->ry_max[ex->n_max+2+cour] = ex->ry_max[ex->n_max-cour];
                ex->iy_max[ex->n_max+2+cour] = ex->iy_max[ex->n_max-cour];
            }
            ex->x_max[ex->n_max+1] = x[n-1];
            ex->ry_max[ex->n_max+1] = CREAL(z[n-1]);
            ex->iy_max[ex->n_max+1] = CIMAG(z[n-1]);
        }
    } else {  /* last is a max */
        #ifdef C99_OK
        if (CREAL(eiphi*z[n-1]) > CREAL(eiphi*z[ex->ind_min[ex->n_min]])) { /* the edge is not a min */
        #else
        if (crealeiphi(phi,z[n-1]) > crealeiphi(phi,z[ex->ind_min[ex->n_min]])) { /* the edge is not a min */
        #endif

            if (2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-nbsym+1] < x[n-1]) { /* symmetrized parts are too short */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
                }
            } else { /* symmetrized parts are long enough */
                for(cour=0;cour<nbsym;cour++) {
                    ex->x_max[ex->n_max+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_max[ex->n_max-1-cour];
                    ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-1-cour];
                    ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-1-cour];
                    ex->x_min[ex->n_min+1+cour] = 2*ex->x_max[ex->n_max]-ex->x_min[ex->n_min-cour];
                    ex->ry_min[ex->n_min+1+cour] = ex->ry_min[ex->n_min-cour];
                    ex->iy_min[ex->n_min+1+cour] = ex->iy_min[ex->n_min-cour];
                }
            }
        } else { /* edge is a min -> sym with respect to the edge*/
            for(cour=0;cour<nbsym;cour++) {
                ex->x_max[ex->n_max+1+cour] = 2*x[n-1]-ex->x_max[ex->n_max-cour];
                ex->ry_max[ex->n_max+1+cour] = ex->ry_max[ex->n_max-cour];
                ex->iy_max[ex->n_max+1+cour] = ex->iy_max[ex->n_max-cour];
            }
            for(cour=0;cour<nbsym-1;cour++) {
                ex->x_min[ex->n_min+2+cour] = 2*x[n-1]-ex->x_min[ex->n_min-cour];
                ex->ry_min[ex->n_min+2+cour] = ex->ry_min[ex->n_min-cour];
                ex->iy_min[ex->n_min+2+cour] = ex->iy_min[ex->n_min-cour];
            }
            ex->x_min[ex->n_min+1] = x[n-1];
            ex->ry_min[ex->n_min+1] = CREAL(z[n-1]);
            ex->iy_min[ex->n_min+1] = CIMAG(z[n-1]);
        }
    }
    
    
    
    (ex->n_min) = ex->n_min + nbsym + 1;
    (ex->n_max) = ex->n_max + nbsym + 1;
}


/************************************************************************/
/*                                                                      */
/* FREE ALLOCATED MEMORY                                                */
/*                                                                      */
/************************************************************************/

void free_extr(extrema_t ex) {
    free(ex.x_max);
    free(ex.x_min);
    free(ex.ind_max);
    free(ex.ind_min);
    free(ex.ry_max);
    free(ex.ry_min);
    free(ex.iy_max);
    free(ex.iy_min);
}
