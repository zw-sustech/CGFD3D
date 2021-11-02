/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

#ifndef CEXTR_H
#define CEXTR_H

/* structure used to store the local extrema of the signal */
typedef struct {
  int n_min;
  int n_max;
  int *ind_min;
  int *ind_max;
  double *x_min;
  double *ry_min;
  double *iy_min;
  double *x_max;
  double *ry_max;
  double *iy_max;
} extrema_t;

extrema_t init_extr(int);
void extr(double *,COMPLEX_T *,double,int,extrema_t *);
void free_extr(extrema_t);
void boundary_conditions(double *,COMPLEX_T *,double,int,extrema_t *);

#endif
