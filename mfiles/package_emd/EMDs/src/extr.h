/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

#ifndef EXTR_H
#define EXTR_H

/* structure used to store the local extrema of the signal */
typedef struct {
  int n_min;
  int n_max;
  double *x_min;
  double *y_min;
  double *x_max;
  double *y_max;
} extrema_t;

extrema_t init_extr(int);
void extr(double *,double *,int,extrema_t *);
void free_extr(extrema_t);
void boundary_conditions(double *,double *,int,extrema_t *);

#endif
