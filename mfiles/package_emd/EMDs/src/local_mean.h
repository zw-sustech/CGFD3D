/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

#ifndef LOCAL_MEAN_H
#define LOCAL_MEAN_H

/* structure used to store envelopes and temporary data */
typedef struct {
  int n;
  double *e_min;
  double *e_max;
  double *tmp1;
  double *tmp2;
} envelope_t;

envelope_t init_local_mean(int);
void free_local_mean(envelope_t);
int mean(double *,double *,double *,int,extrema_t *,envelope_t *);
int mean_and_amplitude(double *,double *,double *,double *,int,extrema_t *,envelope_t *);

#endif
