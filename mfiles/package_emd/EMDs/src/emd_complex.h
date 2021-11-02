/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

#ifndef EMD_COMPLEX_H
#define EMD_COMPLEX_H
typedef struct {
  double r;
  double i;
} complex_data_t;

double CREAL(complex_data_t);
double CIMAG(complex_data_t);
double CABS(complex_data_t);
double crealeiphi(double,complex_data_t);

#endif
