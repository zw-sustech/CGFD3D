/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

double CREAL(complex_data_t c){return c.r;}
double CIMAG(complex_data_t c){return c.i;}
double CABS(complex_data_t c){return sqrt((c.r)*(c.r)+(c.i)*(c.i));}
double crealeiphi(double phi,complex_data_t c){return cos(phi)*c.r-sin(phi)*c.i;}
