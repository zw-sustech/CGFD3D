/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

#ifndef EMD_IO_H
#define EMD_IO_H

#define GREATER(A,B) ((A) >= (B) ? (A) : (B))
#define SMALLER(A,B) ((A) <  (B) ? (A) : (B))

/* structure used to store input data */
typedef struct {
  int n;
  int max_imfs;
  int allocated_x;
  #ifdef _ALT_MEXERRMSGTXT_
  int error_flag;
  #endif
  double *x;
  double *y;
  int nb_iterations;
} input_t;


/* structure used to store an IMF and the associated number of iterations */
typedef struct i {
  int nb_iterations;
  double *pointer;
  struct i *next;
} imf_t;

/* structure of the IMF list */
typedef struct {
  imf_t *first;
  imf_t *last;
  int m; 
  int n; 
} imf_list_t;

input_t get_input(int,int,const mxArray **);
imf_list_t init_imf_list(int);
void add_imf(imf_list_t *,double *,int);
void free_imf_list(imf_list_t);
void write_output(imf_list_t,mxArray **);

#endif
