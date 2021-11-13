/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

/********************************************************/
/* ALLOCATE MEMORY FOR THE ENVELOPES AND TEMPORARY DATA */
/********************************************************/

envelope_t init_local_mean(int n) {
  envelope_t env;
  env.e_min = (double*)malloc(n*sizeof(double));
  env.e_max = (double*)malloc(n*sizeof(double));
  env.tmp1 = (double*)malloc(n*sizeof(double));
  env.tmp2 = (double*)malloc(n*sizeof(double));
  return env;
}

/*************************/
/* FREE ALLOCATED MEMORY */
/*************************/

void free_local_mean(envelope_t env) {
  free(env.e_min);
  free(env.e_max);
  free(env.tmp1);
  free(env.tmp2);
}

/*********************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES OF THE CURRENT IMF */
/*********************************************************/

int mean(double *x,double *z,double *m,int n,extrema_t *ex,envelope_t *env) {
  int i;
  /* detect maxima and minima */
  extr(x,z,n,ex);
  /* if not enough extrema -> stop */
  if (ex->n_min+ex->n_max <7)
    return 1;
  /* add extra points at the edges */
  boundary_conditions(x,z,n,ex);
  /* interpolation - upper envelope */
  interpolation(env->e_max,ex->x_max,ex->y_max,ex->n_max,x,n,env->tmp1,env->tmp2);
  /* interpolation - lower envelope */
  interpolation(env->e_min,ex->x_min,ex->y_min,ex->n_min,x,n,env->tmp1,env->tmp2);
  if ((ex->n_min > LIM_GMP)||(ex->n_min > LIM_GMP)) {
    mexWarnMsgTxt("Too many extrema, interpolation may be unreliable\n");
  }
  /* compute the mean */
  for (i=0;i<n;i++) m[i]=(env->e_max[i]+env->e_min[i])/2;
  return 0;
}

/***************************************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES AND THE AMPLITUDE OF THE CURRENT IMF */
/***************************************************************************/

int mean_and_amplitude(double *x,double *z,double *m,double *a,int n,extrema_t *ex,envelope_t *env) {
  int i;
  /* detect maxima and minima */
  extr(x,z,n,ex);
  /* if not enough extrema -> stop */
  if (ex->n_min+ex->n_max <7)
    return 1;
  /* add extra points at the edges */
  boundary_conditions(x,z,n,ex);
  /* interpolation - upper envelope */
  interpolation(env->e_max,ex->x_max,ex->y_max,ex->n_max,x,n,env->tmp1,env->tmp2);
  /* interpolation - lower envelope */
  interpolation(env->e_min,ex->x_min,ex->y_min,ex->n_min,x,n,env->tmp1,env->tmp2);
  /* compute the mean */
  for (i=0;i<n;i++) m[i]=(env->e_max[i]+env->e_min[i])/2;
  /* compute the amplitude */
  for (i=0;i<n;i++) a[i]=(env->e_max[i]-env->e_min[i])/2;
  return 0;
}
