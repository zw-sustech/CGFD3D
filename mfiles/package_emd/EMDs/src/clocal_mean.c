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
  env.re_min = (double*)malloc(n*sizeof(double));
  env.ie_min = (double*)malloc(n*sizeof(double));
  env.re_max = (double*)malloc(n*sizeof(double));
  env.ie_max = (double*)malloc(n*sizeof(double));
  env.tmp1 = (double*)malloc(n*sizeof(double));
  env.tmp2 = (double*)malloc(n*sizeof(double));
  return env;
}

/*************************/
/* FREE ALLOCATED MEMORY */
/*************************/

void free_local_mean(envelope_t env) {
  free(env.re_min);
  free(env.ie_min);
  free(env.re_max);
  free(env.ie_max);
  free(env.tmp1);
  free(env.tmp2);
}

/***************************************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES AND THE AMPLITUDE OF THE CURRENT IMF */
/***************************************************************************/

int mean_and_amplitude(double *x,COMPLEX_T *z,COMPLEX_T *m,double *a,int n,int nbphases,extrema_t *ex,envelope_t *env) {
  int i,k;
  double phi;
  
  #ifdef C99_OK
  for (i=0;i<n;i++) m[i]=0;
  #else
  for (i=0;i<n;i++) {
    m[i].r = 0;
    m[i].i = 0;
  }
  #endif
  for (i=0;i<n;i++) a[i]=0;
  
  for(k=0;k<nbphases;k++) {
    
    phi = k*M_PI/nbphases;
    /* detect maxima and minima in direction phi=k*M_PI/nbphases*/
    extr(x,z,phi,n,ex);
    if (ex->n_max+ex->n_min <3){ /* not enough extrema in a direction -> stop */
      return 1;
    }

    /* add extra points at the edges */
    boundary_conditions(x,z,phi,n,ex);
    
    /* interpolation - upper envelope */
    interpolation(env->re_max,ex->x_max,ex->ry_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    interpolation(env->ie_max,ex->x_max,ex->iy_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    
    /* interpolation - lower envelope */
    interpolation(env->re_min,ex->x_min,ex->ry_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    interpolation(env->ie_min,ex->x_min,ex->iy_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    if ((ex->n_min > LIM_GMP)||(ex->n_min > LIM_GMP)) {
      mexWarnMsgTxt("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* compute the mean and amplitude*/
    #ifdef C99_OK
    for (i=0;i<n;i++) m[i]+=(env->re_max[i]+env->re_min[i]+I*(env->ie_max[i]+env->ie_min[i]))/(2*nbphases);
    for (i=0;i<n;i++) a[i]+=CABS(env->re_max[i]-env->re_min[i]+I*(env->ie_max[i]-env->ie_min[i]))/(2*nbphases);
    #else
    for (i=0;i<n;i++) {
      m[i].r+=(env->re_max[i]+env->re_min[i])/(2*nbphases);
      m[i].i+=(env->ie_max[i]+env->ie_min[i])/(2*nbphases);
    }
    for (i=0;i<n;i++) a[i]+=sqrt((env->re_max[i]-env->re_min[i])*(env->re_max[i]-env->re_min[i])+(env->ie_max[i]-env->ie_min[i])*(env->ie_max[i]-env->ie_min[i]))/(2*nbphases);
    #endif
  }
 
  return 0;
}

/*********************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES OF THE CURRENT IMF */
/*********************************************************/

int mean(double *x,COMPLEX_T *z,COMPLEX_T *m,int n,int nbphases,extrema_t *ex,envelope_t *env) {
  int i,k;
  double phi;
  
  #ifdef C99_OK
  for (i=0;i<n;i++) m[i]=0;
  #else
  for (i=0;i<n;i++) {
    m[i].r = 0;
    m[i].i = 0;
  }
  #endif
  
  for(k=0;k<nbphases;k++) {
    
    phi = k*M_PI/nbphases;
    /* detect maxima and minima in direction phi=k*M_PI/nbphases*/
    extr(x,z,phi,n,ex);
    if (ex->n_max+ex->n_min <3){ /* not enough extrema in a direction -> stop */
      return 1;
    }
    
    boundary_conditions(x,z,phi,n,ex);
    
    /* interpolation - upper envelope */
    if (ex->n_max < LIM_GMP) {
      interpolation(env->re_max,ex->x_max,ex->ry_max,ex->n_max,x,n,env->tmp1,env->tmp2);
      interpolation(env->ie_max,ex->x_max,ex->iy_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    }
    else {
      mexWarnMsgTxt("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* interpolation - lower envelope */
    if (ex->n_min < LIM_GMP) {
      interpolation(env->re_min,ex->x_min,ex->ry_min,ex->n_min,x,n,env->tmp1,env->tmp2);
      interpolation(env->ie_min,ex->x_min,ex->iy_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    }
    else {
      mexWarnMsgTxt("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* compute the mean*/
    #ifdef C99_OK
    for (i=0;i<n;i++) m[i]+=(env->re_max[i]+env->re_min[i]+I*(env->ie_max[i]+env->ie_min[i]))/(2*nbphases);
    #else
    for (i=0;i<n;i++) {
      m[i].r+=(env->re_max[i]+env->re_min[i])/(2*nbphases);
      m[i].i+=(env->ie_max[i]+env->ie_min[i])/(2*nbphases);
    }
    #endif
    
  }
  return 0;
}

