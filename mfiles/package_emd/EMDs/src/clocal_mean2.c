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

/***************************************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES AND THE AMPLITUDE OF THE CURRENT IMF */
/***************************************************************************/

int mean_and_amplitude(double *x,COMPLEX_T *z,COMPLEX_T *m,double *a,int n,int nbphases,extrema_t *ex,envelope_t *env) {
  int i,k;
  #ifdef C99_OK
  COMPLEX_T eiphi;
  #else
  double phi,cphi,sphi;
  #endif
  
  #ifdef C99_OK
  for (i=0;i<n;i++) m[i]=0;
  #else
  for (i=0;i<n;i++) {
    m[i].r=0;
    m[i].i=0;
  }
  #endif
  for (i=0;i<n;i++) a[i]=0;
  
  for(k=0;k<nbphases;k++) {
    
    #ifdef C99_OK
    eiphi = cexp(-I*k*M_PI/(nbphases));
    for(i=0;i<n;i++) env->tmp1[i] = CREAL(eiphi*z[i]);
    #else
    phi = k*M_PI/(nbphases);
    for(i=0;i<n;i++) env->tmp1[i] = crealeiphi(phi,z[i]);
    #endif

    /* detect maxima and minima in direction phi=k*M_PI/nbphases*/
    extr(x,env->tmp1,n,ex);
    if (ex->n_max+ex->n_min <7){ /* not enough extrema in a direction -> stop */
      return 1;
    }
    
    /* add extra points at the edges */
    boundary_conditions(x,env->tmp1,n,ex);
    
    /* interpolation - upper envelope */
    interpolation(env->e_max,ex->x_max,ex->y_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    
    /* interpolation - lower envelope */
    interpolation(env->e_min,ex->x_min,ex->y_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    
    if ((ex->n_min > LIM_GMP)||(ex->n_min > LIM_GMP)) {
      mexWarnMsgTxt("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* compute the mean and amplitude*/
    #ifdef C99_OK
    eiphi=conj(eiphi);
    for (i=0;i<n;i++) m[i]+=eiphi*(env->e_max[i]+env->e_min[i])/(nbphases);
    #else
    cphi = cos(-phi);
    sphi = sin(-phi);
    for (i=0;i<n;i++) {
      m[i].r+=cphi*(env->e_max[i]+env->e_min[i])/(nbphases);
      m[i].i+=sphi*(env->e_max[i]+env->e_min[i])/(nbphases);
    }
    #endif
    for (i=0;i<n;i++) a[i]+=(env->e_max[i]-env->e_min[i])/(nbphases);
    
  }
  return 0;
}

/*********************************************************/
/* COMPUTES THE MEAN OF THE ENVELOPES OF THE CURRENT IMF */
/*********************************************************/

int mean(double *x,COMPLEX_T *z,COMPLEX_T *m,int n,int nbphases,extrema_t *ex,envelope_t *env) {
  int i,k;
  #ifdef C99_OK
  COMPLEX_T eiphi;
  #else
  double phi,cphi,sphi;
  #endif
  
  #ifdef C99_OK
  for (i=0;i<n;i++) m[i]=0;
  #else
  for (i=0;i<n;i++) {
    m[i].r=0;
    m[i].i=0;
  }
  #endif
  
  for(k=0;k<nbphases;k++) {
    
    #ifdef C99_OK
    eiphi = cexp(-I*k*M_PI/(nbphases));
    for(i=0;i<n;i++) env->tmp1[i] = CREAL(eiphi*z[i]);
    #else
    phi = k*M_PI/(nbphases);
    for(i=0;i<n;i++) env->tmp1[i] = crealeiphi(phi,z[i]);
    #endif

    /* detect maxima and minima in direction phi=k*M_PI/nbphases*/
    extr(x,env->tmp1,n,ex);
    if (ex->n_max+ex->n_min <7){ /* not enough extrema in a direction -> stop */
      return 1;
    }
    
    /* add extra points at the edges */
    boundary_conditions(x,env->tmp1,n,ex);
    
    /* interpolation - upper envelope */
    interpolation(env->e_max,ex->x_max,ex->y_max,ex->n_max,x,n,env->tmp1,env->tmp2);
    
    /* interpolation - lower envelope */
    interpolation(env->e_min,ex->x_min,ex->y_min,ex->n_min,x,n,env->tmp1,env->tmp2);
    
    if ((ex->n_min > LIM_GMP)||(ex->n_min > LIM_GMP)) {
      mexWarnMsgTxt("Too many extrema, interpolation may be unreliable\n");
    }
    
    /* compute the mean*/
    #ifdef C99_OK
    eiphi=conj(eiphi);
    for (i=0;i<n;i++) m[i]+=eiphi*(env->e_max[i]+env->e_min[i])/(nbphases);
    #else
    cphi = cos(-phi);
    sphi = sin(-phi);
    for (i=0;i<n;i++) {
      m[i].r+=cphi*(env->e_max[i]+env->e_min[i])/(nbphases);
      m[i].i+=sphi*(env->e_max[i]+env->e_min[i])/(nbphases);
    }
    #endif
    
  }
  return 0;
}

