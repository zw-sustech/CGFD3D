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
/* GET INPUT DATA                                                       */
/*                                                                      */
/************************************************************************/

input_t get_input(int nlhs,int nrhs,const mxArray *prhs[]) {
  input_t input;
  int n,i;
  double *x,*y,*y_temp,third,fourth;
  
  input.nb_iterations=DEFAULT_NB_ITERATIONS;
  input.allocated_x=0;
  #ifdef _ALT_MEXERRMSGTXT_
  input.error_flag=0;
  #endif
  input.max_imfs=0;

  
  /* argument checking*/
  if (nrhs>4)
    mexErrMsgTxt("Too many arguments");
  if (nrhs<2)
    mexErrMsgTxt("Not enough arguments");
  if (nlhs>2)
    mexErrMsgTxt("Too many output arguments");
  if (!mxIsEmpty(prhs[0]))
    if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) ||
    mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]) ||
    (mxGetNumberOfDimensions(prhs[0]) > 2))
      mexErrMsgTxt("X must be either empty or a double precision real vector.");
  
  if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) ||
  mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1])  ||
  (mxGetNumberOfDimensions(prhs[1]) > 2))
    mexErrMsgTxt("Y must be a double precision real vector.");
  
  /* input reading: x and y */
  n=GREATER(mxGetN(prhs[1]),mxGetM(prhs[1])); /* length of vector x */
  if (mxIsEmpty(prhs[0])) {
    input.allocated_x = 1;
    x = (double *)malloc(n*sizeof(double));
    for(i=0;i<n;i++) x[i] = i;
  }
  else
    x=mxGetPr(prhs[0]);
  y_temp=mxGetPr(prhs[1]);
  
  /* third argument */
  if (nrhs>=3) {
    if (!mxIsEmpty(prhs[2])) { /* if empty -> do nothing */
      if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2])
      || !mxIsDouble(prhs[2]) || mxGetN(prhs[2])!=1 || mxGetM(prhs[2])!=1)
        mexErrMsgTxt("The number of iterations must be a positive integer");
      third=*mxGetPr(prhs[2]);
      if((int)third != third)
        mexErrMsgTxt("The number of iterations must be a positive integer");
      input.nb_iterations = (int)third;
    }
  }
  
  
    /* fourth argument */
  if (nrhs==4) {
    if (!mxIsEmpty(prhs[3])) { /* if empty -> do nothing */
      if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3])
      || !mxIsDouble(prhs[3]) || mxGetN(prhs[3])!=1 || mxGetM(prhs[3])!=1)
        mexErrMsgTxt("NB_IMFS must be a positive integer");
      fourth=*mxGetPr(prhs[3]);
      if ((unsigned int)fourth != fourth)
        mexErrMsgTxt("NB_IMFS must be a positive integer");
      input.max_imfs=(int)fourth;
    }
  }
  
  /* more input checking */
  if (!input.allocated_x && SMALLER(mxGetN(prhs[0]),mxGetM(prhs[0]))!=1 ||
  SMALLER(mxGetN(prhs[1]),mxGetM(prhs[1]))!=1)
    mexErrMsgTxt("X and Y must be vectors");
  if (GREATER(mxGetN(prhs[1]),mxGetM(prhs[1]))!=n)
    mexErrMsgTxt("X and Y must have the same length");
  i=1;
  while (i<n && x[i]>x[i-1]) i++;
  if (i<n) mexErrMsgTxt("Values in X must be non decreasing");
  
  /* copy vector y to avoid erasing input data */
  y=(double *)malloc(n*sizeof(double));
  for (i=0;i<n;i++) y[i]=y_temp[i];
  
  input.n=n;
  input.x=x;
  input.y=y;
  return input;
}


/************************************************************************/
/*                                                                      */
/* INITIALIZATION OF THE LIST                                           */
/*                                                                      */
/************************************************************************/

imf_list_t init_imf_list(int n) {
  imf_list_t list;
  list.first=NULL;
  list.last=NULL;
  list.n=n;
  list.m=0;
  return list;
}


/************************************************************************/
/*                                                                      */
/* ADD AN IMF TO THE LIST                                               */
/*                                                                      */
/************************************************************************/

void add_imf(imf_list_t *list,double *p,int nb_it) {
  double *v=(double *)malloc(list->n*sizeof(double));
  int i;
  imf_t *mode=(imf_t *)malloc(sizeof(imf_t));
  for (i=0;i<list->n;i++) v[i]=p[i];
  mode->pointer=v;
  mode->nb_iterations=nb_it;
  mode->next=NULL;
  if (!list->first) {
    list->first=mode;
  } else {
    (list->last)->next=mode;
  }
  list->last=mode;
  list->m++;
}


/************************************************************************/
/*                                                                      */
/* FREE MEMORY ALLOCATED FOR THE LIST                                   */
/*                                                                      */
/************************************************************************/

void free_imf_list(imf_list_t list) {
  imf_t *current=list.first, *previous;
  while (current) {
    previous=current;
    current=current->next;
    free(previous->pointer);
    free(previous);
  }
}


/************************************************************************/
/*                                                                      */
/* OUTPUT INTO MATLAB ARRAY                                             */
/*                                                                      */
/************************************************************************/

void write_output(imf_list_t list,mxArray *plhs[]) {
  double *out1,*out2;
  imf_t *current;
  int i=0,j,m=list.m,n=list.n;
  plhs[0]=mxCreateDoubleMatrix(m,n,mxREAL);
  out1=mxGetPr(plhs[0]);
  plhs[1]=mxCreateDoubleMatrix(1,m-1,mxREAL);
  out2=mxGetPr(plhs[1]);
  for (current=list.first;current;current=current->next) {
    for (j=0;j<n;j++) *(out1+j*m+i)=current->pointer[j];
    if (i<m-1) *(out2+i)=current->nb_iterations;
    i++;
  }
}
