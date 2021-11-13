/*
* G. Rilling, last modification: 3.2007
* gabriel.rilling@ens-lyon.fr
*
* code based on a student project by T. Boustane and G. Quellec, 11.03.2004
* supervised by P. Chainais (ISIMA - LIMOS - Universite Blaise Pascal - Clermont II
* email : pchainai@isima.fr).
*/

/*************************************************************************/
/*                                                                       */
/* INTERPOLATION                                                         */
/*                                                                       */
/* interpolates the sequence (xs,ys) at instants in x using cubic spline */
/*                                                                       */
/*************************************************************************/

void interpolation(double y[],double xs[],double ys[],int n,double x[], int nx,double *ys2, double *temp) {
  int i,j,jfin,cur,prev;
  double p,sig,a,b,c,d,e,f,g,a0,a1,a2,a3,xc;

  /* Compute second derivatives at the knots */
  ys2[0]=temp[0]=0.0;
  for (i=1;i<n-1;i++) {
    sig=(xs[i]-xs[i-1])/(xs[i+1]-xs[i-1]);
    p=sig*ys2[i-1]+2.0;
    ys2[i]=(sig-1.0)/p;
    temp[i]=(ys[i+1]-ys[i])/(xs[i+1]-xs[i])-(ys[i]-ys[i-1])/(xs[i]-xs[i-1]);
    temp[i]=(6.0*temp[i]/(xs[i+1]-xs[i-1])-sig*temp[i-1])/p;
  }
  ys2[n-1]=0.0;
  for (j=n-2;j>=0;j--) ys2[j]=ys2[j]*ys2[j+1]+temp[j];

  /* Compute the spline coefficients */
  cur=0;
  j=0;
  jfin=n-1;
  while (xs[j+2]<x[0]) j++;
  while (xs[jfin]>x[nx-1]) jfin--;
  for (;j<=jfin;j++) {
    /* Compute the coefficients of the polynomial between two knots */
    a=xs[j];
    b=xs[j+1];
    c=b-a;
    d=ys[j];
    e=ys[j+1];
    f=ys2[j];
    g=ys2[j+1];
    a0=(b*d-a*e+CUBE(b)*f/6-CUBE(a)*g/6)/c+c*(a*g-b*f)/6;
    a1=(e-d-SQUARE(b)*f/2+SQUARE(a)*g/2)/c+c*(f-g)/6;
    a2=(b*f-a*g)/(2*c);
    a3=(g-f)/(6*c);


    prev=cur;
    while ((cur<nx) && ((j==jfin) || (x[cur]<xs[j+1]))) cur++;

    /* Compute the value of the spline at the sampling times x[i] */
    for (i=prev;i<cur;i++) {
      xc=x[i];
      y[i]=a0+a1*xc+a2*SQUARE(xc)+a3*CUBE(xc);
    }
  }
}
