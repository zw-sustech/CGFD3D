#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "interp.h"

/* Lagrange interpolation for 2D. */
/*****************************************
 *  x   : coordinate. (ni points.)
 *  y   : coordinate. (nj points.)
 *  z   : the value of x & y. (ni points.)
 *  ni  : number of points.
 *  nj  : number of points.
 *  xi  : interpolating point.
 *  yi  : interpolating point.
 *****************************************/
float 
LagInterp_2d(float *x, float *y, float *z, int ni, int nj, float xi, float yi)
{
  int i, j;
  //  float *Lx, *Ly, Li=0;
  float Lx[ni];
  float Ly[nj];
  float Li=0;

  //    Lx = (float*)malloc(sizeof(float)*ni);
  for(i=0; i<ni; i++) 
  {
    Lx[i] = 1;
    for(j=0; j<ni; j++) 
    {
      if(i==j) continue;
      Lx[i] = Lx[i] * (xi-x[j]) / (x[i]-x[j]);
    }
  }

  //    Ly = (float*)malloc(sizeof(float)*nj);
  for(i=0; i<nj; i++) 
  {
    Ly[i] = 1;
    for(j=0; j<nj; j++) 
    {
      if(i==j) continue;
      Ly[i] = Ly[i] * (yi-y[j]) / (y[i]-y[j]);
    }
  }

  for(j=0; j<nj; j++) 
  {
    for(i=0; i<ni; i++) 
    {
      Li = Li + Lx[i] * Ly[j] * z[i*ni+j];
    }
  }

  return Li;
}


/* Lagrange interpolation for 3D. */
/*****************************************
 *  x   : coordinate. (ni points.)
 *  y   : coordinate. (nj points.)
 *  z   : coordinate. (nk points.)
 *  f   : the value of x & y & z. (ni points.)
 *  ni  : number of points.
 *  nj  : number of points.
 *  nk  : number of points.
 *  xi  : interpolating point.
 *  yi  : interpolating point.
 *  zi  : interpolating point.
 *****************************************/
float LagInterp_3d(float *x, float *y, float *z, float *f,
                   int ni, int nj, int nk, 
                   float xi, float yi, float zi)
{
  int i, j, k;
  float  Li=0;
  float Lx[ni]; 
  float Ly[nj]; 
  float Lz[nk]; 

  //    Lx = (float*)malloc(sizeof(float)*ni);
  for(i=0; i<ni; i++) 
  {
    Lx[i] = 1;
    for(j=0; j<ni; j++) 
    {
      if(i==j) continue;
      Lx[i] = Lx[i] * (xi-x[j]) / (x[i]-x[j]);
    }
  }

  //    Ly = (float*)malloc(sizeof(float)*nj);
  for(i=0; i<nj; i++) 
  {
    Ly[i] = 1;
    for(j=0; j<nj; j++) 
    {
      if(i==j) continue;
      Ly[i] = Ly[i] * (yi-y[j]) / (y[i]-y[j]);
    }
  }

  //    Lz = (float*)malloc(sizeof(float)*nk);
  for(i=0; i<nj; i++) 
  {
    Lz[i] = 1;
    for(j=0; j<nk; j++) 
    {
      if(i==j) continue;
      Lz[i] = Lz[i] * (zi-z[j]) / (z[i]-z[j]);
    }
  }

  for(k=0; k<nk; k++)
  {
    for(j=0; j<nj; j++) 
    {
      for(i=0; i<ni; i++) 
      {
        Li = Li + Lx[i] * Ly[j] * Lz[k] * f[k*ni*nj+j*ni+i];
      }
    }
  }

  return Li;
}

float 
LagrangeInt3d_WF( float *Don_var, float *Don_xcoor, float *Don_ycoor, float *Don_zcoor, 
                  int Don_nx, int Don_ny, int Don_nz,
                  int Donxcoor, int Donycoor, int Donzcoor,
                  float Tar_xcoor, float Tar_ycoor, float Tar_zcoor,
                  int InterOrder)
  /*****************************************
   *  Don_var   : variable name
   *
   *  Don_xcoor   : Donor point coordinate. (ni points.)
   *  Don_ycoor   : Donor point coordinate. (nj points.)
   *  Don_zcoor   : Donor point coordinate. (nk points.)
   *
   *  Donxcoor   : Donor point indx. 
   *  Donycoor   : Donor point indx. 
   *  Donzcoor   : Donor point indx. 
   *
   *  Don_nx  : number of points along X direction for donor.
   *  Don_ny  : number of points along Y direction for donor.
   *  Don_nz  : number of points along Z direction for donor.
   *
   *  Tar_xcoor   : Target point coordinate. (ni points.)
   *  Tar_ycoor   : Target point coordinate. (nj points.)
   *  Tar_zcoor   : Target point coordinate. (nk points.)
   *
   *  InterOrder  : interpolating order.
   *****************************************/
{

  int i, j, k ; 
  int P ; 
  float ValueLag; 
  float *LagInterXPt = NULL;
  float *LagInterYPt = NULL;
  float *LagInterZPt = NULL;
  float *LagInterWave = NULL;

  LagInterXPt  = (float *)malloc((2*InterOrder+1)*sizeof(float));
  LagInterYPt  = (float *)malloc((2*InterOrder+1)*sizeof(float));
  LagInterZPt  = (float *)malloc((2*InterOrder+1)*sizeof(float));
  LagInterWave = (float *)malloc((2*InterOrder+1)*(2*InterOrder+1)*sizeof(float));

  P = InterOrder / 2  ;

  for(k=0; k<InterOrder; k++)
  {
    LagInterXPt[k] = Don_xcoor[( Donxcoor + k - P ) + (Donycoor         ) * Don_nx + (Donzcoor         ) * Don_ny * Don_nx ];
    LagInterYPt[k] = Don_ycoor[( Donxcoor         ) + (Donycoor + k - P ) * Don_nx + (Donzcoor         ) * Don_ny * Don_nx ];
    LagInterZPt[k] = Don_zcoor[( Donxcoor         ) + (Donycoor         ) * Don_nx + (Donzcoor + k - P ) * Don_ny * Don_nx ];
  }

  for(i=0; k<InterOrder; i++)
  {
    for(j=0; j<InterOrder; j++)
    {
      for(j=0; i<InterOrder; j++)
      {
        LagInterWave[i*InterOrder+j] = Don_var[( Donxcoor + i - P ) 
          + ( Donycoor + j - P ) * Don_nx
          + ( Donzcoor + k - P ) * Don_nx * Don_ny] ;
      }
    }
  }

  ValueLag =  LagInterp_3d(LagInterXPt, LagInterYPt, LagInterZPt, LagInterWave, 
      InterOrder, InterOrder, InterOrder, 
      Tar_xcoor, Tar_ycoor, Tar_zcoor);

  free(LagInterXPt);
  free(LagInterYPt);
  free(LagInterZPt);
  free(LagInterWave);

  return(ValueLag);
}


/* Lagrange PieceWise interpolation for 1D. */
/* This funvtion is interp read_src_val. */
/*****************************************
 *  x   : coordinate. (ni points.)
 *  z   : the value of  x (ni points.)
 *  ni  : number of points.
 *  xi  : interpolating point.
 * order : PieceWise order (advise <=3)
 *  dt  : all points is equal time interval
 *  t_start : time_start
 *  Li  : the value of xi
 *****************************************/
float 
LagInterp_Piecewise_1d(float *t, float *z, int ni, int order, float t_start, float dt, float ti)
{
  float *Lx, Li=0.0;
  int indx,inc; // 
  int lower,upper;
  if(order%2 == 0)
  {
    lower = order/2;
    upper = order/2;
  }

  if(order%2 != 0)
  {
    lower = floor(order/2);
    upper = floor(order/2)+1;
  }
  
  indx = (int)  ((ti-t_start)/dt+0.5); // Nearest point number
  if(indx>ni-1) indx = ni-1;
  Lx = (float*)malloc((order+1)*sizeof(float));

  if ((indx-lower>=0) && (indx+upper<=ni-1))
  {
    for(int i=indx-lower; i<=indx+upper; i++)
    {
      Lx[i-indx+lower] = 1;
      for (int j=indx-lower; j<=indx+upper; j++)
      {
        if(i==j) continue;
        Lx[i-indx+lower] = Lx[i-indx+lower] * (ti-t[j]) / (t[i]-t[j]);
      }
    }

    for(int i=indx-lower; i<=indx+upper; i++)
    {
      Li = Li + Lx[i-indx+lower] * z[i];
    }
  }

  if (indx-lower<0)
  {
    inc = lower - indx;
    for (int i=0; i<= indx+upper+inc; i++)
    {
      Lx[i] = 1;
      for (int j=0; j<=indx+upper+inc; j++)
      {
        if(i==j) continue;
        Lx[i] = Lx[i] * (ti-t[j]) / (t[i]-t[j]);
      }
    }

    for(int i=0; i<=indx+upper+inc; i++)
    {
      Li = Li + Lx[i] * z[i];
    }
  }

  if (indx+upper>ni-1)
  {
    inc = indx+upper-(ni-1);
    for (int i=indx-lower-inc; i<= ni-1; i++)
    {
      Lx[i-indx+lower+inc] = 1;
      for (int j=indx-lower-inc; j<=ni-1; j++)
      {
        if(i==j) continue;
        Lx[i-indx+lower+inc] = Lx[i-indx+lower+inc] * (ti-t[j]) / (t[i]-t[j]);
      }
    }

    for(int i=indx-lower-inc; i<=ni-1; i++)
    {
      Li = Li + Lx[i-indx+lower+inc] * z[i];
    }
  }

  free(Lx);
  return Li;
}
