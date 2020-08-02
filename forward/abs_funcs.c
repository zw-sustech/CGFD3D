/*
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "netcdf.h"
#include "fdlib_math.h"

//#define DEBUG
#define CONSPD 2.0f // power for d
#define CONSPB 2.0f // power for beta
#define CONSPA 1.0f // power for alpha
#define AbsVzero
//#define PI 3.1415926535898

int abs_set_cfspml(size_t nx, size_t ny, size_t nz, 
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2, 
    int *boundary_itype, // input
    int *in_abs_numbers, //
    float *abs_alpha, //
    float *abs_beta, //
    float *abs_velocity, //
    int *abs_numbers, // output
    size_t *abs_coefs_dimpos, 
    float **p_abs_coefs)
{
  int ivar;
  
  float *Ax, *Bx, *Dx;
  float *Ay, *By, *Dy;
  float *Az, *Bz, *Dz;

  int num_of_coefs = 3; // Alpha, Beta, D
  
  size_t abs_ceofs_size = 0;
  
  // copy input to struct
  memcpy(abs_numbers, in_abs_numbers, 2 * 3 * sizeof(int));

  // size
  abs_ceofs_size = num_of_coefs * (nx + ny + nz); 
  *p_abs_coef_size = abs_coefs_size;

  // set pos, to simplify the alloc of abs_coefs, each dim accouts its layer no matter pml or exp
  //  can't mix exp and pml at current
  abs_coefs_dimpos[0] = 0; // x
  abs_coefs_dimpos[1] = abs_coefs_dimpos[0] + num_of_coefs * nx; // y
  abs_coefs_dimpos[2] = abs_coefs_dimpos[1] + num_of_coefs * ny; // z
  
  // vars
  *p_abs_coefs = (float *) fdlib_mem_calloc_1d_float( 
               abs_coefs_size, 0.0, "abs_set_cfspml");
  if (*p_abs_coefs == NULL) {
      fprintf(stderr,"Error: failed to alloc cfspml coefs\n");
      fflush(stderr);
      ierr = -1;
  }
  
  // cal ceofs
  
  // point to
  Ax = *p_abs_coefs + abs_coefs_dimpos[0];
  Bx = Ax + nx;
  Dx = Bx + nx;
  
  Ay = *p_abs_coefs + abs_coefs_dimpos[1];
  By = Ay + ny;
  Dy = By + ny;
  
  Az = *p_abs_coefs + abs_coefs_dimpos[2];
  Bz = Az + nz;
  Dz = Bz + nz;
  
  for (int i = 0; i < nx; i++){
    Ax[i] = 0.0f;
    Bx[i] = 1.0f;
    Dx[i] = 0.0f;
  }

  for (int j = 0; j < ny; j++){
    Ay[j] = 0.0f;
    By[j] = 1.0f;
    Dy[j] = 0.0f;
  }

  for (int k = 0; k < nz; k++){
    Az[k] = 0.0f;
    Bz[k] = 1.0f;
    Dz[k] = 0.0f;
  }

  // x1
  if (boundary_itype[0] == BOUNDARY_TYPE_CFSPML && abs_numbers[0]>0)
  {
    float Rpp  = cal_pml_R(abs_numbers[0]);
    float L0 = abs_numbers[0];
    float amax = cal_pml_amax(PML_fc);
    float dmax = cal_pml_dmax(L0, abs_velocity[0], Rpp);
    float bmax = abs_beta[0];
    // from PML-interior to outer side
    for (int ilay=0; ilay<abs_numbers[0]; ilay++)
    {
      float Lx = ilay + 1;
      // convert to grid index
      i = ni1 + abs_numbers[0] - 1 - ilay;
      Dx[i] = cal_pml_d( Lx, L0, dmax );
      Ax[i] = cal_pml_a( Lx, L0, amax );
      Bx[i] = cal_pml_b( Lx, L0, bmax );
    }
  }

  // need to add x2, y1, y2 etc


  // convert d_x to d_x/beta_x since only d_x/beta_x needed
  for (int i = 0; i < nx; i++) Dx[i] /= Bx[i];
  for (int i = 0; i < ny; i++) Dy[i] /= By[i];
  for (int i = 0; i < nz; i++) Dz[i] /= Bz[i];

  // covert ax = a_x + d_x/beta_x to reduce computational cost
  for (int i = 0; i < nx; i++) Ax[i] += Dx[i];
  for (int i = 0; i < ny; i++) Ay[i] += Dy[i];
  for (int i = 0; i < nz; i++) Az[i] += Dz[i];

#ifdef DEBUG
  FILE *fp;
  char fnm[1000];
  sprintf(fnm, "%s/ABD%03d%03d%03d.txt", OUT, thisid[0], thisid[1], thisid[2]);
  //sprintf(fnm, "output/ABD.txt");
  fp = fopen(fnm, "w");
  for (int i = 0; i < ni; i++)
    fprintf(fp, "ABDx[%05d] = %f %f %f\n", i, P.Ax[i], P.Bx[i], P.Dx[i]);
  for (int i = 0; i < nj; i++)
    fprintf(fp, "ABDy[%05d] = %f %f %f\n", i, P.Ay[i], P.By[i], P.Dy[i]);
  for (int i = 0; i < nk; i++)
    fprintf(fp, "ABDz[%05d] = %f %f %f\n", i, P.Az[i], P.Bz[i], P.Dz[i]);
  fclose(fp);
#endif

  return 0;
}

inline float cal_pml_R(int N){
  return (float) (pow(10, -( (log10((double)N)-1.0)/log10(2.0) + 3.0)));
}

inline float cal_pml_dmax(float L, float Vp, float Rpp){
  return (float) (-Vp / (2.0 * L) * log(Rpp) * (CONSPD + 1.0));
}

inline float cal_pml_amax(float fc){return PI*fc;}

//inline float cal_pml_bmax(){return 3.0;}

inline float cal_pml_d(float x, float L, float dmax){
  return (x<0) ? 0.0f : (float) (dmax * pow(x/L, CONSPD));
}

inline float cal_pml_a(float x, float L, float amax){
  return (x<0) ? 0.0f : (float) (amax * (1.0 - pow(x/L, CONSPA)));
}

inline float cal_pml_b(float x, float L, float bmax){
  return (x<0) ? 1.0f : (float) (1.0 + (bmax-1.0) * pow(x/L, CONSPB));
}

