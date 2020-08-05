/*
 *
 */

// todo:
//  check and complete: abs_set_cfspml
//  check and complete: abs_set_ablexp
//  convert fortrn to c: abs_ablexp_cal_damp

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "netcdf.h"
#include "fdlib_math.h"

//#define DEBUG
#define CONSPD 2.0f // power for d
#define CONSPB 2.0f // power for beta
#define CONSPA 1.0f // power for alpha
//#define PI 3.1415926535898

int abs_set_cfspml(size_t nx, size_t ny, size_t nz, 
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2, 
    int *boundary_itype, // input
    int *in_abs_numbers, //
    float *abs_alpha_max, //
    float *abs_beta_max, //
    float *abs_velocity, //
    int *abs_num_of_layers, // output
    size_t **abs_blk_indx,
    //size_t **abs_blk_coefs_size,
    size_t **abs_blk_vars_size,
    float ***p_abs_blk_coefs,
    float ***p_abs_blk_vars)
{
  int ivar;
  
  float *Ax, *Bx, *Dx;
  float *Ay, *By, *Dy;
  float *Az, *Bz, *Dz;
  
  size_t i,j,k;
  size_t size_of_coef;

  const int num_of_coefs = 3; // Alpha, Beta, D
  
  // copy input to struct
  memcpy(abs_num_of_layers, in_abs_numbers, FD_NDIM_2 * sizeof(int));

  // alloc first level coef size
  **p_abs_blk_coefs = (float **) fdlib_mem_malloc_1d_float( 
               FD_NDIM_2, "abs_set_cfspml");

  // x1
  i_dim = 0;
  if (boundary_itype[i_dim] == FD_BOUNDARY_TYPE_CFSPML && abs_num_of_layers[i_dim]>0)
  {
    // indx
    abs_blk_indx[i_dim][0] = ni1;
    abs_blk_indx[i_dim][1] = ni1 + abs_num_of_layers[i_dim] - 1;
    abs_blk_indx[i_dim][2] = nj1;
    abs_blk_indx[i_dim][3] = nj2;
    abs_blk_indx[i_dim][4] = nk1;
    abs_blk_indx[i_dim][5] = nk2;

    // coef_size
    //abs_blk_coefs_size[i_dim] = (abs_num_of_layers[i_dim] + nj + nk) * num_of_coefs; // for MPML
    
    //abs_blk_coefs_size[i_dim] = (abs_num_of_layers[i_dim]) * num_of_coefs;
    //(**p_abs_blk_coefs)[i_dim] = (float *) fdlib_mem_malloc_1d_float( 
    //           abs_blk_coefs_size[i_dim], "abs_set_cfspml");

    size_of_coef = (abs_num_of_layers[i_dim]) * num_of_coefs;
    (**p_abs_blk_coefs)[i_dim] = (float *) fdlib_mem_malloc_1d_float( 
               size_of_coef, "abs_set_cfspml");

    Ax = (**p_abs_blk_coefs)[i_dim];
    Bx = Ax + abs_num_of_layers[i_dim];
    Dx = Bx + abs_num_of_layers[i_dim];

    // init
    for (i = 0; i < abs_num_of_layers[i_dim]; i++) {
      Ax[i] = 0.0f;
      Bx[i] = 1.0f;
      Dx[i] = 0.0f;
    }

    // cal
    float L0   = abs_num_of_layers[i_dim];
    float Rpp  = cal_pml_R(abs_num_of_layers[i_dim]);
    float dmax = cal_pml_dmax(L0, abs_velocity[i_dim], Rpp);
    //float amax = cal_pml_amax(PML_fc);
    float amax = abs_alpha_max[i_dim];
    float bmax = abs_beta_max[i_dim];
    // from PML-interior to outer side
    for (int ilay=0; ilay<abs_num_of_layers[i_dim]; ilay++)
    {
      float Lx = ilay + 1;
      // convert to grid index from left to right
      i = abs_num_of_layers[i_dim] - 1 - ilay;
      Dx[i] = cal_pml_d( Lx, L0, dmax );
      Ax[i] = cal_pml_a( Lx, L0, amax );
      Bx[i] = cal_pml_b( Lx, L0, bmax );

      // convert d_x to d_x/beta_x since only d_x/beta_x needed
      Dx[i] /= Bx[i];

      // covert ax = a_x + d_x/beta_x to reduce computational cost
      Ax[i] += Dx[i];
    }
  }

  // need to add x2, y1, y2, z1, z2


  return 0;
}

inline float cal_pml_R(int N){
  // use corrected Rpp
  return (float) (pow(10, -( (log10((double)N)-1.0)/log10(2.0) + 4.0)));
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

//
// abl exp type
//
int abs_set_ablexp(size_t nx, size_t ny, size_t nz, 
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

  int num_of_coefs = 1; // damping
  
  size_t abs_ceofs_size = 0;
  
  // copy input to struct
  memcpy(abs_numbers, in_abs_numbers, FD_NDIM_2 * sizeof(int));

  // size
  for (i=0; i<2; i++) { // x1,x2
    abs_coefs_dimpos[i] = abs_ceofs_size;
    abs_ceofs_size += abs_numbers[i] * nj * nk; 
  }
  for (i=2; i<4; i++) { // y1,y2
    abs_coefs_dimpos[i] = abs_ceofs_size;
    abs_ceofs_size += abs_numbers[i] * ni * nk; 
  }
  for (i=4; i<6; i++) { // z1,z2
    abs_coefs_dimpos[i] = abs_ceofs_size;
    abs_ceofs_size += abs_numbers[i] * ni * nj; 
  }

  *p_abs_coef_size = abs_coefs_size;

  // vars
  *p_abs_coefs = (float *) fdlib_mem_calloc_1d_float( 
               abs_coefs_size, 0.0, "abs_set_ablexp");
  if (*p_abs_coefs == NULL) {
      fprintf(stderr,"Error: failed to alloc ablexp coefs\n");
      fflush(stderr);
      ierr = -1;
  }

  // set damping values

  return ierr;
}

int abs_ablexp_cal_damp(i,Vs,ah,nb)
{
  int ierr = 0;

  integer,intent(in) :: i,nb
  real(SP),intent(in) :: Vs,ah
  real(SP) :: d
  real(SP) :: ie
  integer m,n
  ie=i
  !Vs=5000.0_SP
  m=(nb*ah)/(Vs*stept)
  d=0.0_SP
  do n=1,m
     d=d+(n*stept*Vs)**2/(nb*ah)**2
  end do
  d=0.8_SP/d*1.1_SP
  d=exp(-d*(ie/nb)**2)

  return ierr;
}
