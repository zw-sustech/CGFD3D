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

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "abs_funcs.h"

//- may move to par file
#define CONSPD 2.0f // power for d
#define CONSPB 2.0f // power for beta
#define CONSPA 1.0f // power for alpha

#ifndef PI
#define PI 3.1415926535898
#endif

/*
 * set up abs_coefs for cfs-pml
 */

float
abs_cal_pml_R(int N)
{
  // use corrected Rpp
  return (float) (pow(10, -( (log10((double)N)-1.0)/log10(2.0) + 4.0)));
}

float
abs_cal_pml_dmax(float L, float Vp, float Rpp)
{
  return (float) (-Vp / (2.0 * L) * log(Rpp) * (CONSPD + 1.0));
}

float
abs_cal_pml_amax(float fc)
{return PI*fc;}

float
abs_cal_pml_d(float x, float L, float dmax)
{
  return (x<0) ? 0.0f : (float) (dmax * pow(x/L, CONSPD));
}

float
abs_cal_pml_a(float x, float L, float amax)
{
  return (x<0) ? 0.0f : (float) (amax * (1.0 - pow(x/L, CONSPA)));
}

float
abs_cal_pml_b(float x, float L, float bmax)
{
  return (x<0) ? 1.0f : (float) (1.0 + (bmax-1.0) * pow(x/L, CONSPB));
}

void
abs_set_cfspml(
    float *alpha_max, //
    float *beta_max, //
    float *velocity, //
    float dh, // average grid size, need to change use x3d etc to cal
    int    ni1,
    int    ni2,
    int    nj1,
    int    nj2,
    int    nk1,
    int    nk2,
    int *boundary_itype, // input
    int *abs_num_of_layers, // output
    int    *abs_indx,
    int    *abs_coefs_facepos0,
    float **p_abs_coefs,
    int verbose)
{
  const int num_of_coefs = 3; // Alpha, Beta, D

  // coef size
  int    size_of_coef = 0;
  for (int i=0; i<CONST_NDIM_2; i++)
  {
    // facepos starts from previous end
    abs_coefs_facepos0[i] = size_of_coef;

    if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML && abs_num_of_layers[i]>0) {
      size_of_coef += abs_num_of_layers[i] * num_of_coefs;
    }
  }

  // alloc coef
  size_of_coef *= num_of_coefs;
  float *abs_coefs = (float *) fdlib_mem_malloc_1d(size_of_coef*sizeof(float),"abs_set_cfspml");

  // for each face
  for (int i_dim=0; i_dim<CONST_NDIM_2; i_dim++)
  {
    // set default indx to 0
    for (int j=0; j<CONST_NDIM_2; j++)
    {
      abs_indx[i_dim * CONST_NDIM_2 + j] = 0;

      // end size
      if (j%2==1) abs_indx[i_dim * CONST_NDIM_2 + j] = -1;
    }

    if (boundary_itype[i_dim] == FD_BOUNDARY_TYPE_CFSPML && abs_num_of_layers[i_dim]>0)
    {
      int num_lay = abs_num_of_layers[i_dim];

      // coefs for this face
      float *this_coefs = abs_coefs + abs_coefs_facepos0[i_dim];

      // indx
      abs_indx[i_dim*CONST_NDIM_2+0] = ni1;
      abs_indx[i_dim*CONST_NDIM_2+1] = ni2;
      abs_indx[i_dim*CONST_NDIM_2+2] = nj1;
      abs_indx[i_dim*CONST_NDIM_2+3] = nj2;
      abs_indx[i_dim*CONST_NDIM_2+4] = nk1;
      abs_indx[i_dim*CONST_NDIM_2+5] = nk2;

      if        (i_dim==0) {
        abs_indx[i_dim*CONST_NDIM_2+1] = ni1 + num_lay - 1;
      } else if (i_dim==1) {
        abs_indx[i_dim*CONST_NDIM_2+0] = ni2 - num_lay + 1;
      } else if (i_dim==2) {
        abs_indx[i_dim*CONST_NDIM_2+3] = nj1 + num_lay - 1;
      } else if (i_dim==3) {
        abs_indx[i_dim*CONST_NDIM_2+2] = nj2 - num_lay + 1;
      } else if (i_dim==4) {
        abs_indx[i_dim*CONST_NDIM_2+5] = nk1 + num_lay - 1;
      } else { // 5
        abs_indx[i_dim*CONST_NDIM_2+4] = nk2 - num_lay + 1;
      }

      float *A = this_coefs;
      float *B = A + num_lay;
      float *D = B + num_lay;

      // calculate
      float L0   = num_lay * dh;
      float Rpp  = abs_cal_pml_R(num_lay);
      float dmax = abs_cal_pml_dmax(L0, velocity[i_dim], Rpp);
      float amax = alpha_max[i_dim];
      float bmax = beta_max[i_dim];

      // from PML-interior to outer side
      for (int ilay=0; ilay<num_lay; ilay++)
      {
        float L = (ilay + 1) * dh;

        // convert to grid index from left to right, default x1/y1/z2
        int i = num_lay - 1 - ilay;
        if (i_dim % 2 == 1) { // x2/y2/z2
          i = ilay; 
        }

        D[i] = abs_cal_pml_d( L, L0, dmax );
        A[i] = abs_cal_pml_a( L, L0, amax );
        B[i] = abs_cal_pml_b( L, L0, bmax );

        // convert d_x to d_x/beta_x since only d_x/beta_x needed
        D[i] /= B[i];
        // covert ax = a_x + d_x/beta_x 
        A[i] += D[i];
        // covert bx = 1.0/bx 
        B[i] = 1.0 / B[i];
      }
    }
  }

  *p_abs_coefs = abs_coefs;
}

// set cfs pml vars
void
abs_init_vars_cfspml(
    int number_of_levels,
    int number_of_vars,
    int *restrict boundary_itype,
    int    *restrict abs_indx,
    int    *restrict abs_vars_volsiz,
    int    *restrict abs_vars_facepos0,
    int    *abs_vars_size_per_level,
    float **restrict p_abs_vars,
    const int myid, const int verbose)
{
    int    siz_level = 0;

    for (int i=0; i < CONST_NDIM_2; i++)
    {
      abs_vars_facepos0[i] = siz_level;

      // init to 0
      abs_vars_volsiz[i] = 0;

      // set if pml
      if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML)
      {
        abs_vars_volsiz[i] =   (abs_indx[i*CONST_NDIM_2+1] - abs_indx[i*CONST_NDIM_2+0] + 1)
                             * (abs_indx[i*CONST_NDIM_2+3] - abs_indx[i*CONST_NDIM_2+2] + 1)
                             * (abs_indx[i*CONST_NDIM_2+5] - abs_indx[i*CONST_NDIM_2+4] + 1);

        // add to total size
        siz_level += abs_vars_volsiz[i] * number_of_vars;
      }
    }

    *abs_vars_size_per_level = siz_level;

    // vars
    // contain all vars at each side, include rk scheme 4 levels vars
    *p_abs_vars = (float *) fdlib_mem_calloc_1d_float( 
                 siz_level * number_of_levels,
                 0.0, "abs_init_vars_cfspml");
}

//
// abl exp type
//
/*
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
  memcpy(abs_numbers, in_abs_numbers, CONST_NDIM_2 * sizeof(int));

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
*/
