/*
 *
 */

// todo:
//  check : abs_set_ablexp
//  convert fortrn to c: abs_ablexp_cal_damp

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "bdry_pml.h"

//- may move to par file
#define CONSPD 2.0f // power for d
#define CONSPB 2.0f // power for beta
#define CONSPA 1.0f // power for alpha

/*
 * set up abs_coefs for cfs-pml
 */

float
bdry_pml_cal_R(int N)
{
  // use corrected Rpp
  return (float) (pow(10, -( (log10((double)N)-1.0)/log10(2.0) + 4.0)));
}

float
bdry_pml_cal_dmax(float L, float Vp, float Rpp)
{
  return (float) (-Vp / (2.0 * L) * log(Rpp) * (CONSPD + 1.0));
}

float
bdry_pml_cal_amax(float fc)
{return PI*fc;}

float
bdry_pml_cal_d(float x, float L, float dmax)
{
  return (x<0) ? 0.0f : (float) (dmax * pow(x/L, CONSPD));
}

float
bdry_pml_cal_a(float x, float L, float amax)
{
  return (x<0) ? 0.0f : (float) (amax * (1.0 - pow(x/L, CONSPA)));
}

float
bdry_pml_cal_b(float x, float L, float bmax)
{
  return (x<0) ? 1.0f : (float) (1.0 + (bmax-1.0) * pow(x/L, CONSPB));
}

void
bdry_pml_set(gdinfo_t *gdinfo,
             gdcurv_t *gdcurv,
             wfel1st_t *wfel1st,
             bdrypml_t *bdrypml,
             int   *neighid, 
             int   in_is_sides[][2],
             int   in_num_layers[][2],
             float in_alpha_max[][2], //
             float in_beta_max[][2], //
             float in_velocity[][2], //
             int verbose)
{
  float *restrict x3d = gdcurv->x3d;
  float *restrict y3d = gdcurv->y3d;
  float *restrict z3d = gdcurv->z3d;
  int    ni1 = gdinfo->ni1;
  int    ni2 = gdinfo->ni2;
  int    nj1 = gdinfo->nj1;
  int    nj2 = gdinfo->nj2;
  int    nk1 = gdinfo->nk1;
  int    nk2 = gdinfo->nk2;
  int    nx  = gdinfo->nx ;
  int    ny  = gdinfo->ny ;
  int    nz  = gdinfo->nz ;
  int    siz_line = gdinfo->siz_iy;
  int    siz_slice = gdinfo->siz_iz;

  // default disable
  bdrypml->is_enable = 0;

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int ind_1d = iside + idim * 2;

      // default set to input
      bdrypml->is_at_sides  [idim][iside] = in_is_sides[idim][iside];
      bdrypml->num_of_layers[idim][iside] = in_num_layers[idim][iside];

      // reset 0 if not mpi boundary
      if (neighid[ind_1d] != MPI_PROC_NULL)
      {
        bdrypml->is_at_sides  [idim][iside] = 0;
        bdrypml->num_of_layers[idim][iside] = 0;
      }

      // default loop index
      bdrypml->ni1[idim][iside] = ni1;
      bdrypml->ni2[idim][iside] = ni2;
      bdrypml->nj1[idim][iside] = nj1;
      bdrypml->nj2[idim][iside] = nj2;
      bdrypml->nk1[idim][iside] = nk1;
      bdrypml->nk2[idim][iside] = nk2;

      // shrink to actual size
      if (idim == 0 && iside ==0) { // x1
        bdrypml->ni2[idim][iside] = ni1 + bdrypml->num_of_layers[idim][iside] - 1;
      }
      if (idim == 0 && iside ==1) { // x2
        bdrypml->ni1[idim][iside] = ni2 - bdrypml->num_of_layers[idim][iside] + 1;
      }
      if (idim == 1 && iside ==0) { // y1
        bdrypml->nj2[idim][iside] = nj1 + bdrypml->num_of_layers[idim][iside] - 1;
      }
      if (idim == 1 && iside ==1) { // y2
        bdrypml->nj1[idim][iside] = nj2 - bdrypml->num_of_layers[idim][iside] + 1;
      }
      if (idim == 2 && iside ==0) { // z1
        bdrypml->nk2[idim][iside] = nk1 + bdrypml->num_of_layers[idim][iside] - 1;
      }
      if (idim == 2 && iside ==1) { // z2
        bdrypml->nk1[idim][iside] = nk2 - bdrypml->num_of_layers[idim][iside] + 1;
      }

      // enable if any side valid
      if (bdrypml->is_at_sides  [idim][iside] == 1) {
        bdrypml->is_enable = 1;
      }

    } // iside
  } // idim

  // alloc coef
  bdrypml->A = (float ***)malloc(CONST_NDIM * sizeof(float**));
  bdrypml->B = (float ***)malloc(CONST_NDIM * sizeof(float**));
  bdrypml->D = (float ***)malloc(CONST_NDIM * sizeof(float**));

  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    bdrypml->A[idim] = (float **)malloc(2 * sizeof(float*));
    bdrypml->B[idim] = (float **)malloc(2 * sizeof(float*));
    bdrypml->D[idim] = (float **)malloc(2 * sizeof(float*));

    for (int iside=0; iside<2; iside++)
    {
      if (bdrypml->is_at_sides[idim][iside] == 1) {
        int nlay = bdrypml->num_of_layers[idim][iside];
        bdrypml->A[idim][iside] = (float *)malloc( nlay * sizeof(float));
        bdrypml->B[idim][iside] = (float *)malloc( nlay * sizeof(float));
        bdrypml->D[idim][iside] = (float *)malloc( nlay * sizeof(float));
      } else {
        bdrypml->A[idim][iside] = NULL;
        bdrypml->B[idim][iside] = NULL;
        bdrypml->D[idim][iside] = NULL;
      }
    }
  }

  // cal coef for each dim and side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      // skip if not pml
      if (bdrypml->is_at_sides[idim][iside] == 0) continue;

      int num_lay = bdrypml->num_of_layers[idim][iside];

      // estimate length along grid lines for coef calculation
      int i = bdrypml->ni1[idim][iside];
      int j = bdrypml->nj1[idim][iside];
      int k = bdrypml->nk1[idim][iside];
      int iptr = i + j * siz_line + k * siz_slice;

      float x0 = x3d[iptr];
      float y0 = y3d[iptr];
      float z0 = z3d[iptr];
      
      float L0 = 0.0;
      for (int ilay=1; ilay<num_lay; ilay++)
      {
        // along which dim
        if (idim==0) i = bdrypml->ni1[idim][iside] + ilay; 
        if (idim==1) j = bdrypml->nj1[idim][iside] + ilay; 
        if (idim==2) k = bdrypml->nk1[idim][iside] + ilay; 

        iptr = i + j * siz_line + k * siz_slice;

        float x1 = x3d[iptr];
        float y1 = y3d[iptr];
        float z1 = z3d[iptr];

        L0 += sqrtf( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );

        x0 = x1;
        y0 = y1;
        z0 = z1;
      }

      // estimate spacing
      float dh = L0;
      if (num_lay > 1) {
        dh = L0 / (num_lay-1.0);
      }

      float *A = bdrypml->A[idim][iside];
      float *B = bdrypml->B[idim][iside];
      float *D = bdrypml->D[idim][iside];

      // calculate
      float Rpp  = bdry_pml_cal_R(num_lay);
      float dmax = bdry_pml_cal_dmax(L0, in_velocity[idim][iside], Rpp);
      float amax = in_alpha_max[idim][iside];
      float bmax = in_beta_max[idim][iside];

      // from PML-interior to outer side
      for (int ilay=0; ilay<num_lay; ilay++)
      {
        // first point has non-zero value
        float L = (ilay + 0) * dh;

        // convert to grid index from left to right
        if (iside == 0) { // x1/y1/z1
          i = num_lay - 1 - ilay;
        } else { // x2/y2/z2
          i = ilay; 
        }

        D[i] = bdry_pml_cal_d( L, L0, dmax );
        A[i] = bdry_pml_cal_a( L, L0, amax );
        B[i] = bdry_pml_cal_b( L, L0, bmax );

        // convert d_x to d_x/beta_x since only d_x/beta_x needed
        D[i] /= B[i];
        // covert ax = a_x + d_x/beta_x 
        A[i] += D[i];
        // covert bx = 1.0/bx 
        B[i] = 1.0 / B[i];
      }

    } // iside
  } // idim

  // alloc auxvar
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int nx = (bdrypml->ni2[idim][iside] - bdrypml->ni1[idim][iside] + 1);
      int ny = (bdrypml->nj2[idim][iside] - bdrypml->nj1[idim][iside] + 1);
      int nz = (bdrypml->nk2[idim][iside] - bdrypml->nk1[idim][iside] + 1);

      bdry_pml_auxvar_init(nx,ny,nz,wfel1st,
                           &(bdrypml->auxvar[idim][iside]),verbose);
    } // iside
  } // idim

}

// alloc auxvar
void
bdry_pml_auxvar_init(int nx, int ny, int nz, 
                     wfel1st_t *wfel1st,
                     bdrypml_auxvar_t *auxvar,
                     const int verbose)
{
  auxvar->nx   = nx;
  auxvar->ny   = ny;
  auxvar->nz   = nz;
  auxvar->ncmp = wfel1st->ncmp;
  auxvar->nlevel = wfel1st->nlevel;

  auxvar->siz_iy   = auxvar->nx;
  auxvar->siz_iz   = auxvar->nx * auxvar->ny;
  auxvar->siz_icmp = auxvar->nx * auxvar->ny * auxvar->nz;
  auxvar->siz_ilevel = auxvar->siz_icmp * auxvar->ncmp;

  auxvar->Vx_pos  = wfel1st->Vx_seq  * auxvar->siz_icmp;
  auxvar->Vy_pos  = wfel1st->Vy_seq  * auxvar->siz_icmp;
  auxvar->Vz_pos  = wfel1st->Vz_seq  * auxvar->siz_icmp;
  auxvar->Txx_pos = wfel1st->Txx_seq * auxvar->siz_icmp;
  auxvar->Tyy_pos = wfel1st->Tyy_seq * auxvar->siz_icmp;
  auxvar->Tzz_pos = wfel1st->Tzz_seq * auxvar->siz_icmp;
  auxvar->Txz_pos = wfel1st->Txz_seq * auxvar->siz_icmp;
  auxvar->Tyz_pos = wfel1st->Tyz_seq * auxvar->siz_icmp;
  auxvar->Txy_pos = wfel1st->Txy_seq * auxvar->siz_icmp;

  // vars
  // contain all vars at each side, include rk scheme 4 levels vars
  if (auxvar->siz_icmp > 0 ) { // valid pml layer
    auxvar->var = (float *) fdlib_mem_calloc_1d_float( 
                 auxvar->siz_ilevel * auxvar->nlevel,
                 0.0, "bdry_pml_auxvar_init");
  } else { // nx,ny,nz has 0
    auxvar->var = NULL;
  }
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
