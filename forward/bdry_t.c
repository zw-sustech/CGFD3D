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
#include "bdry_t.h"

//- may move to par file
#define CONSPD 2.0f // power for d
#define CONSPB 2.0f // power for beta
#define CONSPA 1.0f // power for alpha

//#define BDRY_T_DEBUG

/*
 * init bdry_t
 */

int
bdry_init(bdry_t *bdry, int nx, int ny, int nz)
{
  bdry->is_enable_pml  = 0;
  bdry->is_enable_mpml = 0;
  bdry->is_enable_ablexp  = 0;
  bdry->is_enable_free = 0;

  bdry->nx = nx;
  bdry->ny = ny;
  bdry->nz = nz;

  return 0;
}

/*
 * matrix for velocity gradient conversion
 *  only implement z2 (top) right now
 */

int
bdry_free_set(gd_t    *gd,
              bdry_t      *bdryfree,
              int   *neighid, 
              int   in_is_sides[][2],
              const int verbose)
{
  int ierr = 0;

  size_t siz_slice  = gd->siz_iz;

  // default disable
  bdryfree->is_enable_free = 0;

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int ind_1d = iside + idim * 2;

      bdryfree->is_sides_free  [idim][iside] = in_is_sides[idim][iside];

      // reset 0 if not mpi boundary
      if (neighid[ind_1d] != MPI_PROC_NULL)
      {
        bdryfree->is_sides_free  [idim][iside] = 0;
      }

      // enable if any side valid
      if (bdryfree->is_sides_free  [idim][iside] == 1) {
        bdryfree->is_enable_free = 1;
      }
    } // iside
  } // idim

  // following only implement z2 (top) right now
  float *matVx2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * CONST_NDIM * CONST_NDIM,
                                      0.0,
                                      "bdry_free_set");

  float *matVy2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * CONST_NDIM * CONST_NDIM,
                                      0.0,
                                      "bdry_free_set");

  bdryfree->matVx2Vz2 = matVx2Vz;
  bdryfree->matVy2Vz2 = matVy2Vz;

  return ierr;
}


void
bdry_pml_set(gd_t *gd,
             wav_t *wav,
             bdry_t *bdrypml,
             int   *neighid, 
             int   in_is_sides[][2],
             int   in_num_layers[][2],
             float in_alpha_max[][2], //
             float in_beta_max[][2], //
             float in_velocity[][2], //
             int verbose)
{
  int    ni1 = gd->ni1;
  int    ni2 = gd->ni2;
  int    nj1 = gd->nj1;
  int    nj2 = gd->nj2;
  int    nk1 = gd->nk1;
  int    nk2 = gd->nk2;
  int    nx  = gd->nx ;
  int    ny  = gd->ny ;
  int    nz  = gd->nz ;
  int    siz_line = gd->siz_iy;
  int    siz_slice = gd->siz_iz;

  // default disable
  bdrypml->is_enable_pml = 0;

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int ind_1d = iside + idim * 2;

      // default set to input
      bdrypml->is_sides_pml  [idim][iside] = in_is_sides[idim][iside];
      bdrypml->num_of_layers[idim][iside] = in_num_layers[idim][iside];

      // reset 0 if not mpi boundary
      if (neighid[ind_1d] != MPI_PROC_NULL)
      {
        bdrypml->is_sides_pml  [idim][iside] = 0;
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
        bdrypml->ni2[idim][iside] = ni1 + bdrypml->num_of_layers[idim][iside];
      }
      if (idim == 0 && iside ==1) { // x2
        bdrypml->ni1[idim][iside] = ni2 - bdrypml->num_of_layers[idim][iside];
      }
      if (idim == 1 && iside ==0) { // y1
        bdrypml->nj2[idim][iside] = nj1 + bdrypml->num_of_layers[idim][iside];
      }
      if (idim == 1 && iside ==1) { // y2
        bdrypml->nj1[idim][iside] = nj2 - bdrypml->num_of_layers[idim][iside];
      }
      if (idim == 2 && iside ==0) { // z1
        bdrypml->nk2[idim][iside] = nk1 + bdrypml->num_of_layers[idim][iside];
      }
      if (idim == 2 && iside ==1) { // z2
        bdrypml->nk1[idim][iside] = nk2 - bdrypml->num_of_layers[idim][iside];
      }

      // enable if any side valid
      if (bdrypml->is_sides_pml  [idim][iside] == 1) {
        bdrypml->is_enable_pml = 1;
      }

    } // iside
  } // idim

  // alloc coef
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      if (bdrypml->is_sides_pml[idim][iside] == 1) {
        int npoints = bdrypml->num_of_layers[idim][iside] + 1;
        bdrypml->A[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->B[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->D[idim][iside] = (float *)malloc( npoints * sizeof(float));
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
      if (bdrypml->is_sides_pml[idim][iside] == 0) continue;

      float *A = bdrypml->A[idim][iside];
      float *B = bdrypml->B[idim][iside];
      float *D = bdrypml->D[idim][iside];

      // esti L0 and dh
      float L0, dh;
      bdry_cal_abl_len_dh(gd,bdrypml->ni1[idim][iside],
                             bdrypml->ni2[idim][iside],
                             bdrypml->nj1[idim][iside],
                             bdrypml->nj2[idim][iside],
                             bdrypml->nk1[idim][iside],
                             bdrypml->nk2[idim][iside],
                             idim,
                             &L0, &dh);

      // para
      int npoints = bdrypml->num_of_layers[idim][iside] + 1;
      float num_lay = npoints - 1;
      float Rpp  = bdry_pml_cal_R(num_lay);
      float dmax = bdry_pml_cal_dmax(L0, in_velocity[idim][iside], Rpp);
      float amax = in_alpha_max[idim][iside];
      float bmax = in_beta_max[idim][iside];

      // from PML-interior to outer side
      for (int n=0; n<npoints; n++)
      {
        // first point has non-zero value
        float L = n * dh;
        int i;

        // convert to grid index from left to right
        if (iside == 0) { // x1/y1/z1
          i = npoints - 1 - n;
        } else { // x2/y2/z2
          i = n; 
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

      bdry_pml_auxvar_init(nx,ny,nz,wav,
                           &(bdrypml->auxvar[idim][iside]),verbose);
    } // iside
  } // idim

}

// alloc auxvar
void
bdry_pml_auxvar_init(int nx, int ny, int nz, 
                     wav_t *wav,
                     bdrypml_auxvar_t *auxvar,
                     const int verbose)
{
  auxvar->nx   = nx;
  auxvar->ny   = ny;
  auxvar->nz   = nz;
  auxvar->ncmp = wav->ncmp;
  auxvar->nlevel = wav->nlevel;

  auxvar->siz_iy   = auxvar->nx;
  auxvar->siz_iz   = auxvar->nx * auxvar->ny;
  auxvar->siz_icmp = auxvar->nx * auxvar->ny * auxvar->nz;
  auxvar->siz_ilevel = auxvar->siz_icmp * auxvar->ncmp;

  auxvar->Vx_pos  = wav->Vx_seq  * auxvar->siz_icmp;
  auxvar->Vy_pos  = wav->Vy_seq  * auxvar->siz_icmp;
  auxvar->Vz_pos  = wav->Vz_seq  * auxvar->siz_icmp;
  auxvar->Txx_pos = wav->Txx_seq * auxvar->siz_icmp;
  auxvar->Tyy_pos = wav->Tyy_seq * auxvar->siz_icmp;
  auxvar->Tzz_pos = wav->Tzz_seq * auxvar->siz_icmp;
  auxvar->Tyz_pos = wav->Tyz_seq * auxvar->siz_icmp;
  auxvar->Txz_pos = wav->Txz_seq * auxvar->siz_icmp;
  auxvar->Txy_pos = wav->Txy_seq * auxvar->siz_icmp;

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

/*
 * set up abs_coefs for cfs-pml
 */

float
bdry_pml_cal_R(float num_lay)
{
  // use corrected Rpp
  return (float) (pow(10, -( (log10((double)num_lay)-1.0)/log10(2.0) + 4.0)));
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
/*
 * esti L and dh along idim damping layers
 */

int
bdry_cal_abl_len_dh(gd_t *gd, 
                    int abs_ni1, int abs_ni2,
                    int abs_nj1, int abs_nj2,
                    int abs_nk1, int abs_nk2,
                    int idim,
                    float *avg_L, float *avg_dh)
{
  int ierr = 0;

  int siz_line  = gd->siz_iy;
  int siz_slice = gd->siz_iz;

  // cartesian grid is simple
  if (gd->type == GD_TYPE_CART)
  {
    if (idim == 0) { // x-axis
      *avg_dh = gd->dx;
      *avg_L  = gd->dx * (abs_ni2 - abs_ni1);
    } else if (idim == 1) { // y-axis
      *avg_dh = gd->dy;
      *avg_L  = gd->dy * (abs_nj2 - abs_nj1);
    } else { // z-axis
      *avg_dh = gd->dz;
      *avg_L  = gd->dz * (abs_nk2 - abs_nk1);
    }
  }
  // curv grid needs avg
  else if (gd->type == GD_TYPE_CURV)
  {
    float *x3d = gd->x3d;
    float *y3d = gd->y3d;
    float *z3d = gd->z3d;

    double L  = 0.0;
    double dh = 0.0;
    int    num = 0;

    if (idim == 0) // x-axis
    {
      for (int k=abs_nk1; k<=abs_nk2; k++)
      {
        for (int j=abs_nj1; j<=abs_nj2; j++)
        {
          int iptr = abs_ni1 + j * siz_line + k * siz_slice;
          double x0 = x3d[iptr];
          double y0 = y3d[iptr];
          double z0 = z3d[iptr];
          for (int i=abs_ni1+1; i<=abs_ni2; i++)
          {
            int iptr = i + j * siz_line + k * siz_slice;

            double x1 = x3d[iptr];
            double y1 = y3d[iptr];
            double z1 = z3d[iptr];

            L += sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );

            x0 = x1;
            y0 = y1;
            z0 = z1;
            num += 1;
          }
        }
      }

      *avg_dh = (float)( L / num );
      *avg_L = (*avg_dh) * (abs_ni2 - abs_ni1);
    } 
    else if (idim == 1) // y-axis
    { 
      for (int k=abs_nk1; k<=abs_nk2; k++)
      {
        for (int i=abs_ni1; i<=abs_ni2; i++)
        {
          int iptr = i + abs_nj1 * siz_line + k * siz_slice;
          double x0 = x3d[iptr];
          double y0 = y3d[iptr];
          double z0 = z3d[iptr];
          for (int j=abs_nj1+1; j<=abs_nj2; j++)
          {
            int iptr = i + j * siz_line + k * siz_slice;

            double x1 = x3d[iptr];
            double y1 = y3d[iptr];
            double z1 = z3d[iptr];

            L += sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );

            x0 = x1;
            y0 = y1;
            z0 = z1;
            num += 1;
          }
        }
      }

      *avg_dh = (float)( L / num );
      *avg_L = (*avg_dh) * (abs_nj2 - abs_nj1);
    }
    else // z-axis
    { 
      for (int j=abs_nj1; j<=abs_nj2; j++)
      {
        for (int i=abs_ni1; i<=abs_ni2; i++)
        {
          int iptr = i + j * siz_line + abs_nk1 * siz_slice;
          double x0 = x3d[iptr];
          double y0 = y3d[iptr];
          double z0 = z3d[iptr];
          for (int k=abs_nk1+1; k<=abs_nk2; k++)
          {
            int iptr = i + j * siz_line + k * siz_slice;

            double x1 = x3d[iptr];
            double y1 = y3d[iptr];
            double z1 = z3d[iptr];

            L += sqrt( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) + (z1-z0)*(z1-z0) );

            x0 = x1;
            y0 = y1;
            z0 = z1;
            num += 1;
          }
        }
      }

      *avg_dh = (float)( L / num );
      *avg_L = (*avg_dh) * (abs_nk2 - abs_nk1);
    } // idim

  } // gd type

  return ierr;
}


/*
 * setup ablexp parameters
 */

int
bdry_ablexp_set(gd_t *gd,
             wav_t *wav,
             bdry_t *bdry,
             int   *neighid, 
             int   in_is_sides[][2],
             int   in_num_layers[][2],
             float in_velocity[][2], //
             float dt,
             int  *topoid,
             int verbose)
{
  int    ierr = 0;

  int    ni1 = gd->ni1;
  int    ni2 = gd->ni2;
  int    nj1 = gd->nj1;
  int    nj2 = gd->nj2;
  int    nk1 = gd->nk1;
  int    nk2 = gd->nk2;
  int    ni  = gd->ni ;
  int    nj  = gd->nj ;
  int    nk  = gd->nk ;
  int    nx  = gd->nx ;
  int    ny  = gd->ny ;
  int    nz  = gd->nz ;
  int    siz_line = gd->siz_iy;
  int    siz_slice = gd->siz_iz;
  int    abs_number[CONST_NDIM][2];

  int n;

  // default disable
  bdry->is_enable_ablexp = 0;

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int ind_1d = iside + idim * 2;

      // default set to input
      bdry->is_sides_ablexp  [idim][iside] = in_is_sides[idim][iside];
      bdry->num_of_layers[idim][iside] = in_num_layers[idim][iside];

      // reset 0 if not mpi boundary
      if (neighid[ind_1d] != MPI_PROC_NULL)
      {
        bdry->is_sides_ablexp  [idim][iside] = 0;
        bdry->num_of_layers[idim][iside] = 0;
      }

      // enable if any side valid
      if (bdry->is_sides_ablexp  [idim][iside] == 1) {
        bdry->is_enable_ablexp = 1;
      }
    } // iside
  } // idim

  // block index for ablexp, default inactive
  bdry_block_t *D = bdry->bdry_blk;
  for (n=0; n < CONST_NDIM_2; n++)
  {
     D[n].enable = 0;
     D[n].ni1 =  0;
     D[n].ni2 = -1;
     D[n].ni  =  0;

     D[n].nj1 =  0;
     D[n].nj2 = -1;
     D[n].nj  =  0;

     D[n].nk1 =  0;
     D[n].nk2 = -1;
     D[n].nk  =  0;
  }

  // alloc coef
  bdry->ablexp_Ex = (float *)malloc( nx * sizeof(float));
  bdry->ablexp_Ey = (float *)malloc( ny * sizeof(float));
  bdry->ablexp_Ez = (float *)malloc( nz * sizeof(float));
  for (int i=0; i<nx; i++) bdry->ablexp_Ex[i] = 1.0;
  for (int j=0; j<ny; j++) bdry->ablexp_Ey[j] = 1.0;
  for (int k=0; k<nz; k++) bdry->ablexp_Ez[k] = 1.0;

  // x1
  n=0;
  D[n].ni=bdry->num_of_layers[0][0]; D[n].ni1=ni1; D[n].ni2=D[n].ni1+D[n].ni-1; 
  D[n].nj=nj                       ; D[n].nj1=nj1; D[n].nj2=D[n].nj1+D[n].nj-1; 
  D[n].nk=nk                       ; D[n].nk1=nk1; D[n].nk2=D[n].nk1+D[n].nk-1; 
  if (D[n].ni>0 && D[n].nj>0 && D[n].nk>0)
  {
     D[n].enable = 1;

     // esti L0 and dh
     float L0, dh;
     bdry_cal_abl_len_dh(gd,D[n].ni1,
                            D[n].ni2,
                            D[n].nj1,
                            D[n].nj2,
                            D[n].nk1,
                            D[n].nk2,
                            0,
                            &L0, &dh);

     for (int i=D[n].ni1; i<=D[n].ni2; i++)
     {
        // the first point of layer is the first dh damping
        bdry->ablexp_Ex[i] = bdry_ablexp_cal_mask(D[n].ni - (i - D[n].ni1),
                                                  in_velocity[0][0], dt,
                                                  D[n].ni, dh);
     }
  }

  // x2
  n += 1;
  D[n].ni=bdry->num_of_layers[0][1];D[n].ni1=ni2 - D[n].ni + 1; D[n].ni2=D[n].ni1+D[n].ni-1; 
  D[n].nj=nj                       ;D[n].nj1=nj1               ; D[n].nj2=D[n].nj1+D[n].nj-1; 
  D[n].nk=nk                       ;D[n].nk1=nk1               ; D[n].nk2=D[n].nk1+D[n].nk-1; 
  if (D[n].ni>0 && D[n].nj>0 && D[n].nk>0)
  {
     D[n].enable = 1;

     // esti L0 and dh
     float L0, dh;
     bdry_cal_abl_len_dh(gd,D[n].ni1,
                            D[n].ni2,
                            D[n].nj1,
                            D[n].nj2,
                            D[n].nk1,
                            D[n].nk2,
                            0,
                            &L0, &dh);

     for (int i=D[n].ni1; i<=D[n].ni2; i++)
     {
        bdry->ablexp_Ex[i] = bdry_ablexp_cal_mask(i - D[n].ni1 + 1,
                                                  in_velocity[0][1], dt,
                                                  D[n].ni, dh);
     }
  }

  int ni_x = bdry->num_of_layers[0][0] + bdry->num_of_layers[0][1];
  int nj_y = bdry->num_of_layers[1][0] + bdry->num_of_layers[1][1];

  //y1
  n += 1;
  D[n].ni = ni - ni_x;
  D[n].ni1= ni1 + bdry->num_of_layers[0][0];
  D[n].ni2= D[n].ni1 + D[n].ni - 1; 
  D[n].nj = bdry->num_of_layers[1][0];
  D[n].nj1= nj1;
  D[n].nj2= D[n].nj1 + D[n].nj - 1; 
  D[n].nk = nk;
  D[n].nk1= nk1;
  D[n].nk2= D[n].nk1 + D[n].nk - 1; 
  if (D[n].ni>0 && D[n].nj>0 && D[n].nk>0)
  {
     D[n].enable = 1;
     // esti L0 and dh
     float L0, dh;
     bdry_cal_abl_len_dh(gd,D[n].ni1,
                            D[n].ni2,
                            D[n].nj1,
                            D[n].nj2,
                            D[n].nk1,
                            D[n].nk2,
                            1,
                            &L0, &dh);

     for (int j=D[n].nj1; j<=D[n].nj2; j++)
     {
        bdry->ablexp_Ey[j] = bdry_ablexp_cal_mask(D[n].nj - (j - D[n].nj1),
                                                  in_velocity[1][0], dt,
                                                  D[n].nj, dh);
     }
  }

  // y2
  n += 1;
  D[n].ni = ni - ni_x;
  D[n].ni1= ni1 + bdry->num_of_layers[0][0];
  D[n].ni2= D[n].ni1 + D[n].ni - 1; 
  D[n].nj = bdry->num_of_layers[1][1];
  D[n].nj1= nj2 - D[n].nj + 1;
  D[n].nj2= D[n].nj1 + D[n].nj - 1; 
  D[n].nk = nk;
  D[n].nk1= nk1;
  D[n].nk2= D[n].nk1 + D[n].nk - 1; 
  if (D[n].ni>0 && D[n].nj>0 && D[n].nk>0)
  {
     D[n].enable = 1;
     // esti L0 and dh
     float L0, dh;
     bdry_cal_abl_len_dh(gd,D[n].ni1,
                            D[n].ni2,
                            D[n].nj1,
                            D[n].nj2,
                            D[n].nk1,
                            D[n].nk2,
                            1,
                            &L0, &dh);

     for (int j=D[n].nj1; j<=D[n].nj2; j++)
     {
        bdry->ablexp_Ey[j] = bdry_ablexp_cal_mask(j - D[n].nj1 + 1,
                                                  in_velocity[1][1], dt,
                                                  D[n].nj, dh);
     }
  }

  // z1
  n += 1;
  D[n].ni = ni - ni_x;
  D[n].ni1= ni1 + bdry->num_of_layers[0][0];
  D[n].ni2= D[n].ni1 + D[n].ni - 1; 
  D[n].nj = nj - nj_y;
  D[n].nj1= nj1 + bdry->num_of_layers[1][0];
  D[n].nj2= D[n].nj1 + D[n].nj - 1; 
  D[n].nk = bdry->num_of_layers[2][0];
  D[n].nk1= nk1;
  D[n].nk2= D[n].nk1 + D[n].nk - 1; 
  if (D[n].ni>0 && D[n].nj>0 && D[n].nk>0)
  {
     D[n].enable = 1;
     // esti L0 and dh
     float L0, dh;
     bdry_cal_abl_len_dh(gd,D[n].ni1,
                            D[n].ni2,
                            D[n].nj1,
                            D[n].nj2,
                            D[n].nk1,
                            D[n].nk2,
                            2,
                            &L0, &dh);

     for (int k=D[n].nk1; k<=D[n].nk2; k++)
     {
        bdry->ablexp_Ez[k] = bdry_ablexp_cal_mask(D[n].nk - (k - D[n].nk1),
                                                  in_velocity[2][0], dt,
                                                  D[n].nk, dh);
     }
  }

  // z2
  n += 1;
  D[n].ni = ni - ni_x;
  D[n].ni1= ni1 + bdry->num_of_layers[0][0];
  D[n].ni2= D[n].ni1 + D[n].ni - 1; 
  D[n].nj = nj - nj_y;
  D[n].nj1= nj1 + bdry->num_of_layers[1][0];
  D[n].nj2= D[n].nj1 + D[n].nj - 1; 
  D[n].nk = bdry->num_of_layers[2][1];
  D[n].nk1= nk2 - D[n].nk + 1;
  D[n].nk2= D[n].nk1 + D[n].nk - 1; 
  if (D[n].ni>0 && D[n].nj>0 && D[n].nk>0)
  {
     D[n].enable = 1;
     // esti L0 and dh
     float L0, dh;
     bdry_cal_abl_len_dh(gd,D[n].ni1,
                            D[n].ni2,
                            D[n].nj1,
                            D[n].nj2,
                            D[n].nk1,
                            D[n].nk2,
                            2,
                            &L0, &dh);

     for (int k=D[n].nk1; k<=D[n].nk2; k++)
     {
        bdry->ablexp_Ez[k] = bdry_ablexp_cal_mask(k - D[n].nk1 + 1,
                                                  in_velocity[2][1], dt,
                                                  D[n].nk, dh);
     }
  }

#ifdef BDRY_T_DEBUG
  fprintf(stderr,"== debug for ablexp at %i %i\n", topoid[0],topoid[1]);
  fprintf(stderr,"--- nx = %d\n", nx);
  for (int i=0; i<nx; i++) {
    fprintf(stderr,"---- Ex[%d] = %g\n", i, bdry->ablexp_Ex[i]);
  }
  fprintf(stderr,"--- ny = %d\n", ny);
  for (int j=0; j<ny; j++) {
    fprintf(stderr,"---- Ey[%d] = %g\n", j, bdry->ablexp_Ey[j]);
  }
  fprintf(stderr,"--- nz = %d\n", nz);
  for (int k=0; k<nz; k++) {
    fprintf(stderr,"---- Ez[%d] = %g\n", k, bdry->ablexp_Ez[k]);
  }

  fprintf(stderr,"--- index for each block\n");
  for (int n=0; n<CONST_NDIM_2; n++)
  {
    fprintf(stderr,"---- block %d enable = %d\n", n, D[n].enable);
    fprintf(stderr,"----- ni1=%d,ni2=%d,ni=%d\n", D[n].ni1, D[n].ni2, D[n].ni);
    fprintf(stderr,"----- nj1=%d,nj2=%d,nj=%d\n", D[n].nj1, D[n].nj2, D[n].nj);
    fprintf(stderr,"----- nk1=%d,nk2=%d,nk=%d\n", D[n].nk1, D[n].nk2, D[n].nk);
  }

  fflush(stderr);
#endif

  return ierr;
}

float
bdry_ablexp_cal_mask(int i, float vel, float dt, int num_lay, float dh)
{
  float len = num_lay * dh;

  int num_step  = (int) (len / vel / dt);

  float total_damp=0.0;
  for (int n=0; n<num_step; n++) {
     total_damp += powf((n*dt*vel)/len, 2.0);
  }

  float alpha = 0.6 / total_damp;

  float mask_val = expf( -alpha * powf((float)i/num_lay, 2.0 ) );

  return mask_val;
}

int
bdry_ablexp_apply(bdry_t *bdry, float *w_end, int ncmp, size_t siz_icmp)
{
  float *Ex = bdry->ablexp_Ex;
  float *Ey = bdry->ablexp_Ey;
  float *Ez = bdry->ablexp_Ez;

  int nx = bdry->nx;
  int ny = bdry->ny;
  int nz = bdry->nz;

  size_t siz_slice = nx * ny;

  bdry_block_t *D = bdry->bdry_blk;
  
  for (int ivar=0; ivar<ncmp; ivar++)
  {
    float *W = w_end + ivar * siz_icmp;

    for (int n=0; n < CONST_NDIM_2; n++)
    {
      if (D[n].enable == 1)
      {
        for (int k = D[n].nk1; k <= D[n].nk2; k++)
        {
          size_t iptr_k = k * siz_slice;
          for (int j = D[n].nj1; j <= D[n].nj2; j++)
          {
            size_t iptr = iptr_k + j * nx + D[n].ni1;
            for (int i = D[n].ni1; i <= D[n].ni2; i++)
            {
              float mask = (Ex[i]<Ey[j]) ? Ex[i] : Ey[j];
              if (mask > Ez[k]) mask = Ez[k];

#ifdef BDRY_T_DEBUG
              fprintf(stderr,"--> n=%d,i=%d,j=%d,k=%d\n", n,i,j,k);
              fprintf(stderr,"----> Ex=%g,Ey=%g,Ez=%g,mask=%g\n",Ex[i],Ey[j],Ez[k],mask);
#endif

              W[iptr] *= mask;

              iptr +=1 ;
            }
          }
        }
      }
    }
  }

  return 0;
}
