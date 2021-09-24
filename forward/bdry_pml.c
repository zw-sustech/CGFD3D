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

void
bdry_pml_set(gdinfo_t *gdinfo,
             gd_t *gd,
             wav_t *wav,
             bdrypml_t *bdrypml,
             int   *neighid, 
             int   in_is_sides[][2],
             int   in_num_layers[][2],
             float in_alpha_max[][2], //
             float in_beta_max[][2], //
             float in_velocity[][2], //
             int verbose)
{
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
      if (bdrypml->is_at_sides  [idim][iside] == 1) {
        bdrypml->is_enable = 1;
      }

    } // iside
  } // idim

  // alloc coef
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      if (bdrypml->is_at_sides[idim][iside] == 1) {
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
      if (bdrypml->is_at_sides[idim][iside] == 0) continue;

      float *A = bdrypml->A[idim][iside];
      float *B = bdrypml->B[idim][iside];
      float *D = bdrypml->D[idim][iside];

      // esti L0 and dh
      float L0, dh;
      bdry_pml_cal_len_dh(gd,bdrypml->ni1[idim][iside],
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

void
bdry_pml_set_stg(gdinfo_t *gdinfo,
                 gd_t *gd,
                 wav_t *wav,
                 bdrypml_t *bdrypml,
                 int   *neighid, 
                 int   in_is_sides[][2],
                 int   in_num_layers[][2],
                 float in_alpha_max[][2], //
                 float in_beta_max[][2], //
                 float in_velocity[][2], //
                 int verbose)
{
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
      if (bdrypml->is_at_sides  [idim][iside] == 1) {
        bdrypml->is_enable = 1;
      }

    } // iside
  } // idim

  // alloc coef
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      if (bdrypml->is_at_sides[idim][iside] == 1) {
        int npoints = bdrypml->num_of_layers[idim][iside] + 1;
        bdrypml->A[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->B[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->D[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->Am[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->Bm[idim][iside] = (float *)malloc( npoints * sizeof(float));
        bdrypml->Dm[idim][iside] = (float *)malloc( npoints * sizeof(float));
      } else {
        bdrypml->A[idim][iside] = NULL;
        bdrypml->B[idim][iside] = NULL;
        bdrypml->D[idim][iside] = NULL;
        bdrypml->Am[idim][iside] = NULL;
        bdrypml->Bm[idim][iside] = NULL;
        bdrypml->Dm[idim][iside] = NULL;
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

      float *A = bdrypml->A[idim][iside];
      float *B = bdrypml->B[idim][iside];
      float *D = bdrypml->D[idim][iside];

      float *Am = bdrypml->Am[idim][iside];
      float *Bm = bdrypml->Bm[idim][iside];
      float *Dm = bdrypml->Dm[idim][iside];

      // esti L0 and dh
      float L0, dh;
      bdry_pml_cal_len_dh(gd,bdrypml->ni1[idim][iside],
                             bdrypml->ni2[idim][iside],
                             bdrypml->nj1[idim][iside],
                             bdrypml->nj2[idim][iside],
                             bdrypml->nk1[idim][iside],
                             bdrypml->nk2[idim][iside],
                             idim,
                             &L0, &dh);

      // adjust L0 to include half length
      L0 += dh / 2.0;

      // para
      int npoints = bdrypml->num_of_layers[idim][iside] + 1;
      float num_lay = npoints - 1 + 0.5; // staggered grid
      float Rpp  = bdry_pml_cal_R(num_lay);
      float dmax = bdry_pml_cal_dmax(L0, in_velocity[idim][iside], Rpp);
      float amax = in_alpha_max[idim][iside];
      float bmax = in_beta_max[idim][iside];

      // from PML-interior to outer side
      for (int n=0; n<npoints; n++)
      {
        float L, Lm;
        int i;

        // convert to grid index from left to right
        if (iside == 0) // x1/y1/z1
        { // first point at middle
          Lm = n * dh;
          L  = (n + 0.5 ) * dh;
          i = npoints - 1 - n;
        }
        else // x2/y2/z2
        { // first point at integer
          L  = n * dh;
          Lm = L + 0.5 * dh;
          i = n; 
        }

        // regular point
        D[i] = bdry_pml_cal_d( L, L0, dmax );
        A[i] = bdry_pml_cal_a( L, L0, amax );
        B[i] = bdry_pml_cal_b( L, L0, bmax );

        // middle point
        Dm[i] = bdry_pml_cal_d( Lm, L0, dmax );
        Am[i] = bdry_pml_cal_a( Lm, L0, amax );
        Bm[i] = bdry_pml_cal_b( Lm, L0, bmax );

        // convert d_x to d_x/beta_x since only d_x/beta_x needed
        D[i] /= B[i];
        // covert ax = a_x + d_x/beta_x 
        A[i] += D[i];
        // covert bx = 1.0/bx 
        B[i] = 1.0 / B[i];

        // for middle point
        Dm[i] /= Bm[i];
        Am[i] += Dm[i];
        Bm[i] = 1.0 / Bm[i];
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
 * esti L and dh along idim damping layers
 */

int
bdry_pml_cal_len_dh(gd_t *gd, 
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
