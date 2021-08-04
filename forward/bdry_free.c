/*
 * prep funcs for free surface condition
 *   implementation of the condition is inside sv_
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "bdry_free.h"

/*
 * converted matrix for velocity gradient
 *  only implement z2 (top) right now
 */

int
bdry_free_set(gdinfo_t        *gdinfo,
              gdcurv_metric_t *metric,
              mdeliso_t       *mdeliso,
              bdryfree_t      *bdryfree,
              int   *neighid, 
              int   in_is_sides[][2],
              const int verbose)
{
  int ierr = 0;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nx  = gdinfo->nx;
  int ny  = gdinfo->ny;
  int nz  = gdinfo->nz;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;
  size_t siz_volume = gdinfo->siz_icmp;

  // point to each var
  float *restrict xi_x = metric->xi_x;
  float *restrict xi_y = metric->xi_y;
  float *restrict xi_z = metric->xi_z;
  float *restrict et_x = metric->eta_x;
  float *restrict et_y = metric->eta_y;
  float *restrict et_z = metric->eta_z;
  float *restrict zt_x = metric->zeta_x;
  float *restrict zt_y = metric->zeta_y;
  float *restrict zt_z = metric->zeta_z;

  float *restrict lam3d = mdeliso->lambda;
  float *restrict  mu3d = mdeliso->mu;

  float e11, e12, e13, e21, e22, e23, e31, e32, e33;
  float lam2mu, lam, mu;
  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];

  // default disable
  bdryfree->is_enable = 0;

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int ind_1d = iside + idim * CONST_NDIM;

      bdryfree->is_at_sides  [idim][iside] = in_is_sides[idim][iside];

      // reset 0 if not mpi boundary
      if (neighid[ind_1d] != MPI_PROC_NULL)
      {
        bdryfree->is_at_sides  [idim][iside] = 0;
      }

      // enable if any side valid
      if (bdryfree->is_at_sides  [idim][iside] == 1) {
        bdryfree->is_enable = 1;
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

  int k = nk2;

  for (size_t j = nj1; j <= nj2; j++)
  {
    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + j * siz_line + k * siz_slice;

      e11 = xi_x[iptr];
      e12 = xi_y[iptr];
      e13 = xi_z[iptr];
      e21 = et_x[iptr];
      e22 = et_y[iptr];
      e23 = et_z[iptr];
      e31 = zt_x[iptr];
      e32 = zt_y[iptr];
      e33 = zt_z[iptr];

      lam    = lam3d[iptr];
      mu     =  mu3d[iptr];
      lam2mu = lam + 2.0f * mu;

      // first dim: irow; sec dim: jcol, as Fortran code
      A[0][0] = lam2mu*e31*e31 + mu*(e32*e32+e33*e33);
      A[0][1] = lam*e31*e32 + mu*e32*e31;
      A[0][2] = lam*e31*e33 + mu*e33*e31;
      A[1][0] = lam*e32*e31 + mu*e31*e32;
      A[1][1] = lam2mu*e32*e32 + mu*(e31*e31+e33*e33);
      A[1][2] = lam*e32*e33 + mu*e33*e32;
      A[2][0] = lam*e33*e31 + mu*e31*e33;
      A[2][1] = lam*e33*e32 + mu*e32*e33;
      A[2][2] = lam2mu*e33*e33 + mu*(e31*e31+e32*e32);
      fdlib_math_invert3x3(A);

      B[0][0] = -lam2mu*e31*e11 - mu*(e32*e12+e33*e13);
      B[0][1] = -lam*e31*e12 - mu*e32*e11;
      B[0][2] = -lam*e31*e13 - mu*e33*e11;
      B[1][0] = -lam*e32*e11 - mu*e31*e12;
      B[1][1] = -lam2mu*e32*e12 - mu*(e31*e11+e33*e13);
      B[1][2] = -lam*e32*e13 - mu*e33*e12;
      B[2][0] = -lam*e33*e11 - mu*e31*e13;
      B[2][1] = -lam*e33*e12 - mu*e32*e13;
      B[2][2] = -lam2mu*e33*e13 - mu*(e31*e11+e32*e12);

      C[0][0] = -lam2mu*e31*e21 - mu*(e32*e22+e33*e23);
      C[0][1] = -lam*e31*e22 - mu*e32*e21;
      C[0][2] = -lam*e31*e23 - mu*e33*e21;
      C[1][0] = -lam*e32*e21 - mu*e31*e22;
      C[1][1] = -lam2mu*e32*e22 - mu*(e31*e21+e33*e23);
      C[1][2] = -lam*e32*e23 - mu*e33*e22;
      C[2][0] = -lam*e33*e21 - mu*e31*e23;
      C[2][1] = -lam*e33*e22 - mu*e32*e23;
      C[2][2] = -lam2mu*e33*e23 - mu*(e31*e21+e32*e22);

      fdlib_math_matmul3x3(A, B, AB);
      fdlib_math_matmul3x3(A, C, AC);

      size_t ij = (j * siz_line + i) * 9;

      // save into mat
      for(int irow = 0; irow < 3; irow++)
        for(int jcol = 0; jcol < 3; jcol++){
          matVx2Vz[ij + irow*3 + jcol] = AB[irow][jcol];
          matVy2Vz[ij + irow*3 + jcol] = AC[irow][jcol];
        }
    }
  }

  bdryfree->matVx2Vz2 = matVx2Vz;
  bdryfree->matVy2Vz2 = matVy2Vz;

  return ierr;
}

