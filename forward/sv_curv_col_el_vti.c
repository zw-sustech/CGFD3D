/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and macdrp schem
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"
#include "sv_curv_col_el.h"
#include "sv_curv_col_el_vti.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_curv_col_el_vti_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gd_t   *gdinfo,
  gdcurv_metric_t  *metric,
  md_t *md,
  bdry_t    *bdry,
  src_t *src,
  // include different order/stentil
  int num_of_fdx_op, fd_op_t *fdx_op,
  int num_of_fdy_op, fd_op_t *fdy_op,
  int num_of_fdz_op, fd_op_t *fdz_op,
  int fdz_max_len, 
  const int myid, const int verbose)
{
  // local pointer get each vars
  float *restrict Vx    = w_cur + wav->Vx_pos ;
  float *restrict Vy    = w_cur + wav->Vy_pos ;
  float *restrict Vz    = w_cur + wav->Vz_pos ;
  float *restrict Txx   = w_cur + wav->Txx_pos;
  float *restrict Tyy   = w_cur + wav->Tyy_pos;
  float *restrict Tzz   = w_cur + wav->Tzz_pos;
  float *restrict Txz   = w_cur + wav->Txz_pos;
  float *restrict Tyz   = w_cur + wav->Tyz_pos;
  float *restrict Txy   = w_cur + wav->Txy_pos;
  float *restrict hVx   = rhs   + wav->Vx_pos ; 
  float *restrict hVy   = rhs   + wav->Vy_pos ; 
  float *restrict hVz   = rhs   + wav->Vz_pos ; 
  float *restrict hTxx  = rhs   + wav->Txx_pos; 
  float *restrict hTyy  = rhs   + wav->Tyy_pos; 
  float *restrict hTzz  = rhs   + wav->Tzz_pos; 
  float *restrict hTxz  = rhs   + wav->Txz_pos; 
  float *restrict hTyz  = rhs   + wav->Tyz_pos; 
  float *restrict hTxy  = rhs   + wav->Txy_pos; 

  float *restrict xi_x  = metric->xi_x;
  float *restrict xi_y  = metric->xi_y;
  float *restrict xi_z  = metric->xi_z;
  float *restrict et_x  = metric->eta_x;
  float *restrict et_y  = metric->eta_y;
  float *restrict et_z  = metric->eta_z;
  float *restrict zt_x  = metric->zeta_x;
  float *restrict zt_y  = metric->zeta_y;
  float *restrict zt_z  = metric->zeta_z;
  float *restrict jac3d = metric->jac;

  float *restrict c11   = md->c11;
  float *restrict c13   = md->c13;
  float *restrict c33   = md->c33;
  float *restrict c55   = md->c55;
  float *restrict c66   = md->c66;
  float *restrict slw3d = md->rho;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;

  int ni  = gdinfo->ni;
  int nj  = gdinfo->nj;
  int nk  = gdinfo->nk;
  int nx  = gdinfo->nx;
  int ny  = gdinfo->ny;
  int nz  = gdinfo->nz;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;
  size_t siz_volume = gdinfo->siz_volume;

  float *matVx2Vz = bdry->matVx2Vz2;
  float *matVy2Vz = bdry->matVy2Vz2;

  // local fd op
  int              fdx_inn_len;
  int    *restrict fdx_inn_indx;
  float *restrict fdx_inn_coef;
  int              fdy_inn_len;
  int    *restrict fdy_inn_indx;
  float  *restrict fdy_inn_coef;
  int              fdz_inn_len;
  int    *restrict fdz_inn_indx;
  float  *restrict fdz_inn_coef;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op
  fdx_inn_len  = fdx_op[num_of_fdx_op-1].total_len;
  fdx_inn_indx = fdx_op[num_of_fdx_op-1].indx;
  fdx_inn_coef = fdx_op[num_of_fdx_op-1].coef;

  fdy_inn_len  = fdy_op[num_of_fdy_op-1].total_len;
  fdy_inn_indx = fdy_op[num_of_fdy_op-1].indx;
  fdy_inn_coef = fdy_op[num_of_fdy_op-1].coef;

  fdz_inn_len  = fdz_op[num_of_fdz_op-1].total_len;
  fdz_inn_indx = fdz_op[num_of_fdz_op-1].indx;
  fdz_inn_coef = fdz_op[num_of_fdz_op-1].coef;

  // inner points
  sv_curv_col_el_vti_rhs_inner(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                    hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                    c11,c13, c33, c55, c66, slw3d,
                                    ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                    fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                    fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                    fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                    myid, verbose);

  // free, abs, source in turn

  // free surface at z2
  if (bdry->is_sides_free[2][1] == 1)
  {
    // tractiong
    sv_curv_col_el_rhs_timg_z2(Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        jac3d, slw3d,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        myid, verbose);

    // velocity: vlow
    sv_curv_col_el_vti_rhs_vlow_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        c11,c13, c33, c55, c66, slw3d,
                                        matVx2Vz,matVy2Vz,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        num_of_fdz_op,fdz_op,fdz_max_len,
                                        myid, verbose);
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_curv_col_el_vti_rhs_cfspml(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                       xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                       c11,c13, c33, c55, c66, slw3d,
                                       nk2, siz_line,siz_slice,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       bdry,
                                       myid, verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_curv_col_el_rhs_src(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    jac3d, slw3d, 
                                    src,
                                    myid, verbose);
  }

  return;
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

void
sv_curv_col_el_vti_rhs_inner(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict c11d, float *restrict c13d, float *restrict c33d,
    float *restrict c55d, float *restrict c66d, float *restrict slw3d,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int    lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float slw;
  float c11,c13,c33,c55,c66;
  float c12;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  float *restrict Vx_ptr;
  float *restrict Vy_ptr;
  float *restrict Vz_ptr;
  float *restrict Txx_ptr;
  float *restrict Txy_ptr;
  float *restrict Txz_ptr;
  float *restrict Tyy_ptr;
  float *restrict Tzz_ptr;
  float *restrict Tyz_ptr;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }


  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
        Vx_ptr = Vx + iptr;
        Vy_ptr = Vy + iptr;
        Vz_ptr = Vz + iptr;
        Txx_ptr = Txx + iptr;
        Tyy_ptr = Tyy + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;
        Tyz_ptr = Tyz + iptr;
        Txy_ptr = Txy + iptr;

        // Vx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyVx, Vx_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVy, Vy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyVy, Vy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVy, Vy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyVz, Vz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxx, Txx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTxx, Txx_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxx, Txx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyy derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTyy, Tyy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTyy, Tyy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTyy, Tyy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tzz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTzz, Tzz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTzz, Tzz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTxz, Txz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTyz, Tyz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTyz, Tyz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTyz, Tyz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txy derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxy, Txy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTxy, Txy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxy, Txy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];
        etx = et_x[iptr];
        ety = et_y[iptr];
        etz = et_z[iptr];
        ztx = zt_x[iptr];
        zty = zt_y[iptr];
        ztz = zt_z[iptr];

        // medium
        slw = slw3d[iptr];
        c11 = c11d[iptr];
        c13 = c13d[iptr];
        c33 = c33d[iptr];
        c55 = c55d[iptr];
        c66 = c66d[iptr];
        c12 = c11 - 2.0 * c66;

        // moment equation
        hVx[iptr] = slw*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz  
                         +etx*DyTxx + ety*DyTxy + etz*DyTxz 
                         +ztx*DzTxx + zty*DzTxy + ztz*DzTxz );
        hVy[iptr] = slw*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz
                         +etx*DyTxy + ety*DyTyy + etz*DyTyz
                         +ztx*DzTxy + zty*DzTyy + ztz*DzTyz );
        hVz[iptr] = slw*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz 
                         +etx*DyTxz + ety*DyTyz + etz*DyTzz
                         +ztx*DzTxz + zty*DzTyz + ztz*DzTzz );

        // Hooke's equatoin

	      hTxx[iptr] = (c11*xix ) * DxVx + ( c12*xiy ) * DxVy  + ( c13*xiz) * DxVz
                   + (c11*etx ) * DyVx + ( c12*ety ) * DyVy  + ( c13*etz) * DyVz
                   + (c11*ztx ) * DzVx + ( c12*zty ) * DzVy  + ( c13*ztz) * DzVz;
      
        hTyy[iptr] = (c12*xix ) * DxVx + ( c11*xiy ) * DxVy + ( c13*xiz) * DxVz
                   + (c12*etx ) * DyVx + ( c11*ety ) * DyVy + ( c13*etz) * DyVz
                   + (c12*ztx ) * DzVx + ( c11*zty ) * DzVy + ( c13*ztz) * DzVz;
     
        hTzz[iptr] = (c13*xix ) * DxVx + ( c13*xiy ) * DxVy + ( c33*xiz) * DxVz
                   + (c13*etx ) * DyVx + ( c13*ety ) * DyVy + ( c33*etz) * DyVz
                   + (c13*ztx ) * DzVx + ( c13*zty ) * DzVy + ( c33*ztz) * DzVz;
  

        hTyz[iptr] = ( c55*xiz) * DxVy + ( c55*xiy ) * DxVz
                   + ( c55*etz) * DyVy + ( c55*ety ) * DyVz
                   + ( c55*ztz) * DzVy + ( c55*zty ) * DzVz;
  
        hTxz[iptr] = ( c55*xiz) * DxVx + (c55*xix ) * DxVz
                   + ( c55*etz) * DyVx + (c55*etx ) * DyVz
                   + ( c55*ztz) * DzVx + (c55*ztx ) * DzVz;
  
        hTxy[iptr] = ( c66*xiy ) * DxVx + (c66*xix ) * DxVy
                   + ( c66*ety ) * DyVx + (c66*etx ) * DyVy
                   + ( c66*zty ) * DzVx + (c66*ztx ) * DzVy;

        iptr += 1;
      }
    }
  }

  return;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement vlow boundary
 */

void
sv_curv_col_el_vti_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict c11d, float *restrict c13d, float *restrict c33d,
    float *restrict c55d, float *restrict c66d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int    lfdy_shift[fdy_len];

  // allocate max_len because fdz may have different lens
  float  lfdz_coef [fdz_max_len];
  int    lfdz_shift[fdz_max_len];

  // local var
  int i,j,k;
  int n_fd; // loop var for fd
  int fdz_len;

  // local var
  float DxVx,DxVy,DxVz;
  float DyVx,DyVy,DyVz;
  float DzVx,DzVy,DzVz;
  float slw;
  float c11,c13,c33,c55,c66;
  float c12;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }

  // loop near surface layers
  //for (size_t n=0; n < 1; n++)
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // get pos and len for this point
    int  lfdz_len  = fdz_op[n].total_len;
    // point to indx/coef for this point
    int   *p_fdz_indx  = fdz_op[n].indx;
    float *p_fdz_coef  = fdz_op[n].coef;
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_shift[n_fd] = p_fdz_indx[n_fd] * siz_slice;
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd];
    }

    // for index
    size_t iptr_k = k * siz_slice;

    for (j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (i=ni1; i<=ni2; i++)
      {
        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];
        etx = et_x[iptr];
        ety = et_y[iptr];
        etz = et_z[iptr];
        ztx = zt_x[iptr];
        zty = zt_y[iptr];
        ztz = zt_z[iptr];

        // medium
        slw = slw3d[iptr];
        c11 = c11d[iptr];
        c13 = c13d[iptr];
        c33 = c33d[iptr];
        c55 = c55d[iptr];
        c66 = c66d[iptr];
        c12 = c11 - 2.0 * c66;

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVx, Vx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT(DxVy, Vy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVy, Vy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVz, Vz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        if (k==nk2) // at surface, convert
        {
          size_t ij = (i + j * siz_line)*9;
          DzVx = matVx2Vz[ij+3*0+0] * DxVx
               + matVx2Vz[ij+3*0+1] * DxVy
               + matVx2Vz[ij+3*0+2] * DxVz
               + matVy2Vz[ij+3*0+0] * DyVx
               + matVy2Vz[ij+3*0+1] * DyVy
               + matVy2Vz[ij+3*0+2] * DyVz;

          DzVy = matVx2Vz[ij+3*1+0] * DxVx
               + matVx2Vz[ij+3*1+1] * DxVy
               + matVx2Vz[ij+3*1+2] * DxVz
               + matVy2Vz[ij+3*1+0] * DyVx
               + matVy2Vz[ij+3*1+1] * DyVy
               + matVy2Vz[ij+3*1+2] * DyVz;

          DzVz = matVx2Vz[ij+3*2+0] * DxVx
               + matVx2Vz[ij+3*2+1] * DxVy
               + matVx2Vz[ij+3*2+2] * DxVz
               + matVy2Vz[ij+3*2+0] * DyVx
               + matVy2Vz[ij+3*2+1] * DyVy
               + matVy2Vz[ij+3*2+2] * DyVz;
        }
        else // lower than surface, lower order
        {
          M_FD_SHIFT(DzVx, Vx, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVy, Vy, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVz, Vz, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
        }

        // Hooke's equatoin
	      hTxx[iptr] = (c11*xix ) * DxVx + ( c12*xiy ) * DxVy  + ( c13*xiz) * DxVz
                   + (c11*etx ) * DyVx + ( c12*ety ) * DyVy  + ( c13*etz) * DyVz
                   + (c11*ztx ) * DzVx + ( c12*zty ) * DzVy  + ( c13*ztz) * DzVz;
      
        hTyy[iptr] = (c12*xix ) * DxVx + ( c11*xiy ) * DxVy + ( c13*xiz) * DxVz
                   + (c12*etx ) * DyVx + ( c11*ety ) * DyVy + ( c13*etz) * DyVz
                   + (c12*ztx ) * DzVx + ( c11*zty ) * DzVy + ( c13*ztz) * DzVz;
     
        hTzz[iptr] = (c13*xix ) * DxVx + ( c13*xiy ) * DxVy + ( c33*xiz) * DxVz
                   + (c13*etx ) * DyVx + ( c13*ety ) * DyVy + ( c33*etz) * DyVz
                   + (c13*ztx ) * DzVx + ( c13*zty ) * DzVy + ( c33*ztz) * DzVz;
  

        hTyz[iptr] = ( c55*xiz) * DxVy + ( c55*xiy ) * DxVz
                   + ( c55*etz) * DyVy + ( c55*ety ) * DyVz
                   + ( c55*ztz) * DzVy + ( c55*zty ) * DzVz;
  
        hTxz[iptr] = ( c55*xiz) * DxVx + (c55*xix ) * DxVz
                   + ( c55*etz) * DyVx + (c55*etx ) * DyVz
                   + ( c55*ztz) * DzVx + (c55*ztx ) * DzVz;
  
        hTxy[iptr] = ( c66*xiy ) * DxVx + (c66*xix ) * DxVy
                   + ( c66*ety ) * DyVx + (c66*etx ) * DyVy
                   + ( c66*zty ) * DzVx + (c66*ztx ) * DzVy;

        iptr += 1;
      }
    }
  }

  return;
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

void
sv_curv_col_el_vti_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict c11d, float *restrict c13d, float *restrict c33d,
    float *restrict c55d, float *restrict c66d, float *restrict slw3d,
    int nk2, size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    bdry_t    *bdry,
    const int myid, const int verbose)
{

  float *matVx2Vz = bdry->matVx2Vz2;
  float *matVy2Vz = bdry->matVy2Vz2;

  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int    lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // val on point
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float slw;
  float c11,c13,c33,c55,c66;
  float c12;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
  float hVx_rhs,hVy_rhs,hVz_rhs;
  float hTxx_rhs,hTyy_rhs,hTzz_rhs,hTxz_rhs,hTyz_rhs,hTxy_rhs;
  // for free surface
  float Dx_DzVx,Dy_DzVx,Dx_DzVy,Dy_DzVy,Dx_DzVz,Dy_DzVz;

  // local
  int i,j,k;
  int iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_B_minus_1;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      // skip to next face if not cfspml
      if (bdry->is_sides_pml[idim][iside] == 0) continue;

      // get index into local var
      int abs_ni1 = bdry->ni1[idim][iside];
      int abs_ni2 = bdry->ni2[idim][iside];
      int abs_nj1 = bdry->nj1[idim][iside];
      int abs_nj2 = bdry->nj2[idim][iside];
      int abs_nk1 = bdry->nk1[idim][iside];
      int abs_nk2 = bdry->nk2[idim][iside];

#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
    //fprintf(stdout," iface=%d,ni1=%d,ni2=%d,nj1=%d,nj2=%d,nk1=%d,nk2=%d\n",
    //        iface,abs_ni1,abs_ni2,abs_nj1,abs_nj2,abs_nk1,abs_nk2);
    //fflush(stdout);
#endif

      // get coef for this face
      float *restrict ptr_coef_A = bdry->A[idim][iside];
      float *restrict ptr_coef_B = bdry->B[idim][iside];
      float *restrict ptr_coef_D = bdry->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);

      // get pml vars
      float *restrict abs_vars_cur = auxvar->cur;
      float *restrict abs_vars_rhs = auxvar->rhs;

      float *restrict pml_Vx   = abs_vars_cur + auxvar->Vx_pos;
      float *restrict pml_Vy   = abs_vars_cur + auxvar->Vy_pos;
      float *restrict pml_Vz   = abs_vars_cur + auxvar->Vz_pos;
      float *restrict pml_Txx  = abs_vars_cur + auxvar->Txx_pos;
      float *restrict pml_Tyy  = abs_vars_cur + auxvar->Tyy_pos;
      float *restrict pml_Tzz  = abs_vars_cur + auxvar->Tzz_pos;
      float *restrict pml_Txz  = abs_vars_cur + auxvar->Txz_pos;
      float *restrict pml_Tyz  = abs_vars_cur + auxvar->Tyz_pos;
      float *restrict pml_Txy  = abs_vars_cur + auxvar->Txy_pos;

      float *restrict pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
      float *restrict pml_hVy  = abs_vars_rhs + auxvar->Vy_pos;
      float *restrict pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
      float *restrict pml_hTxx = abs_vars_rhs + auxvar->Txx_pos;
      float *restrict pml_hTyy = abs_vars_rhs + auxvar->Tyy_pos;
      float *restrict pml_hTzz = abs_vars_rhs + auxvar->Tzz_pos;
      float *restrict pml_hTxz = abs_vars_rhs + auxvar->Txz_pos;
      float *restrict pml_hTyz = abs_vars_rhs + auxvar->Tyz_pos;
      float *restrict pml_hTxy = abs_vars_rhs + auxvar->Txy_pos;

      // for each dim
      if (idim == 0 ) // x direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_slice;
          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // pml coefs
              int abs_i = i - abs_ni1;
              coef_D = ptr_coef_D[abs_i];
              coef_A = ptr_coef_A[abs_i];
              coef_B = ptr_coef_B[abs_i];
              coef_B_minus_1 = coef_B - 1.0;

              // metric
              xix = xi_x[iptr];
              xiy = xi_y[iptr];
              xiz = xi_z[iptr];

              // medium
              slw = slw3d[iptr];
              c11 = c11d[iptr];
              c13 = c13d[iptr];
              c33 = c33d[iptr];
              c55 = c55d[iptr];
              c66 = c66d[iptr];
              c12 = c11 - 2.0 * c66;

              // xi derivatives
              M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVy , Vy , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTyy, Tyy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTzz, Tzz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTyz, Tyz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxTxy, Txy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( xix*DxTxx + xiy*DxTxy + xiz*DxTxz );
               hVy_rhs = slw * ( xix*DxTxy + xiy*DxTyy + xiz*DxTyz );
               hVz_rhs = slw * ( xix*DxTxz + xiy*DxTyz + xiz*DxTzz );
              hTxx_rhs = (c11*xix)*DxVx + (c12*xiy)*DxVy + (c13*xiz)*DxVz; 
              hTyy_rhs = (c12*xix)*DxVx + (c11*xiy)*DxVy + (c13*xiz)*DxVz;
              hTzz_rhs = (c13*xix)*DxVx + (c13*xiy)*DxVy + (c33*xiz)*DxVz;
              hTyz_rhs = (c55*xiz)*DxVy + (c55*xiy)*DxVz;
              hTxz_rhs = (c55*xiz)*DxVx + (c55*xix)*DxVz;
              hTxy_rhs = (c66*xiy)*DxVx + (c66*xix)*DxVy;

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
              hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
              hTyz[iptr] += coef_B_minus_1 * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
              hTxy[iptr] += coef_B_minus_1 * hTxy_rhs - coef_B * pml_Txy[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
              pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
              pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
              pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
              pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
              pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

              // add contributions from free surface condition
              //  not consider timg because conflict with main cfspml,
              //     need to revise in the future if required
              if (bdry->is_sides_free[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i + j * siz_line)*9;
                Dx_DzVx = matVx2Vz[ij+3*0+0] * DxVx
                        + matVx2Vz[ij+3*0+1] * DxVy
                        + matVx2Vz[ij+3*0+2] * DxVz;

                Dx_DzVy = matVx2Vz[ij+3*1+0] * DxVx
                        + matVx2Vz[ij+3*1+1] * DxVy
                        + matVx2Vz[ij+3*1+2] * DxVz;

                Dx_DzVz = matVx2Vz[ij+3*2+0] * DxVx
                        + matVx2Vz[ij+3*2+1] * DxVy
                        + matVx2Vz[ij+3*2+2] * DxVz;

                // metric
                ztx = zt_x[iptr];
                zty = zt_y[iptr];
                ztz = zt_z[iptr];

                // keep xi derivative terms, including free surface convered
                hTxx_rhs = (c11*ztx)*Dx_DzVx + (c12*zty)*Dx_DzVy + (c13*ztz)*Dx_DzVz; 
                hTyy_rhs = (c12*ztx)*Dx_DzVx + (c11*zty)*Dx_DzVy + (c13*ztz)*Dx_DzVz;
                hTzz_rhs = (c13*ztx)*Dx_DzVx + (c13*zty)*Dx_DzVy + (c33*ztz)*Dx_DzVz;
                hTyz_rhs = (c55*ztz)*Dx_DzVy + (c55*zty)*Dx_DzVz;
                hTxz_rhs = (c55*ztz)*Dx_DzVx + (c55*ztx)*Dx_DzVz;
                hTxy_rhs = (c66*zty)*Dx_DzVx + (c66*ztx)*Dx_DzVy;

                // make corr to Hooke's equatoin
                hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
                hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs;
                hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
                hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;
                hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs;
                hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs;

                // aux var
                //   a1 = alpha + d / beta, dealt in abs_set_cfspml
                pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
                pml_hTyy[iptr_a] += coef_D * hTyy_rhs;
                pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
                pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
                pml_hTyz[iptr_a] += coef_D * hTyz_rhs;
                pml_hTxy[iptr_a] += coef_D * hTxy_rhs;
              } // if nk2

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else if (idim == 1) // y direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_slice;
          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;

            // pml coefs
            int abs_j = j - abs_nj1;
            coef_D = ptr_coef_D[abs_j];
            coef_A = ptr_coef_A[abs_j];
            coef_B = ptr_coef_B[abs_j];
            coef_B_minus_1 = coef_B - 1.0;

            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // metric
              etx = et_x[iptr];
              ety = et_y[iptr];
              etz = et_z[iptr];

              // medium
              slw = slw3d[iptr];
              c11 = c11d[iptr];
              c13 = c13d[iptr];
              c33 = c33d[iptr];
              c55 = c55d[iptr];
              c66 = c66d[iptr];
              c12 = c11 - 2.0 * c66;

              // et derivatives
              M_FD_SHIFT(DyVx , Vx , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyVy , Vy , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyVz , Vz , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyTxx, Txx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyTyy, Tyy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyTzz, Tzz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyTxz, Txz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyTyz, Tyz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyTxy, Txy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( etx*DyTxx + ety*DyTxy + etz*DyTxz );
               hVy_rhs = slw * ( etx*DyTxy + ety*DyTyy + etz*DyTyz );
               hVz_rhs = slw * ( etx*DyTxz + ety*DyTyz + etz*DyTzz );
              hTxx_rhs = (c11*etx)*DyVx + (c12*ety)*DyVy + (c13*etz)*DyVz; 
              hTyy_rhs = (c12*etx)*DyVx + (c11*ety)*DyVy + (c13*etz)*DyVz;
              hTzz_rhs = (c13*etx)*DyVx + (c13*ety)*DyVy + (c33*etz)*DyVz;
              hTyz_rhs = (c55*etz)*DyVy + (c55*ety)*DyVz;
              hTxz_rhs = (c55*etz)*DyVx + (c55*etx)*DyVz;
              hTxy_rhs = (c66*ety)*DyVx + (c66*etx)*DyVy;

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
              hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
              hTyz[iptr] += coef_B_minus_1 * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
              hTxy[iptr] += coef_B_minus_1 * hTxy_rhs - coef_B * pml_Txy[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
              pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
              pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
              pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
              pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
              pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

              // add contributions from free surface condition
              if (bdry->is_sides_free[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i + j * siz_line)*9;
                Dy_DzVx = matVy2Vz[ij+3*0+0] * DyVx
                        + matVy2Vz[ij+3*0+1] * DyVy
                        + matVy2Vz[ij+3*0+2] * DyVz;

                Dy_DzVy = matVy2Vz[ij+3*1+0] * DyVx
                        + matVy2Vz[ij+3*1+1] * DyVy
                        + matVy2Vz[ij+3*1+2] * DyVz;

                Dy_DzVz = matVy2Vz[ij+3*2+0] * DyVx
                        + matVy2Vz[ij+3*2+1] * DyVy
                        + matVy2Vz[ij+3*2+2] * DyVz;

                // metric
                ztx = zt_x[iptr];
                zty = zt_y[iptr];
                ztz = zt_z[iptr];

                // keep eta derivative terms, including free surface convered
                hTxx_rhs = (c11*ztx)*Dy_DzVx + (c12*zty)*Dy_DzVy + (c13*ztz)*Dy_DzVz; 
                hTyy_rhs = (c12*ztx)*Dy_DzVx + (c11*zty)*Dy_DzVy + (c13*ztz)*Dy_DzVz;
                hTzz_rhs = (c13*ztx)*Dy_DzVx + (c13*zty)*Dy_DzVy + (c33*ztz)*Dy_DzVz;
                hTyz_rhs = (c55*ztz)*Dy_DzVy + (c55*zty)*Dy_DzVz;
                hTxz_rhs = (c55*ztz)*Dy_DzVx + (c55*ztx)*Dy_DzVz;
                hTxy_rhs = (c66*zty)*Dy_DzVx + (c66*ztx)*Dy_DzVy;

                // make corr to Hooke's equatoin
                hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
                hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs;
                hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
                hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;
                hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs;
                hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs;

                // aux var
                //   a1 = alpha + d / beta, dealt in abs_set_cfspml
                pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
                pml_hTyy[iptr_a] += coef_D * hTyy_rhs;
                pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
                pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
                pml_hTyz[iptr_a] += coef_D * hTyz_rhs;
                pml_hTxy[iptr_a] += coef_D * hTxy_rhs;
              }

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else // z direction
      {
        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_slice;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_B_minus_1 = coef_B - 1.0;

          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // metric
              ztx = zt_x[iptr];
              zty = zt_y[iptr];
              ztz = zt_z[iptr];

              // medium
              slw = slw3d[iptr];
              c11 = c11d[iptr];
              c13 = c13d[iptr];
              c33 = c33d[iptr];
              c55 = c55d[iptr];
              c66 = c66d[iptr];
              c12 = c11 - 2.0 * c66;

              // zt derivatives
              M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVy , Vy , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxx, Txx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTyy, Tyy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTyz, Tyz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzTxy, Txy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = slw * ( ztx*DzTxx + zty*DzTxy + ztz*DzTxz );
               hVy_rhs = slw * ( ztx*DzTxy + zty*DzTyy + ztz*DzTyz );
               hVz_rhs = slw * ( ztx*DzTxz + zty*DzTyz + ztz*DzTzz );
              hTxx_rhs = (c11*ztx)*DzVx + (c12*zty)*DzVy + (c13*ztz)*DzVz; 
              hTyy_rhs = (c12*ztx)*DzVx + (c11*zty)*DzVy + (c13*ztz)*DzVz;
              hTzz_rhs = (c13*ztx)*DzVx + (c13*zty)*DzVy + (c33*ztz)*DzVz;
              hTyz_rhs = (c55*ztz)*DzVy + (c55*zty)*DzVz;
              hTxz_rhs = (c55*ztz)*DzVx + (c55*ztx)*DzVz;
              hTxy_rhs = (c66*zty)*DzVx + (c66*ztx)*DzVy;

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
              hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
              hTyz[iptr] += coef_B_minus_1 * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
              hTxy[iptr] += coef_B_minus_1 * hTxy_rhs - coef_B * pml_Txy[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
              pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
              pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
              pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
              pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
              pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      } // if which dim
    } // iside
  } // idim

  return;
}

/*******************************************************************************
 * free surface coef
 * converted matrix for velocity gradient
 *  only implement z2 (top) right now
 ******************************************************************************/

int
sv_curv_col_el_vti_dvh2dvz(gd_t        *gdinfo,
                                   gdcurv_metric_t *metric,
                                   md_t       *md,
                                   bdry_t      *bdryfree,
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

  float *restrict c11d = md->c11;
  float *restrict c13d = md->c13;
  float *restrict c33d = md->c33;
  float *restrict c55d = md->c55;
  float *restrict c66d = md->c66;

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;

  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];

  float c11,c13,c33,c55,c66;
  float c12;
  float xix, xiy ,xiz, etx, ety, etz, ztx, zty, ztz;
 
  int k = nk2;

  for (size_t j = nj1; j <= nj2; j++)
  {
    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + j * siz_line + k * siz_slice;

      xix = xi_x[iptr];
      xiy = xi_y[iptr];
      xiz = xi_z[iptr];
      etx = et_x[iptr];
      ety = et_y[iptr];
      etz = et_z[iptr];
      ztx = zt_x[iptr];
      zty = zt_y[iptr];
      ztz = zt_z[iptr];
      
      c11 = c11d[iptr];
      c13 = c13d[iptr];
      c33 = c33d[iptr];
      c55 = c55d[iptr];
      c66 = c66d[iptr];
      c12 = c11 - 2.0 * c66;

      // first dim: irow; sec dim: jcol, as Fortran code
      A[0][0] = (c11*ztx)*ztx + (c66*zty)*zty + (c55*ztz)*ztz;
      A[0][1] = (c12*zty)*ztx + (c66*ztx)*zty;
      A[0][2] = (c13*ztz)*ztx + (c55*ztx)*ztz; 
      A[1][0] = (c66*zty)*ztx + (c12*ztx)*zty; 
      A[1][1] = (c66*ztx)*ztx + (c11*zty)*zty + (c55*ztz)*ztz; 
      A[1][2] = (c13*ztz)*zty + (c55*zty)*ztz;
      A[2][0] = (c55*ztz)*ztx + (c13*ztx)*ztz;
      A[2][1] = (c55*ztz)*zty + (c13*zty)*ztz;
      A[2][2] = (c55*ztx)*ztx + (c55*zty)*zty + (c33*ztz)*ztz; 
      fdlib_math_invert3x3(A);
                                                       
      B[0][0] = (c11*xix)*ztx + (c66*xiy)*zty + (c55*xiz)*ztz;
      B[0][1] = (c12*xiy)*ztx + (c66*xix)*zty;
      B[0][2] = (c13*xiz)*ztx + (c55*xix)*ztz; 
      B[1][0] = (c66*xiy)*ztx + (c12*xix)*zty; 
      B[1][1] = (c66*xix)*ztx + (c11*xiy)*zty + (c55*xiz)*ztz; 
      B[1][2] = (c13*xiz)*zty + (c55*xiy)*ztz;
      B[2][0] = (c55*xiz)*ztx + (c13*xix)*ztz;
      B[2][1] = (c55*xiz)*zty + (c13*xiy)*ztz;
      B[2][2] = (c55*xix)*ztx + (c55*xiy)*zty + (c33*xiz)*ztz; 
       
      C[0][0] = (c11*etx)*ztx + (c66*ety)*zty + (c55*etz)*ztz;
      C[0][1] = (c12*ety)*ztx + (c66*etx)*zty;
      C[0][2] = (c13*etz)*ztx + (c55*etx)*ztz; 
      C[1][0] = (c66*ety)*ztx + (c12*etx)*zty; 
      C[1][1] = (c66*etx)*ztx + (c11*ety)*zty + (c55*etz)*ztz; 
      C[1][2] = (c13*etz)*zty + (c55*ety)*ztz;
      C[2][0] = (c55*etz)*ztx + (c13*etx)*ztz;
      C[2][1] = (c55*etz)*zty + (c13*ety)*ztz;
      C[2][2] = (c55*etx)*ztx + (c55*ety)*zty + (c33*etz)*ztz; 
      fdlib_math_matmul3x3(A, B, AB);
      fdlib_math_matmul3x3(A, C, AC);

      size_t ij = (j * siz_line + i) * 9;

      // save into mat
      for(int irow = 0; irow < 3; irow++)
        for(int jcol = 0; jcol < 3; jcol++){
          matVx2Vz[ij + irow*3 + jcol] = -1.0f * AB[irow][jcol];
          matVy2Vz[ij + irow*3 + jcol] = -1.0f * AC[irow][jcol];
        }
    }
  }

  return ierr;
}

