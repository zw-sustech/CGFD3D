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
#include "sv_curv_col_el_iso.h"
#include "sv_curv_col_el_aniso.h"

//#define SV_EQ1ST_CURV_COLGRD_ISO_DEBUG

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_curv_col_el_aniso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gdinfo_t   *gdinfo,
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
  float *restrict c12   = md->c12;
  float *restrict c13   = md->c13;
  float *restrict c14   = md->c14;
  float *restrict c15   = md->c15;
  float *restrict c16   = md->c16;
  float *restrict c22   = md->c22;
  float *restrict c23   = md->c23;
  float *restrict c24   = md->c24;
  float *restrict c25   = md->c25;
  float *restrict c26   = md->c26;
  float *restrict c33   = md->c33;
  float *restrict c34   = md->c34;
  float *restrict c35   = md->c35;
  float *restrict c36   = md->c36;
  float *restrict c44   = md->c44;
  float *restrict c45   = md->c45;
  float *restrict c46   = md->c46;
  float *restrict c55   = md->c55;
  float *restrict c56   = md->c56;
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
  sv_curv_col_el_aniso_rhs_inner(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                    hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                    c11,c12,c13,c14,c15,c16,
                                        c22,c23,c24,c25,c26,
                                            c33,c34,c35,c36,
                                                c44,c45,c46,
                                                    c55,c56,
                                                        c66, slw3d,
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
    sv_curv_col_el_aniso_rhs_timg_z2(Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        jac3d, slw3d,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        myid, verbose);

    // velocity: vlow
    sv_curv_col_el_aniso_rhs_vlow_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        c11,c12,c13,c14,c15,c16,
                                            c22,c23,c24,c25,c26,
                                                c33,c34,c35,c36,
                                                    c44,c45,c46,
                                                        c55,c56,
                                                            c66, slw3d,
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
    sv_curv_col_el_aniso_rhs_cfspml(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                       xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                       c11,c12,c13,c14,c15,c16,
                                           c22,c23,c24,c25,c26,
                                               c33,c34,c35,c36,
                                                   c44,c45,c46,
                                                       c55,c56,
                                                           c66, slw3d,
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
    sv_curv_col_el_iso_rhs_src(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
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
sv_curv_col_el_aniso_rhs_inner(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict c11d, float *restrict c12d, float *restrict c13d,
    float *restrict c14d, float *restrict c15d, float *restrict c16d,
                          float *restrict c22d, float *restrict c23d,
    float *restrict c24d, float *restrict c25d, float *restrict c26d,
                                                float *restrict c33d,
    float *restrict c34d, float *restrict c35d, float *restrict c36d,
    float *restrict c44d, float *restrict c45d, float *restrict c46d,
                          float *restrict c55d, float *restrict c56d,
                                                float *restrict c66d,
    float *restrict slw3d,
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
  float c11,c12,c13,c14,c15,c16;
  float     c22,c23,c24,c25,c26;
  float         c33,c34,c35,c36;
  float             c44,c45,c46;
  float                 c55,c56;
  float                     c66;
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

#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
  /*
  for (int i=0; i < fdx_len; i++) {
    fprintf(stdout," %d", fdx_indx [i]);
  }
  fprintf(stdout,"\n");
  for (int i=0; i < fdx_len; i++) {
    fprintf(stdout," %f", fdx_coef [i]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);

  for (int j=0; j < fdy_len; j++) {
    fprintf(stdout," %d", fdy_indx [j]);
  }
  fprintf(stdout,"\n");
  for (int j=0; j < fdy_len; j++) {
    fprintf(stdout," %f", fdy_coef [j]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);

  for (int k=0; k < fdz_len; k++) {
    fprintf(stdout," %d", fdz_indx [k]);
  }
  fprintf(stdout,"\n");
  for (int k=0; k < fdz_len; k++) {
    fprintf(stdout," %f", fdz_coef [k]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);
  */
#endif

  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
    if (myid==0 && verbose>900) fprintf(stdout,"----> nk1=%d,nk2=%d,k=%d\n", nk1,nk2,k);
#endif
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
      if (myid==0 && verbose>910) fprintf(stdout,"-----> nj1=%d,nj2=%d,j=%d\n", nj1,nj2,j);
#endif
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
        if (myid==0 && verbose>920) fprintf(stdout,"-----> ni1=%d,ni2=%d,i=%d\n", ni1,ni2,i);
#endif
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
        c12 = c12d[iptr];
        c13 = c13d[iptr];
        c14 = c14d[iptr];
        c15 = c15d[iptr];
        c16 = c16d[iptr];
        c22 = c22d[iptr];
        c23 = c23d[iptr];
        c24 = c24d[iptr];
        c25 = c25d[iptr];
        c26 = c26d[iptr];
        c33 = c33d[iptr];
        c34 = c34d[iptr];
        c35 = c35d[iptr];
        c36 = c36d[iptr];
        c44 = c44d[iptr];
        c45 = c45d[iptr];
        c46 = c46d[iptr];
        c55 = c55d[iptr];
        c56 = c56d[iptr];
        c66 = c66d[iptr];

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

	      hTxx[iptr] = (c11*xix + c16*xiy + c15*xiz) * DxVx + (c16*xix + c12*xiy + c14*xiz) * DxVy + (c15*xix + c14*xiy + c13*xiz) * DxVz
                   + (c11*etx + c16*ety + c15*etz) * DyVx + (c16*etx + c12*ety + c14*etz) * DyVy + (c15*etx + c14*ety + c13*etz) * DyVz
                   + (c11*ztx + c16*zty + c15*ztz) * DzVx + (c16*ztx + c12*zty + c14*ztz) * DzVy + (c15*ztx + c14*zty + c13*ztz) * DzVz;
      
        hTyy[iptr] = (c12*xix + c26*xiy + c25*xiz) * DxVx + (c26*xix + c22*xiy + c24*xiz) * DxVy + (c25*xix + c24*xiy + c23*xiz) * DxVz
                   + (c12*etx + c26*ety + c25*etz) * DyVx + (c26*etx + c22*ety + c24*etz) * DyVy + (c25*etx + c24*ety + c23*etz) * DyVz
                   + (c12*ztx + c26*zty + c25*ztz) * DzVx + (c26*ztx + c22*zty + c24*ztz) * DzVy + (c25*ztx + c24*zty + c23*ztz) * DzVz;
     
        hTzz[iptr] = (c13*xix + c36*xiy + c35*xiz) * DxVx + (c36*xix + c23*xiy + c34*xiz) * DxVy + (c35*xix + c34*xiy + c33*xiz) * DxVz
                   + (c13*etx + c36*ety + c35*etz) * DyVx + (c36*etx + c23*ety + c34*etz) * DyVy + (c35*etx + c34*ety + c33*etz) * DyVz
                   + (c13*ztx + c36*zty + c35*ztz) * DzVx + (c36*ztx + c23*zty + c34*ztz) * DzVy + (c35*ztx + c34*zty + c33*ztz) * DzVz;
  

        hTyz[iptr] = (c14*xix + c46*xiy + c45*xiz) * DxVx + (c46*xix + c24*xiy + c44*xiz) * DxVy + (c45*xix + c44*xiy + c34*xiz) * DxVz
                   + (c14*etx + c46*ety + c45*etz) * DyVx + (c46*etx + c24*ety + c44*etz) * DyVy + (c45*etx + c44*ety + c34*etz) * DyVz
                   + (c14*ztx + c46*zty + c45*ztz) * DzVx + (c46*ztx + c24*zty + c44*ztz) * DzVy + (c45*ztx + c44*zty + c34*ztz) * DzVz;
  
        hTxz[iptr] = (c15*xix + c56*xiy + c55*xiz) * DxVx + (c56*xix + c25*xiy + c45*xiz) * DxVy + (c55*xix + c45*xiy + c35*xiz) * DxVz
                   + (c15*etx + c56*ety + c55*etz) * DyVx + (c56*etx + c25*ety + c45*etz) * DyVy + (c55*etx + c45*ety + c35*etz) * DyVz
                   + (c15*ztx + c56*zty + c55*ztz) * DzVx + (c56*ztx + c25*zty + c45*ztz) * DzVy + (c55*ztx + c45*zty + c35*ztz) * DzVz;
  
        hTxy[iptr] = (c16*xix + c66*xiy + c56*xiz) * DxVx + (c66*xix + c26*xiy + c46*xiz) * DxVy + (c56*xix + c46*xiy + c36*xiz) * DxVz
                   + (c16*etx + c66*ety + c56*etz) * DyVx + (c66*etx + c26*ety + c46*etz) * DyVy + (c56*etx + c46*ety + c36*etz) * DyVz
                   + (c16*ztx + c66*zty + c56*ztz) * DzVx + (c66*ztx + c26*zty + c46*ztz) * DzVy + (c56*ztx + c46*zty + c36*ztz) * DzVz;

        iptr += 1;
      }
    }
  }
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_curv_col_el_aniso_rhs_timg_z2(
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTx,DyTy,DzTz;
  float slwjac;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // to save traction and other two dir force var
  float vecxi[fdx_len];
  float vecet[fdy_len];
  float veczt[fdz_len];
  int n, iptr4vec;

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

  // last indx, free surface force Tx/Ty/Tz to 0 in cal
  int k_min = nk2 - fdz_indx[fdz_len-1];

  // point affected by timg
  for (size_t k=k_min; k <= nk2; k++)
  {
    // k corresponding to 0 index of the fd op

    // index of free surface
    int n_free = nk2 - k - fdz_indx[0]; // first indx is negative

    // for 1d index
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
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

          // slowness and jac
          slwjac = slw3d[iptr] / jac3d[iptr];

          //
          // for hVx
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txx[iptr4vec]
                                          + xi_y[iptr4vec] * Txy[iptr4vec]
                                          + xi_z[iptr4vec] * Txz[iptr4vec] );
          }
          for (n=0; n<fdy_len; n++) {
            iptr4vec = iptr + fdy_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txx[iptr4vec]
                                          + et_y[iptr4vec] * Txy[iptr4vec]
                                          + et_z[iptr4vec] * Txz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                          + zt_y[iptr4vec] * Txy[iptr4vec]
                                          + zt_z[iptr4vec] * Txz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++)
          {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            //int n_img = n_free-(n-n_free);
            iptr4vec = iptr + n_img * siz_slice;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                           + zt_y[iptr4vec] * Txy[iptr4vec]
                                           + zt_z[iptr4vec] * Txz[iptr4vec] );
            //veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVx[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
            /*
            if (hVx[iptr] > 1.0e-25) {
              fprintf(stderr, "WARNING: nonzero value for zero input:\n");
              fprintf(stderr, "  i=%d,j=%d,k=%d: DxTx=%f, DyTy=%f, DzTz=%f\n",
                      i,j,k,DxTx,DyTy,DzTz);
              fprintf(stderr, "  %e,%e,%e,%e,%e\n",
                                  veczt[0],
                                  veczt[1],
                                  veczt[2],
                                  veczt[3],
                                  veczt[4]);
              fprintf(stderr, "  n_free=%d,siz_slice=%d,fdz_len=%d,fdz_indx=%d,%d,%d,%d,%d\n",
                          n_free,siz_slice,fdz_len,
                          fdz_indx[0],
                          fdz_indx[1],
                          fdz_indx[2],
                          fdz_indx[3],
                          fdz_indx[4]);

              for (n=0; n<n_free; n++) {
                iptr4vec = iptr + fdz_indx[n] * siz_slice;
                fprintf(stderr,"   n=%d,jac=%e,zt_x=%e,Txx=%e,zt_y=%e,Txy=%e,zt_z=%e,Txz=%e\n",
                           jac3d[iptr4vec] ,    zt_x[iptr4vec] , Txx[iptr4vec],
                                                zt_y[iptr4vec] , Txy[iptr4vec],
                                                zt_z[iptr4vec] , Txz[iptr4vec]);
              }

              fflush(stdout);
              exit(2);
            }
            */
#endif

          //
          // for hVy
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txy[iptr4vec]
                                          + xi_y[iptr4vec] * Tyy[iptr4vec]
                                          + xi_z[iptr4vec] * Tyz[iptr4vec] );
          }
          for (n=0; n<fdy_len; n++) {
            iptr4vec = iptr + fdy_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txy[iptr4vec]
                                          + et_y[iptr4vec] * Tyy[iptr4vec]
                                          + et_z[iptr4vec] * Tyz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txy[iptr4vec]
                                          + zt_y[iptr4vec] * Tyy[iptr4vec]
                                          + zt_z[iptr4vec] * Tyz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            iptr4vec = iptr + n_img * siz_slice;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txy[iptr4vec]
                                           + zt_y[iptr4vec] * Tyy[iptr4vec]
                                           + zt_z[iptr4vec] * Tyz[iptr4vec] );
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVy[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;

          //
          // for hVz
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txz[iptr4vec]
                                          + xi_y[iptr4vec] * Tyz[iptr4vec]
                                          + xi_z[iptr4vec] * Tzz[iptr4vec] );
          }
          for (n=0; n<fdy_len; n++) {
            iptr4vec = iptr + fdy_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txz[iptr4vec]
                                          + et_y[iptr4vec] * Tyz[iptr4vec]
                                          + et_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                          + zt_y[iptr4vec] * Tyz[iptr4vec]
                                          + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // at surface -> set to 0
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            iptr4vec = iptr + n_img * siz_slice;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                           + zt_y[iptr4vec] * Tyz[iptr4vec]
                                           + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // for hVx 
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVz[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;

          // next
          iptr += 1;
      }
    }
  }
}

/*
 * implement vlow boundary
 */

void
sv_curv_col_el_aniso_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict c11d, float *restrict c12d, float *restrict c13d,
    float *restrict c14d, float *restrict c15d, float *restrict c16d,
                          float *restrict c22d, float *restrict c23d,
    float *restrict c24d, float *restrict c25d, float *restrict c26d,
                                                float *restrict c33d,
    float *restrict c34d, float *restrict c35d, float *restrict c36d,
    float *restrict c44d, float *restrict c45d, float *restrict c46d,
                          float *restrict c55d, float *restrict c56d,
                                                float *restrict c66d,
                                                 float *restrict slw3d,
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
  float c11,c12,c13,c14,c15,c16;
  float     c22,c23,c24,c25,c26;
  float         c33,c34,c35,c36;
  float             c44,c45,c46;
  float                 c55,c56;
  float                     c66;
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
        c12 = c12d[iptr];
        c13 = c13d[iptr];
        c14 = c14d[iptr];
        c15 = c15d[iptr];
        c16 = c16d[iptr];
        c22 = c22d[iptr];
        c23 = c23d[iptr];
        c24 = c24d[iptr];
        c25 = c25d[iptr];
        c26 = c26d[iptr];
        c33 = c33d[iptr];
        c34 = c34d[iptr];
        c35 = c35d[iptr];
        c36 = c36d[iptr];
        c44 = c44d[iptr];
        c45 = c45d[iptr];
        c46 = c46d[iptr];
        c55 = c55d[iptr];
        c56 = c56d[iptr];
        c66 = c66d[iptr];

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
	      hTxx[iptr] = (c11*xix + c16*xiy + c15*xiz) * DxVx + (c16*xix + c12*xiy + c14*xiz) * DxVy + (c15*xix + c14*xiy + c13*xiz) * DxVz
                   + (c11*etx + c16*ety + c15*etz) * DyVx + (c16*etx + c12*ety + c14*etz) * DyVy + (c15*etx + c14*ety + c13*etz) * DyVz
                   + (c11*ztx + c16*zty + c15*ztz) * DzVx + (c16*ztx + c12*zty + c14*ztz) * DzVy + (c15*ztx + c14*zty + c13*ztz) * DzVz;
      
        hTyy[iptr] = (c12*xix + c26*xiy + c25*xiz) * DxVx + (c26*xix + c22*xiy + c24*xiz) * DxVy + (c25*xix + c24*xiy + c23*xiz) * DxVz
                   + (c12*etx + c26*ety + c25*etz) * DyVx + (c26*etx + c22*ety + c24*etz) * DyVy + (c25*etx + c24*ety + c23*etz) * DyVz
                   + (c12*ztx + c26*zty + c25*ztz) * DzVx + (c26*ztx + c22*zty + c24*ztz) * DzVy + (c25*ztx + c24*zty + c23*ztz) * DzVz;
     
        hTzz[iptr] = (c13*xix + c36*xiy + c35*xiz) * DxVx + (c36*xix + c23*xiy + c34*xiz) * DxVy + (c35*xix + c34*xiy + c33*xiz) * DxVz
                   + (c13*etx + c36*ety + c35*etz) * DyVx + (c36*etx + c23*ety + c34*etz) * DyVy + (c35*etx + c34*ety + c33*etz) * DyVz
                   + (c13*ztx + c36*zty + c35*ztz) * DzVx + (c36*ztx + c23*zty + c34*ztz) * DzVy + (c35*ztx + c34*zty + c33*ztz) * DzVz;
  

        hTyz[iptr] = (c14*xix + c46*xiy + c45*xiz) * DxVx + (c46*xix + c24*xiy + c44*xiz) * DxVy + (c45*xix + c44*xiy + c34*xiz) * DxVz
                   + (c14*etx + c46*ety + c45*etz) * DyVx + (c46*etx + c24*ety + c44*etz) * DyVy + (c45*etx + c44*ety + c34*etz) * DyVz
                   + (c14*ztx + c46*zty + c45*ztz) * DzVx + (c46*ztx + c24*zty + c44*ztz) * DzVy + (c45*ztx + c44*zty + c34*ztz) * DzVz;
  
        hTxz[iptr] = (c15*xix + c56*xiy + c55*xiz) * DxVx + (c56*xix + c25*xiy + c45*xiz) * DxVy + (c55*xix + c45*xiy + c35*xiz) * DxVz
                   + (c15*etx + c56*ety + c55*etz) * DyVx + (c56*etx + c25*ety + c45*etz) * DyVy + (c55*etx + c45*ety + c35*etz) * DyVz
                   + (c15*ztx + c56*zty + c55*ztz) * DzVx + (c56*ztx + c25*zty + c45*ztz) * DzVy + (c55*ztx + c45*zty + c35*ztz) * DzVz;
  
        hTxy[iptr] = (c16*xix + c66*xiy + c56*xiz) * DxVx + (c66*xix + c26*xiy + c46*xiz) * DxVy + (c56*xix + c46*xiy + c36*xiz) * DxVz
                   + (c16*etx + c66*ety + c56*etz) * DyVx + (c66*etx + c26*ety + c46*etz) * DyVy + (c56*etx + c46*ety + c36*etz) * DyVz
                   + (c16*ztx + c66*zty + c56*ztz) * DzVx + (c66*ztx + c26*zty + c46*ztz) * DzVy + (c56*ztx + c46*zty + c36*ztz) * DzVz;

        iptr += 1;
      }
    }
  }
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

void
sv_curv_col_el_aniso_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict c11d, float *restrict c12d, float *restrict c13d,
    float *restrict c14d, float *restrict c15d, float *restrict c16d,
                          float *restrict c22d, float *restrict c23d,
    float *restrict c24d, float *restrict c25d, float *restrict c26d,
                                                float *restrict c33d,
    float *restrict c34d, float *restrict c35d, float *restrict c36d,
    float *restrict c44d, float *restrict c45d, float *restrict c46d,
                          float *restrict c55d, float *restrict c56d,
                                                float *restrict c66d,
                                                  float *restrict slw3d,
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
  float c11,c12,c13,c14,c15,c16;
  float     c22,c23,c24,c25,c26;
  float         c33,c34,c35,c36;
  float             c44,c45,c46;
  float                 c55,c56;
  float                     c66;
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
              c12 = c12d[iptr];
              c13 = c13d[iptr];
              c14 = c14d[iptr];
              c15 = c15d[iptr];
              c16 = c16d[iptr];
              c22 = c22d[iptr];
              c23 = c23d[iptr];
              c24 = c24d[iptr];
              c25 = c25d[iptr];
              c26 = c26d[iptr];
              c33 = c33d[iptr];
              c34 = c34d[iptr];
              c35 = c35d[iptr];
              c36 = c36d[iptr];
              c44 = c44d[iptr];
              c45 = c45d[iptr];
              c46 = c46d[iptr];
              c55 = c55d[iptr];
              c56 = c56d[iptr];
              c66 = c66d[iptr];

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
              hTxx_rhs = (c11*xix+c16*xiy+c15*xiz)*DxVx + (c16*xix+c12*xiy+c14*xiz)*DxVy + (c15*xix+c14*xiy+c13*xiz)*DxVz; 
              hTyy_rhs = (c12*xix+c26*xiy+c25*xiz)*DxVx + (c26*xix+c22*xiy+c24*xiz)*DxVy + (c25*xix+c24*xiy+c23*xiz)*DxVz;
              hTzz_rhs = (c13*xix+c36*xiy+c35*xiz)*DxVx + (c36*xix+c23*xiy+c34*xiz)*DxVy + (c35*xix+c34*xiy+c33*xiz)*DxVz;
              hTyz_rhs = (c14*xix+c46*xiy+c45*xiz)*DxVx + (c46*xix+c24*xiy+c44*xiz)*DxVy + (c45*xix+c44*xiy+c34*xiz)*DxVz;
              hTxz_rhs = (c15*xix+c56*xiy+c55*xiz)*DxVx + (c56*xix+c25*xiy+c45*xiz)*DxVy + (c55*xix+c45*xiy+c35*xiz)*DxVz;
              hTxy_rhs = (c16*xix+c66*xiy+c56*xiz)*DxVx + (c66*xix+c26*xiy+c46*xiz)*DxVy + (c56*xix+c46*xiy+c36*xiz)*DxVz;

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
                hTxx_rhs = (c11*ztx+c16*zty+c15*ztz)*Dx_DzVx + (c16*ztx+c12*zty+c14*ztz)*Dx_DzVy + (c15*ztx+c14*zty+c13*ztz)*Dx_DzVz; 
                hTyy_rhs = (c12*ztx+c26*zty+c25*ztz)*Dx_DzVx + (c26*ztx+c22*zty+c24*ztz)*Dx_DzVy + (c25*ztx+c24*zty+c23*ztz)*Dx_DzVz;
                hTzz_rhs = (c13*ztx+c36*zty+c35*ztz)*Dx_DzVx + (c36*ztx+c23*zty+c34*ztz)*Dx_DzVy + (c35*ztx+c34*zty+c33*ztz)*Dx_DzVz;
                hTyz_rhs = (c14*ztx+c46*zty+c45*ztz)*Dx_DzVx + (c46*ztx+c24*zty+c44*ztz)*Dx_DzVy + (c45*ztx+c44*zty+c34*ztz)*Dx_DzVz;
                hTxz_rhs = (c15*ztx+c56*zty+c55*ztz)*Dx_DzVx + (c56*ztx+c25*zty+c45*ztz)*Dx_DzVy + (c55*ztx+c45*zty+c35*ztz)*Dx_DzVz;
                hTxy_rhs = (c16*ztx+c66*zty+c56*ztz)*Dx_DzVx + (c66*ztx+c26*zty+c46*ztz)*Dx_DzVy + (c56*ztx+c46*zty+c36*ztz)*Dx_DzVz;

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
              c12 = c12d[iptr];
              c13 = c13d[iptr];
              c14 = c14d[iptr];
              c15 = c15d[iptr];
              c16 = c16d[iptr];
              c22 = c22d[iptr];
              c23 = c23d[iptr];
              c24 = c24d[iptr];
              c25 = c25d[iptr];
              c26 = c26d[iptr];
              c33 = c33d[iptr];
              c34 = c34d[iptr];
              c35 = c35d[iptr];
              c36 = c36d[iptr];
              c44 = c44d[iptr];
              c45 = c45d[iptr];
              c46 = c46d[iptr];
              c55 = c55d[iptr];
              c56 = c56d[iptr];
              c66 = c66d[iptr];

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
              hTxx_rhs = (c11*etx+c16*ety+c15*etz)*DyVx + (c16*etx+c12*ety+c14*etz)*DyVy + (c15*etx+c14*ety+c13*etz)*DyVz; 
              hTyy_rhs = (c12*etx+c26*ety+c25*etz)*DyVx + (c26*etx+c22*ety+c24*etz)*DyVy + (c25*etx+c24*ety+c23*etz)*DyVz;
              hTzz_rhs = (c13*etx+c36*ety+c35*etz)*DyVx + (c36*etx+c23*ety+c34*etz)*DyVy + (c35*etx+c34*ety+c33*etz)*DyVz;
              hTyz_rhs = (c14*etx+c46*ety+c45*etz)*DyVx + (c46*etx+c24*ety+c44*etz)*DyVy + (c45*etx+c44*ety+c34*etz)*DyVz;
              hTxz_rhs = (c15*etx+c56*ety+c55*etz)*DyVx + (c56*etx+c25*ety+c45*etz)*DyVy + (c55*etx+c45*ety+c35*etz)*DyVz;
              hTxy_rhs = (c16*etx+c66*ety+c56*etz)*DyVx + (c66*etx+c26*ety+c46*etz)*DyVy + (c56*etx+c46*ety+c36*etz)*DyVz;

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
                hTxx_rhs = (c11*ztx+c16*zty+c15*ztz)*Dy_DzVx + (c16*ztx+c12*zty+c14*ztz)*Dy_DzVy + (c15*ztx+c14*zty+c13*ztz)*Dy_DzVz; 
                hTyy_rhs = (c12*ztx+c26*zty+c25*ztz)*Dy_DzVx + (c26*ztx+c22*zty+c24*ztz)*Dy_DzVy + (c25*ztx+c24*zty+c23*ztz)*Dy_DzVz;
                hTzz_rhs = (c13*ztx+c36*zty+c35*ztz)*Dy_DzVx + (c36*ztx+c23*zty+c34*ztz)*Dy_DzVy + (c35*ztx+c34*zty+c33*ztz)*Dy_DzVz;
                hTyz_rhs = (c14*ztx+c46*zty+c45*ztz)*Dy_DzVx + (c46*ztx+c24*zty+c44*ztz)*Dy_DzVy + (c45*ztx+c44*zty+c34*ztz)*Dy_DzVz;
                hTxz_rhs = (c15*ztx+c56*zty+c55*ztz)*Dy_DzVx + (c56*ztx+c25*zty+c45*ztz)*Dy_DzVy + (c55*ztx+c45*zty+c35*ztz)*Dy_DzVz;
                hTxy_rhs = (c16*ztx+c66*zty+c56*ztz)*Dy_DzVx + (c66*ztx+c26*zty+c46*ztz)*Dy_DzVy + (c56*ztx+c46*zty+c36*ztz)*Dy_DzVz;

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
              c12 = c12d[iptr];
              c13 = c13d[iptr];
              c14 = c14d[iptr];
              c15 = c15d[iptr];
              c16 = c16d[iptr];
              c22 = c22d[iptr];
              c23 = c23d[iptr];
              c24 = c24d[iptr];
              c25 = c25d[iptr];
              c26 = c26d[iptr];
              c33 = c33d[iptr];
              c34 = c34d[iptr];
              c35 = c35d[iptr];
              c36 = c36d[iptr];
              c44 = c44d[iptr];
              c45 = c45d[iptr];
              c46 = c46d[iptr];
              c55 = c55d[iptr];
              c56 = c56d[iptr];
              c66 = c66d[iptr];

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
              hTxx_rhs = (c11*ztx+c16*zty+c15*ztz)*DzVx + (c16*ztx+c12*zty+c14*ztz)*DzVy + (c15*ztx+c14*zty+c13*ztz)*DzVz; 
              hTyy_rhs = (c12*ztx+c26*zty+c25*ztz)*DzVx + (c26*ztx+c22*zty+c24*ztz)*DzVy + (c25*ztx+c24*zty+c23*ztz)*DzVz;
              hTzz_rhs = (c13*ztx+c36*zty+c35*ztz)*DzVx + (c36*ztx+c23*zty+c34*ztz)*DzVy + (c35*ztx+c34*zty+c33*ztz)*DzVz;
              hTyz_rhs = (c14*ztx+c46*zty+c45*ztz)*DzVx + (c46*ztx+c24*zty+c44*ztz)*DzVy + (c45*ztx+c44*zty+c34*ztz)*DzVz;
              hTxz_rhs = (c15*ztx+c56*zty+c55*ztz)*DzVx + (c56*ztx+c25*zty+c45*ztz)*DzVy + (c55*ztx+c45*zty+c35*ztz)*DzVz;
              hTxy_rhs = (c16*ztx+c66*zty+c56*ztz)*DzVx + (c66*ztx+c26*zty+c46*ztz)*DzVy + (c56*ztx+c46*zty+c36*ztz)*DzVz;

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
sv_curv_col_el_aniso_dvh2dvz(gdinfo_t        *gdinfo,
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
  float *restrict c12d = md->c12;
  float *restrict c13d = md->c13;
  float *restrict c14d = md->c14;
  float *restrict c15d = md->c15;
  float *restrict c16d = md->c16;
  float *restrict c22d = md->c22;
  float *restrict c23d = md->c23;
  float *restrict c24d = md->c24;
  float *restrict c25d = md->c25;
  float *restrict c26d = md->c26;
  float *restrict c33d = md->c33;
  float *restrict c34d = md->c34;
  float *restrict c35d = md->c35;
  float *restrict c36d = md->c36;
  float *restrict c44d = md->c44;
  float *restrict c45d = md->c45;
  float *restrict c46d = md->c46;
  float *restrict c55d = md->c55;
  float *restrict c56d = md->c56;
  float *restrict c66d = md->c66;

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;

  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];

  float c11,c12,c13,c14,c15,c16;
  float     c22,c23,c24,c25,c26;
  float         c33,c34,c35,c36;
  float             c44,c45,c46;
  float                 c55,c56;
  float                     c66;
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
      c12 = c12d[iptr];
      c13 = c13d[iptr];
      c14 = c14d[iptr];
      c15 = c15d[iptr];
      c16 = c16d[iptr];
      c22 = c22d[iptr];
      c23 = c23d[iptr];
      c24 = c24d[iptr];
      c25 = c25d[iptr];
      c26 = c26d[iptr];
      c33 = c33d[iptr];
      c34 = c34d[iptr];
      c35 = c35d[iptr];
      c36 = c36d[iptr];
      c44 = c44d[iptr];
      c45 = c45d[iptr];
      c46 = c46d[iptr];
      c55 = c55d[iptr];
      c56 = c56d[iptr];
      c66 = c66d[iptr];

      // first dim: irow; sec dim: jcol, as Fortran code
      A[0][0] = (c11*ztx+c16*zty+c15*ztz)*ztx + (c16*ztx+c66*zty+c56*ztz)*zty + (c15*ztx+c56*zty+c55*ztz)*ztz;
      A[0][1] = (c16*ztx+c12*zty+c14*ztz)*ztx + (c66*ztx+c26*zty+c46*ztz)*zty + (c56*ztx+c25*zty+c45*ztz)*ztz;
      A[0][2] = (c15*ztx+c14*zty+c13*ztz)*ztx + (c56*ztx+c46*zty+c36*ztz)*zty + (c55*ztx+c45*zty+c35*ztz)*ztz; 
      A[1][0] = (c16*ztx+c66*zty+c56*ztz)*ztx + (c12*ztx+c26*zty+c25*ztz)*zty + (c14*ztx+c46*zty+c45*ztz)*ztz; 
      A[1][1] = (c66*ztx+c26*zty+c46*ztz)*ztx + (c26*ztx+c22*zty+c24*ztz)*zty + (c46*ztx+c24*zty+c44*ztz)*ztz; 
      A[1][2] = (c56*ztx+c46*zty+c36*ztz)*ztx + (c25*ztx+c24*zty+c23*ztz)*zty + (c45*ztx+c44*zty+c34*ztz)*ztz;
      A[2][0] = (c15*ztx+c56*zty+c55*ztz)*ztx + (c14*ztx+c46*zty+c45*ztz)*zty + (c13*ztx+c36*zty+c35*ztz)*ztz;
      A[2][1] = (c56*ztx+c25*zty+c45*ztz)*ztx + (c46*ztx+c24*zty+c44*ztz)*zty + (c36*ztx+c23*zty+c34*ztz)*ztz;
      A[2][2] = (c55*ztx+c45*zty+c35*ztz)*ztx + (c45*ztx+c44*zty+c34*ztz)*zty + (c35*ztx+c34*zty+c33*ztz)*ztz; 
      fdlib_math_invert3x3(A);
                                                       
      B[0][0] = (c11*xix+c16*xiy+c15*xiz)*ztx + (c16*xix+c66*xiy+c56*xiz)*zty + (c15*xix+c56*xiy+c55*xiz)*ztz;
      B[0][1] = (c16*xix+c12*xiy+c14*xiz)*ztx + (c66*xix+c26*xiy+c46*xiz)*zty + (c56*xix+c25*xiy+c45*xiz)*ztz;
      B[0][2] = (c15*xix+c14*xiy+c13*xiz)*ztx + (c56*xix+c46*xiy+c36*xiz)*zty + (c55*xix+c45*xiy+c35*xiz)*ztz; 
      B[1][0] = (c16*xix+c66*xiy+c56*xiz)*ztx + (c12*xix+c26*xiy+c25*xiz)*zty + (c14*xix+c46*xiy+c45*xiz)*ztz; 
      B[1][1] = (c66*xix+c26*xiy+c46*xiz)*ztx + (c26*xix+c22*xiy+c24*xiz)*zty + (c46*xix+c24*xiy+c44*xiz)*ztz; 
      B[1][2] = (c56*xix+c46*xiy+c36*xiz)*ztx + (c25*xix+c24*xiy+c23*xiz)*zty + (c45*xix+c44*xiy+c34*xiz)*ztz;
      B[2][0] = (c15*xix+c56*xiy+c55*xiz)*ztx + (c14*xix+c46*xiy+c45*xiz)*zty + (c13*xix+c36*xiy+c35*xiz)*ztz;
      B[2][1] = (c56*xix+c25*xiy+c45*xiz)*ztx + (c46*xix+c24*xiy+c44*xiz)*zty + (c36*xix+c23*xiy+c34*xiz)*ztz;
      B[2][2] = (c55*xix+c45*xiy+c35*xiz)*ztx + (c45*xix+c44*xiy+c34*xiz)*zty + (c35*xix+c34*xiy+c33*xiz)*ztz; 
       
      C[0][0] = (c11*etx+c16*ety+c15*etz)*ztx + (c16*etx+c66*ety+c56*etz)*zty + (c15*etx+c56*ety+c55*etz)*ztz;
      C[0][1] = (c16*etx+c12*ety+c14*etz)*ztx + (c66*etx+c26*ety+c46*etz)*zty + (c56*etx+c25*ety+c45*etz)*ztz;
      C[0][2] = (c15*etx+c14*ety+c13*etz)*ztx + (c56*etx+c46*ety+c36*etz)*zty + (c55*etx+c45*ety+c35*etz)*ztz; 
      C[1][0] = (c16*etx+c66*ety+c56*etz)*ztx + (c12*etx+c26*ety+c25*etz)*zty + (c14*etx+c46*ety+c45*etz)*ztz; 
      C[1][1] = (c66*etx+c26*ety+c46*etz)*ztx + (c26*etx+c22*ety+c24*etz)*zty + (c46*etx+c24*ety+c44*etz)*ztz; 
      C[1][2] = (c56*etx+c46*ety+c36*etz)*ztx + (c25*etx+c24*ety+c23*etz)*zty + (c45*etx+c44*ety+c34*etz)*ztz;
      C[2][0] = (c15*etx+c56*ety+c55*etz)*ztx + (c14*etx+c46*ety+c45*etz)*zty + (c13*etx+c36*ety+c35*etz)*ztz;
      C[2][1] = (c56*etx+c25*ety+c45*etz)*ztx + (c46*etx+c24*ety+c44*etz)*zty + (c36*etx+c23*ety+c34*etz)*ztz;
      C[2][2] = (c55*etx+c45*ety+c35*etz)*ztx + (c45*etx+c44*ety+c34*etz)*zty + (c35*etx+c34*ety+c33*etz)*ztz; 
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

/*******************************************************************************
 * ABL-EXP boundary
 ******************************************************************************/

/* 
// apply on wavefield, not on derivatives
int sv_eliso1st_curv_macdrp_apply_ablexp_pointcoef(float *restrict w_cur, 
    int number_of_vars, int *restrict vars_pos,
    // grid size
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    int ni, int nj, int nk, int nx, int ny, int nz,
    // ablexp info
    int *restrict abs_indx, float *restrict abs_damp
    )
{
  int ierr = 0;

  int iptr,iptr_abs, iptr_var, iptr_k, iptr_j;
  size_t siz_line, siz_slice;
  int i,j,k,ivar;

  siz_line  = nx;
  siz_slice = nx * ny;

  for (ivar=0; ivar<number_of_vars; ivar++)
  {
    iptr_var = vars_pos[ivar];
    iptr_abs = 0;

    for (k=abs_indx[4]; k<=abs_indx[5]; k++)
    {
      iptr_k = iptr_var + k * siz_slice;
      for (j=abs_indx[2]; j<=abs_indx[3]; j++)
      {
        iptr_j = j * siz_line + iptr_k;
        for (i=abs_indx[0]; i<=abs_indx[2]; i++)
        {
          iptr = i + iptr_j;

          w_cur[iptr] *= abs_damp[iptr_abs];

          // next
          iptr_abs += 1;
        }
      }
    }
  }

  return ierr;
}

int sv_eliso1st_curv_macdrp_apply_ablexp(float *restrict w_cur, 
    int number_of_vars, int *restrict vars_pos,
    // grid size
    int nx, int ny, int nz,
    // ablexp info
    int *restrict abs_indx, float *restrict abs_damp
    )
{
  int ierr = 0;

  int iptr,iptr_abs, iptr_var, iptr_k, iptr_j;
  size_t siz_line, siz_slice;
  int i,j,k,ivar;
  float *restrict Dx, Dy, Dz;
  float dampx, dampy, dampz, dampmin;

  siz_line  = nx;
  siz_slice = nx * ny;

  Dx = abs_damp;
  Dy = Dx + nx;
  Dz = Dy + ny; 

  for (ivar=0; ivar<number_of_vars; ivar++)
  {
    iptr_var = vars_pos[ivar];

    for (k=abs_indx[4]; k<=abs_indx[5]; k++)
    {
      iptr_k = iptr_var + k * siz_slice;
      dampz = Dz[k];

      for (j=abs_indx[2]; j<=abs_indx[3]; j++)
      {
        iptr_j = j * siz_line + iptr_k;
        dampy = Dy[j];
        dampmin = (dampy < dampz) ? dampy : dampz;

        for (i=abs_indx[0]; i<=abs_indx[2]; i++)
        {
          iptr = i + iptr_j;
          dampx = Dx[i];
          dampmin = (dampx < dampy) ? dampx : dampy;

          w_cur[iptr] *= dampmin;

        }
      }
    }
  }

  return ierr;
}
*/
