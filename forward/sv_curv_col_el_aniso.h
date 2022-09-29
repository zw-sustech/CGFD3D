#ifndef SV_CURV_COL_EL_ANISO_H
#define SV_CURV_COL_EL_ANISO_H

#include "fd_t.h"
#include "mympi_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

void
sv_curv_col_el_aniso_onestage(
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
  const int myid, const int verbose);

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
    const int myid, const int verbose);

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
    const int myid, const int verbose);

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
    const int myid, const int verbose);

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
    const int myid, const int verbose);

int
sv_curv_col_el_aniso_dvh2dvz(gd_t        *gdinfo,
                                   gdcurv_metric_t *metric,
                                   md_t       *md,
                                   bdry_t      *bdryfree,
                                   const int verbose);

#endif
