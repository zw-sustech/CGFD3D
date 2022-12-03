#ifndef SV_CURV_COL_VIS_ISO_H
#define SV_CURV_COL_VIS_ISO_H

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
sv_curv_col_vis_iso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gd_t   *gdinfo,
  gdcurv_metric_t  *metric,
  md_t *mdeliso,
  bdry_t  *bdry,
  src_t *src,
  // include different order/stentil
  int num_of_fdx_op, fd_op_t *fdx_op,
  int num_of_fdy_op, fd_op_t *fdy_op,
  int num_of_fdz_op, fd_op_t *fdz_op,
  int fdz_max_len, 
  const int myid, const int verbose);

void
sv_curv_col_vis_iso_rhs_vlow_z2(
    float *restrict  Vx ,  float *restrict  Vy ,  float *restrict  Vz ,
    float *restrict hTxx,  float *restrict hTyy,  float *restrict hTzz,
    float *restrict hTxz,  float *restrict hTyz,  float *restrict hTxy,
    float **restrict Jxx,  float **restrict Jyy,  float **restrict Jzz,
    float **restrict Jxy,  float **restrict Jxz,  float **restrict Jyz,
    float **restrict hJxx, float **restrict hJyy, float **restrict hJzz,
    float **restrict hJxy, float **restrict hJxz, float **restrict hJyz,
    float *restrict xi_x,  float *restrict xi_y,  float *restrict xi_z,
    float *restrict et_x,  float *restrict et_y,  float *restrict et_z,
    float *restrict zt_x,  float *restrict zt_y,  float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d,  float *restrict slw3d,
    float *restrict wl, float **restrict Ylam, float **restrict Ymu,
    float *restrict matVx2Vz, float *restrict matVy2Vz, float *restrict matA,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
    int nmaxwell,
    const int myid, const int verbose);

void
sv_curv_col_vis_iso_atten(
    float *restrict hTxx,  float *restrict hTyy,  float *restrict hTzz,
    float *restrict hTxz,  float *restrict hTyz,  float *restrict hTxy,
    float **restrict Jxx,  float **restrict Jyy,  float **restrict Jzz,
    float **restrict Jxy,  float **restrict Jxz,  float **restrict Jyz,
    float **restrict hJxx, float **restrict hJyy, float **restrict hJzz,
    float **restrict hJxy, float **restrict hJxz, float **restrict hJyz,
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    float *restrict wl, float **restrict Ylam, float **restrict Ymu,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int nmaxwell,
    const int myid, const int verbose);

int
sv_curv_col_vis_iso_dvh2dvz(gd_t        *gdinfo,
                                 gdcurv_metric_t *metric,
                                 md_t       *md,
                                 bdry_t      *bdryfree,
                                 const int verbose);
#endif
