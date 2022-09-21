#ifndef SV_CURV_COL_ISO_H
#define SV_CURV_COL_ISO_H

#include "gd_info.h"
#include "gd_t.h"
#include "md_t.h"
#include "src_t.h"

/*************************************************
 * function prototype
 *************************************************/

void
sv_curv_col_el_rhs_timg_z2(
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

int
sv_curv_col_el_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    src_t *src, // short nation for reference member
    const int myid, const int verbose);

int
sv_curv_col_el_rhs_srcdd(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    src_t *src, // short nation for reference member
    const int myid, const int verbose);

int
sv_curv_col_el_graves_Qs(float *w, int ncmp, float dt, gdinfo_t *gdinfo, md_t *md);

#endif
