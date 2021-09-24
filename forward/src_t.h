#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

#include <mpi.h>

#include "constants.h"
#include "gd_info.h"
#include "gd_t.h"

// cal force_vec_stf/moment_ten_rate 1d index for icmp,it,istage
//  with respect to the start pointer of this source point
#define M_SRC_IND(icmp,it,istage,nt,num_stage) \
  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int total_number;
  int max_nt; // max nt of stf and mrf per src
  int max_stage; // max number of rk stages
  int max_ext; // max extened points

  // for getting value in calculation
  int it;
  int istage;

  // for output
  char evtnm[CONST_MAX_STRLEN];

  // time independent
  int *si; // local i index 
  int *sj; // local j index 
  int *sk; // local k index 
  int *it_begin; // start t index
  int *it_end;   // end   t index
  int   *ext_num; // valid extend points for this src
  int   *ext_indx; // max_ext * total_number
  float *ext_coef;

  // force and/or moment
  int force_actived;
  int moment_actived;

  // time dependent
  // force stf
  float *Fx; // max_stage * max_nt * total_number;
  float *Fy;
  float *Fz;
  // moment rate
  float *Mxx; // max_stage *max_nt * total_number;
  float *Myy;
  float *Mzz;
  float *Mxz;
  float *Myz;
  float *Mxy;
} src_t;

/*************************************************
 * function prototype
 *************************************************/

int
src_coord_to_glob_indx(gdinfo_t *gdinfo,
                       gd_t *gdcurv,
                       float sx,
                       float sy,
                       float sz,
                       MPI_Comm comm,
                       int myid,
                       int   *ou_si, int *ou_sj, int *ou_sk,
                       float *ou_sx_inc, float *ou_sy_inc, float *ou_sz_inc,
                       float *restrict wrk3d);

int
src_glob_ext_ishere(int si, int sj, int sk, int half_ext, gdinfo_t *gdinfo);

int
src_coord_to_local_indx(gdinfo_t *gdinfo,
                        gd_t *gdcurv,
                        float sx, float sy, float sz,
                        int *si, int *sj, int *sk,
                        float *sx_inc, float *sy_inc, float *sz_inc,
                        float *restrict wrk3d);

int
src_set_by_par(gdinfo_t *gdinfo,
               gd_t *gdcurv,
               src_t    *src,
               float t0,
               float dt,
               int   max_stage,
               float *rk_stage_time,
               int   npoint_half_ext,
               char  *in_source_name,
               int   in_num_of_src,
               int   **source_index,
               float **source_inc,
               float **source_coords,
               float **force_vector, 
               int   *source_force_actived,
               float **moment_tensor,
               int   *source_moment_actived,
               char  **wavelet_name,
               float **wavelet_coefs,
               float *wavelet_tstart,
               float *wavelet_tend,
               MPI_Comm comm, 
               int myid,
               int verbose);

int
src_read_locate_valsrc(gdinfo_t *gdinfo,
                       gd_t *gdcurv,
                       src_t    *src,
                       char *pfilepath,
                       float t0,
                       float dt,
                       int   max_stages,
                       float *rk_stage_time,
                       int   npoint_half_ext,
                       MPI_Comm comm,
                       int myid,
                       int verbose);

int
src_read_locate_anasrc(gdinfo_t *gdinfo,
                       gd_t *gdcurv,
                       src_t    *src,
                       char *pfilepath,
                       float t0,
                       float dt,
                       int   max_stages,
                       float *rk_stage_time,
                       int   npoint_half_ext,
                       MPI_Comm comm,
                       int myid,
                       int verbose);

float 
fun_ricker(float t, float fc, float t0);

float 
fun_ricker_deriv(float t, float fc, float t0);

float
fun_gauss(float t, float a, float t0);

float
fun_gauss_deriv(float t, float a, float t0);

void 
angle2moment(float strike, float dip, float rake, float* source_moment_tensor);

int
src_coord2index(float sx, float sy, float sz,
                int nx, int ny, int nz,
                int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
                float *restrict x3d,
                float *restrict y3d,
                float *restrict z3d,
                float *restrict wrk3d,
                int *si, int *sj, int *sk,
                float *sx_inc, float *sy_inc, float *sz_inc);

int
src_cart2curv_rdinterp(float sx, float sy, float sz, 
                       int num_points,
                       float *points_x, // x coord of all points
                       float *points_y,
                       float *points_z,
                       float *points_i, // curv coord of all points
                       float *points_j,
                       float *points_k,
                       float *si_curv, // interped curv coord
                       float *sj_curv,
                       float *sk_curv);

int
src_cart2curv_sample(float sx, float sy, float sz, 
                     int num_points,
                     float *points_x, // x coord of all points
                     float *points_y,
                     float *points_z,
                     float *points_i, // curv coord of all points
                     float *points_j,
                     float *points_k,
                     int    nx_sample,
                     int    ny_sample,
                     int    nz_sample,
                     float *si_curv, // interped curv coord
                     float *sj_curv,
                     float *sk_curv);

int
src_set_time(src_t *src, int it, int istage);

void
src_cal_norm_delt3d(float *delt, float x0, float y0, float z0,
                    float rx0, float ry0, float rz0, int LenDelt);

#endif
