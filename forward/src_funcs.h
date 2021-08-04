#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

#include <mpi.h>

#include "constants.h"
#include "gd_info.h"
#include "gd_curv.h"

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
src_gen_single_point_gauss(gdinfo_t *gdinfo,
                           gdcurv_t *gdcurv,
                           src_t    *src,
                           float t0,
                           float dt,
                           int   num_of_stages,
                           float *rk_stage_time,
                           int   npoint_half_ext,
                           char  *source_name,
                           int   *source_gridindex,
                           float *source_coords,
                           float *force_vector,
                           float *moment_tensor,
                           char  *wavelet_name,
                           float *wavelet_coefs,
                           float wavelet_tstart,
                           float wavelet_tend,
                           MPI_Comm comm, 
                           int myid,
                           int verbose);

int
src_read_locate_valsrc(char *pfilepath,
                      size_t siz_line,
                      size_t siz_slice,
                      float t0,
                      float dt,
                      int   num_of_stages,
                      float *rk_stage_time,
                      int   glob_phys_ix1, // gloabl start index along x this thread
                      int   glob_phys_ix2, // gloabl end index along x
                      int   glob_phys_iy1,
                      int   glob_phys_iy2,
                      int   glob_phys_iz1,
                      int   glob_phys_iz2,
                      int   ni1,
                      int   ni2,
                      int   nj1,
                      int   nj2,
                      int   nk1,
                      int   nk2,
                      int   npoint_half_ext,
                      int   npoint_ghosts,
                      float *x3d,
                      float *y3d,
                      float *z3d,
                      MPI_Comm comm,
                      int myid,
                      // following output
                      char **p_event_name,
                      int  *num_of_force, // inout: if force source, if in this thread
                      int **restrict p_force_info,
                      float  **restrict p_force_vec_stf,
                      int    **restrict p_force_ext_indx,
                      float  **restrict p_force_ext_coef,
                      int  *num_of_moment, // inout: if moment source, if in this thread
                      int    **restrict p_moment_info,
                      float  **restrict p_moment_ten_rate,
                      int    **restrict p_moment_ext_indx,
                      float  **restrict p_moment_ext_coef,
                      int verbose);

int
src_read_locate_anasrc(char *pfilepath,
                      size_t siz_line,
                      size_t siz_slice,
                      float t0,
                      float dt,
                      int   num_of_stages,
                      float *rk_stage_time,
                      int   glob_phys_ix1, // gloabl start index along x this thread
                      int   glob_phys_ix2, // gloabl end index along x
                      int   glob_phys_iy1,
                      int   glob_phys_iy2,
                      int   glob_phys_iz1,
                      int   glob_phys_iz2,
                      int   ni1,
                      int   ni2,
                      int   nj1,
                      int   nj2,
                      int   nk1,
                      int   nk2,
                      int   npoint_half_ext,
                      int   npoint_ghosts,
                      float *x3d,
                      float *y3d,
                      float *z3d,
                      MPI_Comm comm,
                      int myid,
                      // following output
                      char **p_event_name,
                      int   *num_of_force, // force in this thread
                      int **restrict p_force_info,
                      float  **restrict p_force_vec_stf,
                      int    **restrict p_force_ext_indx,
                      float  **restrict p_force_ext_coef,
                      int   *num_of_moment, // moment in this thread
                      int    **restrict p_moment_info,
                      float  **restrict p_moment_ten_rate,
                      int    **restrict p_moment_ext_indx,
                      float  **restrict p_moment_ext_coef,
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
cal_norm_delt3d(float *delt, float x0, float y0, float z0, float rx0, float ry0, float rz0, int LenDelt);

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

#endif
