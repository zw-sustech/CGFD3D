#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

#include <mpi.h>

#ifndef PI
#define PI 3.1415926535898
#endif

#define M_SRC_INFO_SEQ_SI 0
#define M_SRC_INFO_SEQ_SJ 1
#define M_SRC_INFO_SEQ_SK 2
// position in stf or mvf vector for this source
#define M_SRC_INFO_SEQ_POS 3
#define M_SRC_INFO_SEQ_ITBEG 4
#define M_SRC_INFO_SEQ_ITEND 5
#define M_SRC_INFO_SEQ_NEXT_MAX 6
#define M_SRC_INFO_SEQ_NEXT_THIS 7

#define M_SRC_INFO_NVAL 8

// cal force_vec_stf/moment_ten_rate 1d index for icmp,it,istage
//  with respect to the start pointer of this source point
#define M_SRC_IND(icmp,it,istage,nt,num_stage) \
  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

void
src_gen_single_point_gauss(size_t siz_line,
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
                           // following output
                           struct fd_src_t *src,
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
src_get_stage_stf(
                  int num_of_force,
                  int  *restrict force_info, // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
                  float *restrict force_vec_stf,
                  int num_of_moment,
                  int  *restrict moment_info, // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
                  float *restrict moment_ten_rate,
                  int it, int istage, int num_of_stages,
                  float *restrict force_vec_value,
                  float *restrict moment_ten_value,
                  const int myid, const int verbose);

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

#endif
