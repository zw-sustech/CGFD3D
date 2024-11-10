#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

#include <mpi.h>

#include "constants.h"
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

  /*
   * add vars for distributed and discrete (dd) sources
   *  reload stf during time marching
   *  all points have same start and end t to simplify io and apply
   */

  // if dd should be added for this it step
  int dd_is_valid;

  int dd_is_add_at_point;
  int dd_smo_hlen;

  // force and/or moment
  int dd_vi_actived;
  int dd_mij_actived;

  int dd_total_number;
  int dd_max_nt;
  int dd_nt_per_read;  // time block size

  int dd_nt_this_read;  // time block size of this read
  int dd_it_here; // cur it relative to it0 of this read 

  // time independent
  int *dd_si; // local i index 
  int *dd_sj; // local j index 
  int *dd_sk; // local k index 
  size_t *dd_indx; // 1d index

  float *dd_si_inc; // local i shift 
  float *dd_sj_inc; // local j shift 
  float *dd_sk_inc; // local k shift 

  // file id
  FILE *fp_vi;
  FILE *fp_mij;

  // vars
  float *dd_vi; // 3 * max_stage * nt_per_read * total_number;
  float *dd_mij;

} src_t;

/*************************************************
 * function prototype
 *************************************************/

int
src_coord_to_glob_indx(
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
src_glob_ext_ishere(int si, int sj, int sk, int half_ext, gd_t *gdinfo);

int
src_glob_ishere(int si, int sj, int sk, int half_ext, gd_t *gdinfo);

int
src_coord_to_local_indx(
                        gd_t *gdcurv,
                        float sx, float sy, float sz,
                        int *si, int *sj, int *sk,
                        float *sx_inc, float *sy_inc, float *sz_inc,
                        float *restrict wrk3d);

int
src_read_locate_file(
                     gd_t     *gd,
                     src_t    *src,
                     float    *restrict mu3d,
                     char     *in_src_file,
                     float     t0,
                     float     dt,
                     int       max_stage,
                     float    *rk_stage_time,
                     int       npoint_half_ext,
                     MPI_Comm  comm,
                     int       myid,
                     int       verbose);

int
src_dd_read2local(
                     gd_t     *gd,
                     src_t    *src,
                     char     *in_file,
                     char     *tmp_dir,
                     int       dd_is_add_at_point,
                     int       dd_nt_per_read,
                     float     t0,
                     float     dt,
                     int       nt_total,
                     int       max_stage,
                     float    *rk_stage_time,
                     int       npoint_half_ext,
                     MPI_Comm  comm,
                     int       myid,
                     int*      topoid,
                     int       verbose);


int
src_read_meta_src(FILE *fp, char *evtnm, int *ns, 
                 int *is_stf_given, float *stf_length, float *stf_dt, int *stf_nt,
                 int *cmp_type, int *mechanism_type, int *coord_type, int *z_type);

int
src_read_loc_all_src(FILE *fp, int num_source, float *sx_all, float *sy_all, float *sz_all);

int
src_skip_one_stfmech_src(FILE *fp, int is_stf_by_value, int in_stf_nt);

int
src_read_one_stf_name(FILE *fp, float *tstart, char *wavelet_name, float *wavelet_coefs);

int
src_read_one_mech_name(FILE *fp, gd_t *gd, float *mu3d, 
                        int in_cmp_type, int moment_actived, int in_mechanism_type,
                        int *evt_index, int in_this_thread,
                        float *evt_fx, float *evt_fy, float *evt_fz,
                        float *evt_mxx, float *evt_myy, float *evt_mzz,
                        float *evt_myz, float *evt_mxz, float *evt_mxy, float *evt_M0);

int
src_read_one_mech_value(FILE *fp, gd_t *gd, float *mu3d, int in_stf_nt, float in_stf_dt,
                        int in_cmp_type, int moment_actived, int in_mechanism_type,
                        int *evt_index,  int in_this_thread,
                        float *tstart, float *f1, float *f2, float *f3,
                        float *m11, float *m22, float *m33, float *m23, float *m13, float *m12, 
                        float *evt_M0);

int
src_put_into_struct(src_t *src, gd_t *gd,
                    float     t0,
                    float     dt,
                    int max_stage,
                    float *rk_stage_time,
                    int npoint_half_ext, 
                    int is_local, 
                    int in_stf_by_value,
                    int in_stf_nt, float in_stf_dt, float *t_in, int max_nt,
                    float wavelet_tstart, char *wavelet_name, float *wavelet_coefs,
                    int force_actived, int moment_actived,
                    int *evt_index, float *evt_inc, 
                    float fx, float fy, float fz,
                    float mxx, float myy, float mzz, float myz, float mxz, float mxy,
                    float *f1, float *f2, float *f3, 
                    float *m11, float *m22, float *m33, float *m23, float *m13, float *m12);

int
src_dd_free(src_t *src);

int
src_dd_accit_loadstf(src_t *src, int it, int myid);

float
src_cal_wavelet(float t, float dt, char *wavelet_name, float *wavelet_coefs);

float 
fun_ricker(float t, float fc, float t0);

float 
fun_ricker_deriv(float t, float fc, float t0);

float
fun_gauss(float t, float a, float t0);

float
fun_gauss_deriv(float t, float a, float t0);

float
fun_liuetal2006(float t, float tau);

float
fun_bell(float t, float riset);

float
fun_bell_deriv(float t, float riset);

float
fun_step(float t);

float
fun_delta(float t, float dt);

float
fun_klauder(float t, float tshift, float f1, float f2, float T);

float
fun_klauder_blackman(float t, float tshift, float f1, float f2, float T, float dt);

float
blackman_window(float t, float dt, float tshift);

void 
angle2moment(float strike, float dip, float rake, float* source_moment_tensor);

int
src_muDA_to_moment(float strike, float dip, float rake, float mu, float D, float A,
          float *mxx, float *myy, float *mzz, float *myz, float *mxz, float *mxy);

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

int
src_print(src_t *src, int verbose);

#endif
