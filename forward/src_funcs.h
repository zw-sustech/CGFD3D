#ifndef SRC_FUNCS_H
#define SRC_FUNCS_H

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

struct CubicPt 
{
float coordx; 
float coordy; 
float coordz;

int xindx; 
int yindx;
int zindx;
} ;

struct SrcIndx
 {
  int si; 
  int sj; 
  int sk;
  float sx_inc;
  float sy_inc;
  float sz_inc;
 };

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
                           int   *source_gridindex,
                           float *source_coords,
                           float *force_vector,
                           float *moment_tensor,
                           char  *wavelet_name,
                           float *wavelet_coefs,
                           float wavelet_tstart,
                           float wavelet_tend,
                           // following output
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

void
src_gen_multiple_point_gauss(size_t siz_line,
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
                             char  *pfilepath,
                             // following output
                             int *num_of_force, // force in this thread
                             int **restrict p_force_info,
                             float  **restrict p_force_vec_stf,
                             int    **restrict p_force_ext_indx,
                             float  **restrict p_force_ext_coef,
                             int *num_of_moment, // moment in this thread
                             int    **restrict p_moment_info,
                             float  **restrict p_moment_ten_rate,
                             int    **restrict p_moment_ext_indx,
                             float  **restrict p_moment_ext_coef,
                             int verbose);

float 
fun_ricker(float t, float fc, float t0);

float
fun_gauss(float t, float a, float t0);

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

struct CubicPt *
LocaSrc(float sx, float sy, float sz,
        int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
        size_t siz_line, size_t siz_slice, size_t siz_volume, 
        float *restrict c3d,
        size_t *restrict c3d_pos,
        struct CubicPt *Pt);

struct SrcIndx 
CoorMap(float sx, float sy, float sz, 
        struct CubicPt  *Pt, 
        struct SrcIndx SrcInfro);

#endif
