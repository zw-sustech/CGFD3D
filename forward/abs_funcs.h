#ifndef ABS_FUNCS_H
#define ABS_FUNCS_H

// cal  1d index for icmp,it,istage
//  with respect to the start pointer of this source point
//#define M_ABS_IND(icmp,it,istage,nt,num_stage) \
//  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

float
abs_cal_pml_R(int N);

float
abs_cal_pml_dmax(float L, float Vp, float Rpp);

float
abs_cal_pml_amax(float fc);

float
abs_cal_pml_d(float x, float L, float dmax);

float
abs_cal_pml_a(float x, float L, float amax);

float
abs_cal_pml_b(float x, float L, float bmax);

void
abs_set_cfspml(
    float *alpha_max, //
    float *beta_max, //
    float *velocity, //
    float dh, // average grid size, need to change use x3d etc to cal
    size_t ni1,
    size_t ni2,
    size_t nj1,
    size_t nj2,
    size_t nk1,
    size_t nk2,
    int *boundary_itype, // input
    int *abs_num_of_layers, // output
    size_t *abs_indx;
    size_t *abs_coefs_facepos0,
    float **p_abs_coefs,
    int verbose);

// set cfs pml vars
void
abs_init_vars_cfspml(
    int number_of_levels,
    int number_of_vars,
    int *restrict boundary_itype,
    size_t *restrict abs_indx,
    size_t *restrict abs_vars_volsiz,
    size_t *restrict abs_vars_facepos0,
    size_t *abs_vars_size_per_level,
    float *restrict p_abs_blk_vars
    const int myid, const int verbose);

#endif
