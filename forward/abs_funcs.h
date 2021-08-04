#ifndef ABS_FUNCS_H
#define ABS_FUNCS_H

// cal  1d index for icmp,it,istage
//  with respect to the start pointer of this source point
//#define M_ABS_IND(icmp,it,istage,nt,num_stage) \
//  ((icmp) * (nt) * (num_stage) + (it) * (num_stage) + (istage))

struct bdrypml_t
{
  // only support one type each run, may combine in the future
  //  which requires different vars in blk_t
  int  abs_itype;
  
  // only for this block, may diff with global values from input par
  int  abs_num_of_layers[ CONST_NDIM_2 ];
  
  // grid index of each face
  int    abs_indx[CONST_NDIM_2 * CONST_NDIM_2];
  
  int      abs_coefs_facepos0[CONST_NDIM_2];  // 
  float   *abs_coefs; // all coefs all faces 
  
  // may add abs_numb_of_vars later for mpml
  //   may not need it since 
  
  int     abs_vars_volsiz[CONST_NDIM_2]; // size of single var in abs_blk_vars for each block
  int     abs_vars_facepos0[CONST_NDIM_2]; // start pos of each blk in first level; other level similar
  int     abs_vars_size_per_level; // size of vars per level
  float   *abs_vars; //  order: vars_one_block -> all block -> one level -> 4 levels
};

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
    int    ni1,
    int    ni2,
    int    nj1,
    int    nj2,
    int    nk1,
    int    nk2,
    int *boundary_itype, // input
    int *abs_num_of_layers, // output
    int    *abs_indx,
    int    *abs_coefs_facepos0,
    float **p_abs_coefs,
    int verbose);

// set cfs pml vars
void
abs_init_vars_cfspml(
    int number_of_levels,
    int number_of_vars,
    int *restrict boundary_itype,
    int    *restrict abs_indx,
    int    *restrict abs_vars_volsiz,
    int    *restrict abs_vars_facepos0,
    int    *abs_vars_size_per_level,
    float **restrict p_abs_vars,
    const int myid, const int verbose);

#endif
