#ifndef FD_T_H
#define FD_T_H

#include <stdlib.h>

// use siz_shift to find adjacent point of the stentil for 3d var
#define M_FD_SHIFT(deriv, var, iptr, fd_length, fd_shift, fd_coef, n) \
   deriv = 0.0; \
   for (n=0; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n]]; \
   }

// assume var has the same size as fd_coef, ordered one by one, thus no index needed
#define M_FD_NOINDX(deriv, var, fd_length, fd_coef, n) \
   deriv = 0.0; \
   for (n=0; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[n]; \
   }

// use indx relative to cur point as (-1,0,1), need to multiply siz_shift for 3d array
#define M_FD_INDX(deriv, var, iptr, fd_length, fd_indx, fd_coef, shift, n) \
   deriv = 0.0; \
   for (n=0; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n] * shift]; \
   }

#define FD_NDIM   3
#define FD_NDIM_2 6 // 2 * ndim
#define FD_MAX_STRLEN 1000

#define FD_INFO_POS_OF_INDX  0
#define FD_INFO_LENGTH_TOTAL 1
#define FD_INFO_LENGTH_HALF  2
#define FD_INFO_LENGTH_LEFT  3
#define FD_INFO_LENGTH_RIGTH 4
#define FD_INFO_SIZE         5

#define FD_BOUNDARY_TYPE_FREE   0
#define FD_BOUNDARY_TYPE_CFSPML 1
#define FD_BOUNDARY_TYPE_ABLEXP 1

/*
 * for different fd schemes
 */

struct fd_t{

  //
  // MacCormack-type scheme
  //

  // 1d scheme for different points to surface, or different half length
  
  const int mac_max_half_stentil = 3;
  const int mac_all_coef_size    = 10;

  // forw/back, free at 0/1/2/3 point to free, indx/total/half/left/right
  // half, left, right number will be used to pack message
  size_t mac_all_info[2][mac_max_half_stentil+1][FD_INFO_SIZE] =
  {
    { // forw
      // POS, TOTAL, HALF, LEFT, RIGHT
      { 0, 0, 0, 0, 0 }, // 0 free, not used
      { 0, 2, 1, 0, 1 }, // 1 to free
      { 2, 3, 2, 0, 2 }, // 2
      { 5, 5, 3, 1, 3 }  // 3, normal op for cur scheme
    },
    { // back
      { 0, 0, 0, 0, 0 }, // 0 free, not used
      { 0, 2, 1, 1, 0 }, // 1 to free
      { 2, 3, 2, 2, 0 }, // 2
      { 5, 5, 3, 3, 1 }  // 3, normal op for cur scheme
    }
  };

  size_t mac_all_indx[2][mac_all_coef_size] =
  {
    { // fowrd
      0, 1,
      0, 1, 2,
      -1, 0, 1, 2, 3
    },
    { // back
      -1, 0, 
      -2, -1, 0,
      -3,-2,-1,0,1
    } 
  };

  float mac_all_coef[2][mac_all_coef_size] = 
  {
    {
      -1.0, 1.0,
      -7.0/6.0, 8.0/6.0, -1.0/6.0, 
      -0.30874,-0.6326 ,1.233 ,-0.3334,0.04168 
    },
    {
      1.0, -1.0,
      1.0/6.0, -8.0/6.0, 7.0/6.0, 
      -0.04168 , 0.3334, -1.233 , 0.6326 , 0.30874
    } 
  };

  // center scheme for macdrp
  size_t mac_center_dhs_indx[3]   = { 0, 3, 8 };
  int    mac_center_dhs_lentot[3] = { 3, 5, 7 }; // total len
  int    mac_center_dhs_lenhal[3] = { 2, 3, 4 }; // left len
  float  mac_center_dhs_opindx[10] =
  {
    -1, 0, 1,
    -2,-1,0, 1, 2,
    -3,-2,-1, 0, 1, 2, 3
  };
  float  mac_center_dhs_opcoef[10] =
  {
    -0.5, 0.0, 0.5,
    -7.0/6.0, 8.0/6.0, -1.0/6.0,  // not finish
    -0.02084, 0.1667, -0.7709, 0.0, 0.7709, -0.1667, 0.02084 // need to check
  };

  //
	// pairs for 3d space for MacCormack-type schemes
  //

	int num_of_pairs;

  size_t ****pair_fdx_all_info;  // [pair][stage][k2free][pos, total, half, left, right]
  size_t  ***pair_fdx_all_indx;  // [pair][stage][1-10 ele], not include [hlen] due to different length
  float   ***pair_fdx_all_coef;

  size_t ****pair_fdy_all_info; 
  size_t  ***pair_fdy_all_indx; 
  float   ***pair_fdy_all_coef;

  size_t ****pair_fdz_all_info;
  size_t  ***pair_fdz_all_indx;
  float   ***pair_fdz_all_coef;

  //
  // center schemes at different points to boundary
  //
  int    **fdx_all_info; // [k2free][pos, total, half, left, right]
  size_t **fdx_all_indx;
  float  **fdx_all_coef;

  int    **fdy_all_info;
  size_t **fdy_all_indx;
  float  **fdy_all_coef;

  int    **fdz_all_info;
  size_t **fdz_all_indx;
  float  **fdz_all_coef;

  //
  // filter schemes at different points to boundary
  //
  int    **filtx_all_info; // [k2free][pos, total, half, left, right] 
  size_t **filtx_all_indx;
  float  **filtx_all_coef;

  int    **filty_all_info;
  size_t **filty_all_indx;
  float  **filty_all_coef;

  int    **filtz_all_info;
  size_t **filtz_all_indx;
  float  **filtz_all_coef;

  //
  // Runge-Kutta time scheme
  //
  int num_rk_stages = 4;
  float rk_a[4] = { 0., 0.5, 0.5, 1.0 };
  float rk_b[4] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
};

//
// block structure
//
struct fd_blk_t{

    // size of a single var
    size_t siz_line;
    size_t siz_slice;
    size_t siz_volume; // number of points per var

    //size_t siz_vars; // volume * num_of_vars, not easy for understand, may named with w3d and aux

    //
    // grid index
    //
    int ni1, ni2, nj1, nj2, nk1, nk2, ni, nj, nk;
    int nx1, nx2, ny1, ny2, nz1, nz2, nx, ny, nz;

    //
    // grid metrics
    //  x3d, y3d, z3d, jac, xi_x, etc
    float  *g3d; // grid 3d vars
    int     g3d_num_of_vars;
    size_t *g3d_pos; // for each var
    char  **g3d_name;

    //
    // media
    //  rho, lambda, mu etc
    float  *m3d; // media 3d vars
    int     m3d_num_of_vars;
    size_t *m3d_pos;
    char  **m3d_name;

    // Wavefield
    //
    float  *w3d; // wavefield 3d vars
    //size_t pos_vx, pos_vz, pos_txx, pos_tyy
    int     w3d_num_of_vars;
    int     w3d_num_of_levels;;
    size_t  w3d_size_per_level;
    size_t *w3d_pos; // for vars in a single level
    char  **w3d_name;

    // auxiliary vars
    //
    float  *a3d; // auxiliary vars
    int     a3d_num_of_vars;
    int     a3d_num_of_levels;
    size_t  a3d_size_per_level;
    size_t *a3d_pos;
    char  **a3d_name;

    // boundary 
    int *boundary_itype[ FD_NDIM_2 ];

    // free surface
    float *matVx2Vz, *matVy2Vz;

    // source term
    int num_of_force;
    int num_of_moment;

    //
    // abs
    //

    // only support one type each run, may combine in the future
    //  which requires different vars in blk_t
    int  abs_itype[3];

    int     abs_num_of_layers[ FD_NDIM_2 ];

    //  _blk means all dims, eliminate _dimpos var
    size_t  abs_blk_indx[FD_NDIM_2][FD_NDIM_2];

    // may add abs_numb_of_vars later for mpml

    //size_t  abs_blk_coefs_size[FD_NDIM_2]; // size of all coef for one block, seems no use
    //size_t  abs_coefs_dimpos[2*3];

    // abs_blk_coefs[idim] + abs_num_of_layers[idim] for next coef
    float   **abs_blk_coefs; // [dim][A,B,D etc]

    //size_t  abs_blk_vars_size_per_level[FD_NDIM_2]; // 

    size_t  abs_blk_vars_level_size; // size of vars per level
    size_t  abs_blk_vars_siz_volume[FD_NDIM_2]; // size of single var in abs_blk_vars for each block
    size_t  abs_blk_vars_blkpos[FD_NDIM_2]; // start pos of each blk in first level; other level similar
    float   *abs_blk_vars; //  order: vars_one_block -> all block -> one level -> 4 levels
    //size_t  abs_vars_dimpos[2*3];

    //
    // connection to other blk or mpi thread
    //
    size_t size_of_buff;
    int *inn_bdry_blk_id; // blk id
    size_t *inn_bdry_blk_indx_pair; // this blk indx to bdry_blk indx

    int out_mpi_num_of_neig; //
    int    *out_mpi_neig_ids; //
    size_t *out_mpi_neig_buff_size; //
    int *out_mpi_neigh_blk_ids; // mpi id and blk id
    size_t *out_bdry_blk


    // mem usage
    size_t number_of_float;
    size_t number_of_btye;
};

#endif
