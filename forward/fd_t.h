#ifndef FD_T_H
#define FD_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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

#endif
