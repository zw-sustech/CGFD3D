#ifndef SCHEME_STRUCT_H
#define SCHEME_STRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct scheme_struct{

  //
  // MacCormack-type scheme
  //

  // forw/back op
  //
  // dhs: different half stentil
  // save different op into a 1d array, use _indx[forw_dhs_lenhal-1] to retrieve the op in funcs
  int    forw_dhs_lentota[3] = { 2, 3, 5 }; // total len
  int    forw_dhs_lenhalf[3] = { 1, 2, 3 }; // half len exclude cur point
  int    forw_dhs_lenleft[3] = { 0, 0, 1 }; // left len
  int    forw_dhs_lenrigh[3] = { 1, 2, 3 }; // right len
  size_t forw_dhs_indx   [3] = { 0, 2, 5 }; 
  size_t forw_dhs_opindx[10] =
  {
    0, 1,
    0, 1, 2,
    -1, 0, 1, 2, 3
  };
  float  forw_dhs_opcoef[10] =
  {
    -1.0, 1.0,
    -7.0/6.0, 8.0/6.0, -1.0/6.0, 
    -0.30874,-0.6326 ,1.233 ,-0.3334,0.04168 
  };

  int    back_dhs_lentota[3] = { 2, 3, 5 }; // total len
  int    back_dhs_lenhalf[3] = { 1, 2, 3 }; // half len including cur point
  int    back_dhs_lenleft[3] = { 1, 2, 3 }; // left len
  int    back_dhs_lenrigh[3] = { 0, 0, 1 }; // right len
  size_t back_dhs_indx   [3] = { 0, 2, 5 };
  size_t back_dhs_opindx[10] =
  {
    -1, 0, 
    -2, -1, 0,
    -3,-2,-1,0,1
  };
  float  back_dhs_opcoef[10] =
  {
    1.0, -1.0,
    1.0/6.0, -8.0/6.0, 7.0/6.0, 
    -0.04168 , 0.3334, -1.233 , 0.6326 , 0.30874
  };

  //
	// pairs for 3d space
  //

	int num_of_pairs;

  int    *pair_fdx_length; // total, half, left, right
  size_t *pair_fdx_indx;  // 
  float  *pair_fdx_coef ;

  int    *pair_fdy_length; //
  size_t *pair_fdy_indx;  // 
  float  *pair_fdy_coef ;

  int    *pair_fdz_length; //
  size_t *pair_fdz_indx;  // 
  float  *pair_fdz_coef ;

  int    *pair_fdz_dhs_lenmax; // number of pair
  int    *pair_fdz_dhs_lentota; //
  int    *pair_fdz_dhs_lenhalf; //
  int    *pair_fdz_dhs_lenleft; //
  int    *pair_fdz_dhs_lenrigh; //
  size_t *pair_fdz_dhs_indx;  // 
  size_t *pair_fdz_dhs_opindx ;
  float  *pair_fdz_dhs_opcoef ;

  // center scheme
  int    *fdx_length; // total, helf length
  size_t *fdx_pos_shift;
  float  *fdx_coef;

  int    *fdy_length; // total, helf length
  size_t *fdy_pos_shift;
  float  *fdy_coef;

  int    *fdz_length; // total, helf length
  size_t *fdz_pos_shift;
  float  *fdz_coef;

  float mac_coef_forw_all[3][5] = 
  {
    { -1.0    , 1.0    , 0.0     , 0.0   , 0.0    },
    { -7.0/6.0, 8.0/6.0, -1.0/6.0, 0.0   , 0.0    },
    { -0.30874,-0.6326 ,1.233    ,-0.3334,0.04168 }
  }
  int mac_indx_forw_all[3][5] = 
  {
    { 0    , 1      , 0     , 0   , 0    },
    { 0    , 1,  2, 0   , 0    },
    { -1,0,1,2,3 }
  }

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
  // cent scheme plus explicit filter
  //
  float filt_dhs_opcoef;

  //
  // Runge-Kutta time scheme
  //
  int num_rk_stages = 4;
  float rk_a[4] = { 0., 0.5, 0.5, 1.0 };
  float rk_b[4] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };
};

#endif
