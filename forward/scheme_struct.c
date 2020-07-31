
/*********************************************************************
 *
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "scheme_t.h"

//
// set grid size
//
int scheme_set_macdrp(struct scheme_t *fd)
{
  int ierr=0;

  fd->num_of_pairs = 8;
  // BBB FFB FFF BBF
  // BFB FBB FBF BFF
  int FD_Flags[8][3] = // 8 pairs, x/y/z 3 dim
  {
    { 0,  0,  0},
    { 1,  1,  0},
    { 1,  1,  1},
    { 0,  0,  1},
    { 0,  1,  0},
    { 1,  0,  0},
    { 1,  0,  1},
    { 0,  1,  1}
  }; //  0(B) 1(F)

  // alloc
  fd->pair_fdx_len = (int *) fdlib_mem_calloc_3d_int( 
                 fd->num_of_pairs, fd->num_rk_stages, 5 * fd->op_max_half_stentil,
                 0, "scheme_set_macdrp");
  fd->pair_fdx_indx = (int *) fdlib_mem_calloc_3d_size_t( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->op_all_coef_size, 0, "scheme_set_macdrp");
  fd->pair_fdx_coef = (int *) fdlib_mem_calloc_3d_flat( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->op_all_coef_size, 0.0, "scheme_set_macdrp");

  fd->pair_fdy_len = (int *) fdlib_mem_calloc_3d_int( 
                 fd->num_of_pairs , fd->num_rk_stages , 5 * fd->op_max_half_stentil,
                 0, "scheme_set_macdrp");
  fd->pair_fdy_indx = (int *) fdlib_mem_calloc_3d_size_t( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->op_all_coef_size, 0, "scheme_set_macdrp");
  fd->pair_fdy_coef = (int *) fdlib_mem_calloc_3d_flat( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->op_all_coef_size, 0.0, "scheme_set_macdrp");

  fd->pair_fdz_len = (int *) fdlib_mem_calloc_3d_int( 
                 fd->num_of_pairs , fd->num_rk_stages , 5 * fd->op_max_half_stentil,
                 0, "scheme_set_macdrp");
  fd->pair_fdz_indx = (int *) fdlib_mem_calloc_3d_size_t( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->op_all_coef_size, 0, "scheme_set_macdrp");
  fd->pair_fdz_coef = (int *) fdlib_mem_calloc_3d_flat( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->op_all_coef_size, 0.0, "scheme_set_macdrp");

  // set
  for (int ipair=0; ipair < fd->num_of_pairs; ipair++)
  {
    int idir0 = FD_Flags[ipair][0];
    int jdir0 = FD_Flags[ipair][1];
    int kdir0 = FD_Flags[ipair][2];

    for (int istage=0; istage < fd->num_rk_stages; istage++)
    {
      int idir = (idir0 + istage) % 2;
      int jdir = (jdir0 + istage) % 2;
      int kdir = (kdir0 + istage) % 2;

      // len
      for (int ih=0; ih < fd->op_max_half_stentil; ih++)
      {
        for (int i=0; i < 5; i++)
        {
          pair_fdx_len[ipar][istage][ih * t + i] = fd_dhs_len[idir][i][ih];
          pair_fdy_len[ipar][istage][ih * t + i] = fd_dhs_len[jdir][i][ih];
          pair_fdz_len[ipar][istage][ih * t + i] = fd_dhs_len[kdir][i][ih];
        }
      }

      // indx and coef
      for (int i=0; i < fd->op_all_coef_size; i++)
      {
        pair_fdx_len[ipar][istage][i] = fd_dhs_opindx[idir][i];
        pair_fdy_len[ipar][istage][i] = fd_dhs_opindx[jdir][i];
        pair_fdz_len[ipar][istage][i] = fd_dhs_opindx[kdir][i];

        pair_fdx_coef[ipar][istage][i] = fd_dhs_opcoef[idir][i];
        pair_fdy_coef[ipar][istage][i] = fd_dhs_opcoef[jdir][i];
        pair_fdz_coef[ipar][istage][i] = fd_dhs_opcoef[kdir][i];
      }
    }
  }

  pair_fdx_length = 
  
  return ierr;
}

