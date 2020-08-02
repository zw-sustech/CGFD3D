/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "fd_t.h"
#include "fdlib_mem.h"

//
// set grid size
//
int fd_set_macdrp(struct fd_t *fd)
{
  int ierr=0;

  fd->num_of_pairs = 8;
  // BBB FFB FFF BBF
  // BFB FBB FBF BFF
  int FD_Flags[8][FD_NDIM] = // 8 pairs, x/y/z 3 dim
  {
    { 0,  0,  0 },
    { 1,  1,  0 },
    { 1,  1,  1 },
    { 0,  0,  1 },
    { 0,  1,  0 },
    { 1,  0,  0 },
    { 1,  0,  1 },
    { 0,  1,  1 }
  }; //  0(B) 1(F)

  // alloc
  fd->pair_fdx_all_info = (size_t ****) fdlib_mem_calloc_4d_size_t( 
                 fd->num_of_pairs, fd->num_rk_stages, fd->mac_max_half_stentil+1,FD_INFO_SIZE,
                 0, "fd_set_macdrp");
  fd->pair_fdx_all_indx = (size_t ***) fdlib_mem_calloc_3d_size_t( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");
  fd->pair_fdx_all_coef = (int ***) fdlib_mem_calloc_3d_float( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");

  fd->pair_fdy_all_info = (size_t ****) fdlib_mem_calloc_4d_size_t( 
                 fd->num_of_pairs, fd->num_rk_stages, fd->mac_max_half_stentil+1,FD_INFO_SIZE,
                 0, "fd_set_macdrp");
  fd->pair_fdy_all_indx = (size_t ***) fdlib_mem_calloc_3d_size_t( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");
  fd->pair_fdy_all_coef = (int ***) fdlib_mem_calloc_3d_float( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");

  fd->pair_fdz_all_info = (size_t ****) fdlib_mem_calloc_4d_size_t( 
                 fd->num_of_pairs, fd->num_rk_stages, fd->mac_max_half_stentil+1,FD_INFO_SIZE,
                 0, "fd_set_macdrp");
  fd->pair_fdz_all_indx = (size_t ***) fdlib_mem_calloc_3d_size_t( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");
  fd->pair_fdz_all_coef = (int ***) fdlib_mem_calloc_3d_float( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");

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

      // info
      for (int ipoint=0; ipoint <= fd->op_max_half_stentil; ipoint++)
      {
        for (int i=0; i < FD_INFO_SIZE; i++)
        {
          pair_fdx_all_info[ipair][istage][ipoint][i] = mac_all_info[idir][ipoint][i];
          pair_fdy_all_info[ipair][istage][ipoint][i] = mac_all_info[jdir][ipoint][i];
          pair_fdz_all_info[ipair][istage][ipoint][i] = mac_all_info[kdir][ipoint][i];
          //*(*(*(pair_fdx_len+ipair)+istage)+ih * 5 + i) = fd_dhs_len[idir][i][ih];
          //*(*(*(pair_fdy_len+jpair)+istage)+ih * 5 + i) = fd_dhs_len[idir][j][ih];
          //*(*(*(pair_fdz_len+kpair)+istage)+ih * 5 + i) = fd_dhs_len[idir][k][ih];
        }
      }

      // indx and coef
      for (int i=0; i < fd->mac_all_coef_size; i++)
      {
        pair_fdx_all_indx[ipair][istage][i] = mac_all_indx[idir][i];
        pair_fdy_all_indx[ipair][istage][i] = mac_all_indx[jdir][i];
        pair_fdz_all_indx[ipair][istage][i] = mac_all_indx[kdir][i];

        pair_fdx_all_coef[ipair][istage][i] = mac_all_coef[idir][i];
        pair_fdy_all_coef[ipair][istage][i] = mac_all_coef[jdir][i];
        pair_fdz_all_coef[ipair][istage][i] = mac_all_coef[kdir][i];
      }
    }
  }

  return ierr;
}

