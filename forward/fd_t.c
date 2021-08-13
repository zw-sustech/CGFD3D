/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "fdlib_mem.h"
#include "fd_t.h"

/*
 * set MacCormack DRP scheme with rk 4
 */

int 
fd_set_macdrp(fd_t *fd)
{
  int ierr = 0;

  //----------------------------------------------------------------------------
  // 4th rk scheme
  //----------------------------------------------------------------------------

  fd->num_rk_stages = 4;

  fd->rk_a = (float *) fdlib_mem_malloc_1d(
                          fd->num_rk_stages*sizeof(float),"fd_set_macdrp");
  fd->rk_b = (float *) fdlib_mem_malloc_1d(
                          fd->num_rk_stages*sizeof(float),"fd_set_macdrp");
  fd->rk_rhs_time = (float *) fdlib_mem_malloc_1d(
                          fd->num_rk_stages*sizeof(float), "fd_set_macdrp");

  fd->rk_a[0] = 0.5;
  fd->rk_a[1] = 0.5;
  fd->rk_a[2] = 1.0;

  fd->rk_b[0] = 1.0/6.0;
  fd->rk_b[1] = 1.0/3.0;
  fd->rk_b[2] = 1.0/3.0;
  fd->rk_b[3] = 1.0/6.0;

  // rhs side time in terms of percentage time step, useful for source stf
  fd->rk_rhs_time[0] = 0.0;
  fd->rk_rhs_time[1] = 0.5;
  fd->rk_rhs_time[2] = 0.5;
  fd->rk_rhs_time[3] = 1.0;

  //----------------------------------------------------------------------------
  // MacCormack-type scheme
  //----------------------------------------------------------------------------

  // set max
  fd->fdx_max_len = 5;
  fd->fdy_max_len = 5;
  fd->fdz_max_len = 5;
  fd->fdx_max_half_len = 3;
  fd->fdy_max_half_len = 3;
  fd->fdz_max_half_len = 3;
  fd->fdx_nghosts = 3;
  fd->fdy_nghosts = 3;
  fd->fdz_nghosts = 3;

  //----------------------------------------------------------------------------
  // 1d scheme for different points to surface, or different half length
  //----------------------------------------------------------------------------
#define  m_mac_max_len   5
#define  m_mac_num_lay   4

  // op index at all layers, zero not used point
  int  mac_all_indx[m_mac_num_lay][2][m_mac_max_len] =
  {
    { // forw/back at free, no use
      {    0,1,0,0,0 }, // 0, at free, no use
      { -1,0,0,0,0 }
    },
    { // 1, 2nd
      {    0,1,0,0,0 },
      { -1,0,0,0,0 }
    }, 
    { // 2, 4th
      {        0,1,2,0,0 }, 
      { -2,-1, 0,0,0 }
    },
    { // 3, drp, normal op
      {       -1,0,1,2,3 },
      { -3,-2,-1,0,1 } 
    }
  };

  // op coef at all layers
  float mac_all_coef[m_mac_num_lay][2][m_mac_max_len] =
  {
    { // at free
      {     -1.0, 1.0, 0.0, 0.0, 0.0 },
      {-1.0, 1.0, 0.0, 0.0, 0.0}
    },
    { // 1, 2nd
      {     -1.0, 1.0, 0.0, 0.0, 0.0 },
      {-1.0, 1.0, 0.0, 0.0, 0.0}
    }, 
    { // 2, 4th
      {                   -7.0/6.0, 8.0/6.0, -1.0/6.0, 0.0, 0.0 }, 
      { 1.0/6.0, -8.0/6.0, 7.0/6.0, 0.0, 0.0 }
    },
    { // 3, drp, normal op
      {                   -0.30874,-0.6326 ,1.233 ,-0.3334,0.04168 },
      { -0.04168 , 0.3334, -1.233 , 0.6326 , 0.30874 }
    } 
  };

  int  mac_all_total_len[m_mac_num_lay][2] = 
  { 
    { 2, 2 },
    { 2, 2 },
    { 3, 3 },
    { 5, 5 }
  };

  // half len without cur point
  int  mac_all_half_len[m_mac_num_lay][2] = 
  { 
    { 1, 1 },
    { 1, 1 },
    { 2, 2 },
    { 3, 3 }
  };

  int  mac_all_left_len[m_mac_num_lay][2] = 
  { 
    { 0, 1 },
    { 0, 1 },
    { 0, 2 },
    { 1, 3 }
  };

  int  mac_all_right_len[m_mac_num_lay][2] = 
  { 
    { 1, 0 },
    { 1, 0 },
    { 2, 0 },
    { 3, 1 }
  };

  // center scheme for macdrp at diff layers, which is used in metric calculation
#define  m_mac_center_max_len   7

  int  mac_center_all_indx[m_mac_num_lay][m_mac_center_max_len] =
  {
    // forw/back at free, no use
    { -1,0,1,0,0,0,0 },
    // 1, 2nd
    { -1,0,1,0,0,0,0 },
    // 2, 4th
    { -2,-1,0,1,-2,0,0 },
    // 3, drp, normal op
    { -3,-2,-1,0,1,2,3 } 
  };

  // op coef at all layers
  float mac_center_all_coef[m_mac_num_lay][m_mac_center_max_len] =
  {
    { -0.5, 0.0, 0.5, 0,0,0,0 },
    { -0.5, 0.0, 0.5, 0,0,0,0 },
    { 0.08333333, -0.6666666, 0.0, 0.6666666, -0.08333333,0,0 },
    {-0.02084, 0.1667, -0.7709, 0.0, 0.7709, -0.1667, 0.02084 }
  };

  int  mac_center_all_total_len[m_mac_num_lay] = { 3, 3, 5, 7 };

  // half len without cur point
  int  mac_center_all_half_len[m_mac_num_lay] = { 1,1,2,3 };

  //----------------------------------------------------------------------------
  // combine to 3d op
  //----------------------------------------------------------------------------

  fd->num_of_pairs = 8;
  // BBB FFB FFF BBF
  // BFB FBB FBF BFF
  int FD_Flags[8][CONST_NDIM] = // 8 pairs, x/y/z 3 dim
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

  fd->num_of_fdx_op = 1;
  fd->num_of_fdy_op = 1;
  fd->num_of_fdz_op = 4;

  // alloc

  fd->pair_fdx_op = (fd_op_t ***)malloc(fd->num_of_pairs * sizeof(fd_op_t **));
  fd->pair_fdy_op = (fd_op_t ***)malloc(fd->num_of_pairs * sizeof(fd_op_t **));
  fd->pair_fdz_op = (fd_op_t ***)malloc(fd->num_of_pairs * sizeof(fd_op_t **));

  for (int ipair = 0; ipair < fd->num_of_pairs; ipair++)
  {
    fd->pair_fdx_op[ipair] = (fd_op_t **)malloc(fd->num_rk_stages * sizeof(fd_op_t *));
    fd->pair_fdy_op[ipair] = (fd_op_t **)malloc(fd->num_rk_stages * sizeof(fd_op_t *));
    fd->pair_fdz_op[ipair] = (fd_op_t **)malloc(fd->num_rk_stages * sizeof(fd_op_t *));

    for (int istage = 0; istage < fd->num_rk_stages; istage++)
    {
      fd->pair_fdx_op[ipair][istage] = (fd_op_t *)malloc(fd->num_of_fdx_op
                                                           * sizeof(fd_op_t));
      fd->pair_fdy_op[ipair][istage] = (fd_op_t *)malloc(fd->num_of_fdy_op
                                                           * sizeof(fd_op_t));
      fd->pair_fdz_op[ipair][istage] = (fd_op_t *)malloc(fd->num_of_fdz_op 
                                                           * sizeof(fd_op_t));
    }
  }

  // set
  for (int ipair=0; ipair < fd->num_of_pairs; ipair++)
  {
    int idir0 = FD_Flags[ipair][0];
    int jdir0 = FD_Flags[ipair][1];
    int kdir0 = FD_Flags[ipair][2];

    for (int istage=0; istage < fd->num_rk_stages; istage++)
    {
      // switch forw/back based on mod value
      int idir = (idir0 + istage) % 2;
      int jdir = (jdir0 + istage) % 2;
      int kdir = (kdir0 + istage) % 2;

      // for fdx
      fd_op_t *fdx_op = fd->pair_fdx_op[ipair][istage];
      fdx_op->total_len = mac_all_total_len[m_mac_num_lay-1][idir];
      fdx_op->half_len  = mac_all_half_len[m_mac_num_lay-1][idir];
      fdx_op->left_len  = mac_all_left_len[m_mac_num_lay-1][idir];
      fdx_op->right_len = mac_all_right_len[m_mac_num_lay-1][idir];

      fdx_op->indx = (int   *)malloc(fdx_op->total_len * sizeof(int));
      fdx_op->coef = (float *)malloc(fdx_op->total_len * sizeof(float));

      for (int n=0; n < fdx_op->total_len; n++)
      {
        fdx_op->indx[n] = mac_all_indx[m_mac_num_lay-1][idir][n];
        fdx_op->coef[n] = mac_all_coef[m_mac_num_lay-1][idir][n];
      }

      // for fdy
      fd_op_t *fdy_op = fd->pair_fdy_op[ipair][istage];
      fdy_op->total_len = mac_all_total_len[m_mac_num_lay-1][jdir];
      fdy_op->half_len  = mac_all_half_len[m_mac_num_lay-1][jdir];
      fdy_op->left_len  = mac_all_left_len[m_mac_num_lay-1][jdir];
      fdy_op->right_len = mac_all_right_len[m_mac_num_lay-1][jdir];

      fdy_op->indx = (int   *)malloc(fdy_op->total_len * sizeof(int));
      fdy_op->coef = (float *)malloc(fdy_op->total_len * sizeof(float));
      for (int n=0; n < fdy_op->total_len; n++)
      {
        fdy_op->indx[n] = mac_all_indx[m_mac_num_lay-1][jdir][n];
        fdy_op->coef[n] = mac_all_coef[m_mac_num_lay-1][jdir][n];
      }

      // for fdz at each layer near surface
      for (int ilay=0; ilay < fd->num_of_fdz_op; ilay++)
      {
        fd_op_t *fdz_op = fd->pair_fdz_op[ipair][istage] + ilay;

        fdz_op->total_len = mac_all_total_len[ilay][kdir];
        fdz_op->half_len  = mac_all_half_len[ilay][kdir];
        fdz_op->left_len  = mac_all_left_len[ilay][kdir];
        fdz_op->right_len = mac_all_right_len[ilay][kdir];

        fdz_op->indx = (int   *)malloc(fdz_op->total_len * sizeof(int));
        fdz_op->coef = (float *)malloc(fdz_op->total_len * sizeof(float));
        for (int n=0; n < fdz_op->total_len; n++)
        {
          fdz_op->indx[n] = mac_all_indx[ilay][kdir][n];
          fdz_op->coef[n] = mac_all_coef[ilay][kdir][n];
        }
      }

    } // istage
  } // ipair

  // centered op for metrics, only use normal order
  fd->fdc_len      = mac_center_all_total_len[m_mac_num_lay-1];
  fd->fdc_half_len = mac_center_all_half_len[m_mac_num_lay-1];
  fd->fdc_nghosts  = fd->fdc_half_len;

  fd->fdc_indx = (int   *) fdlib_mem_malloc_1d(fd->fdc_len*sizeof(int),"fdc_set_macdrp");
  fd->fdc_coef = (float *) fdlib_mem_malloc_1d(fd->fdc_len*sizeof(float),"fdc_set_macdrp");
  for (size_t i=0; i < fd->fdc_len; i++) {
    fd->fdc_indx[i] = mac_center_all_indx[m_mac_num_lay-1][i];
    fd->fdc_coef[i] = mac_center_all_coef[m_mac_num_lay-1][i];
  }

  return ierr;
}

void
fd_print(fd_t *fd)
{    
  fprintf(stdout, "\n-------------------------------------------------------\n");
  fprintf(stdout, "print fd scheme info:\n");
  fprintf(stdout, "-------------------------------------------------------\n\n");

  fprintf(stdout, " num_rk_stages = %d\n", fd->num_rk_stages);

  fprintf(stdout, "  rk_a = ");
  for (int i=0; i<fd->num_rk_stages-1; i++) {
    fprintf(stdout, " %e", fd->rk_a[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "  rk_b = ");
  for (int i=0; i<fd->num_rk_stages; i++) {
    fprintf(stdout, " %e", fd->rk_b[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, " fdc_len = %d\n", fd->fdc_len);
  fprintf(stdout, " fdc_half_len = %d\n", fd->fdc_half_len);
  fprintf(stdout, " fdc_nghosts = %d\n", fd->fdc_nghosts);
  fprintf(stdout, "  fdc_indx = ");
  for (int i=0; i<fd->fdc_len; i++) {
    fprintf(stdout, " %d", fd->fdc_indx[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "  fdc_coef = ");
  for (int i=0; i<fd->fdc_len; i++) {
    fprintf(stdout, " %e", fd->fdc_coef[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, " fdx_max_len = %d\n", fd->fdx_max_len);

  fprintf(stdout, "  others info need to be coded\n");
  // rest
}

/*
 * set staggered-grid scheme 
 */

int 
fd_set_stg4(fdstg_t *fd)
{
  int ierr = 0;

  // set max
  fd->fdx_nghosts = 3;
  fd->fdy_nghosts = 3;
  fd->fdz_nghosts = 3;
  fd->fdx_max_len = 5;
  fd->fdy_max_len = 5;
  fd->fdz_max_len = 5;
  fd->fdx_max_half_len = 3;
  fd->fdy_max_half_len = 3;
  fd->fdz_max_half_len = 3;

  //----------------------------------------------------------------------------
  // 1d scheme for different points to surface, or different half length
  //----------------------------------------------------------------------------
#define  m_stg_max_len   4
#define  m_stg_num_lay   2

  // op index at all layers, zero not used point
  int  stg_all_indx[m_stg_num_lay][m_stg_max_len] =
  {
    { -1, 0, 0,0 },
    { -2,-1, 0,1 }
  };

  // op coef at all layers
  float stg_all_coef[m_stg_num_lay][m_stg_max_len] =
  {
    { -1.0, 1.0, 0,0 },
    { 0.0416666667, -1.125, 1.125, -0.0416666667 }
  };

  int  stg_all_total_len[m_stg_num_lay] = 
  { 2, 4 };

  // half len without cur point
  int  stg_all_half_len[m_stg_num_lay] = 
  { 1, 2 };

  int  stg_all_left_len[m_stg_num_lay] = 
  { 1, 2 }; 

  int  stg_all_right_len[m_stg_num_lay] = 
  { 0, 1 };

  //----------------------------------------------------------------------------
  // combine to fd op
  //----------------------------------------------------------------------------

  fd->num_of_fdx_op = 1;
  fd->num_of_fdy_op = 1;
  fd->num_of_fdz_op = 2;

  // alloc

  fd->lay_fdx_op = (fd_op_t *)malloc(fd->num_of_fdx_op * sizeof(fd_op_t));
  fd->lay_fdy_op = (fd_op_t *)malloc(fd->num_of_fdy_op * sizeof(fd_op_t));
  fd->lay_fdz_op = (fd_op_t *)malloc(fd->num_of_fdz_op * sizeof(fd_op_t));

  // set fdx
  fd_op_t *fdx_op = fd->lay_fdx_op;
  fdx_op->total_len = stg_all_total_len[m_stg_num_lay-1];
  fdx_op->half_len  = stg_all_half_len [m_stg_num_lay-1];
  fdx_op->left_len  = stg_all_left_len [m_stg_num_lay-1];
  fdx_op->right_len = stg_all_right_len[m_stg_num_lay-1];

  fdx_op->indx = (int   *)malloc(fdx_op->total_len * sizeof(int));
  fdx_op->coef = (float *)malloc(fdx_op->total_len * sizeof(float));

  for (int n=0; n < fdx_op->total_len; n++)
  {
    fdx_op->indx[n] = stg_all_indx[m_stg_num_lay-1][n];
    fdx_op->coef[n] = stg_all_coef[m_stg_num_lay-1][n];
  }

  // set fdy
  fd_op_t *fdy_op = fd->lay_fdy_op;
  fdy_op->total_len = stg_all_total_len[m_stg_num_lay-1];
  fdy_op->half_len  = stg_all_half_len [m_stg_num_lay-1];
  fdy_op->left_len  = stg_all_left_len [m_stg_num_lay-1];
  fdy_op->right_len = stg_all_right_len[m_stg_num_lay-1];

  fdy_op->indx = (int   *)malloc(fdy_op->total_len * sizeof(int));
  fdy_op->coef = (float *)malloc(fdy_op->total_len * sizeof(float));

  for (int n=0; n < fdy_op->total_len; n++)
  {
    fdy_op->indx[n] = stg_all_indx[m_stg_num_lay-1][n];
    fdy_op->coef[n] = stg_all_coef[m_stg_num_lay-1][n];
  }

  // for fdz at each layer near surface
  for (int ilay=0; ilay < fd->num_of_fdz_op; ilay++)
  {
    fd_op_t *fdz_op = fd->lay_fdz_op + ilay;

    fdz_op->total_len = stg_all_total_len[ilay];
    fdz_op->half_len  = stg_all_half_len [ilay];
    fdz_op->left_len  = stg_all_left_len [ilay];
    fdz_op->right_len = stg_all_right_len[ilay];

    fdz_op->indx = (int   *)malloc(fdz_op->total_len * sizeof(int));
    fdz_op->coef = (float *)malloc(fdz_op->total_len * sizeof(float));
    for (int n=0; n < fdz_op->total_len; n++)
    {
      fdz_op->indx[n] = stg_all_indx[ilay][n];
      fdz_op->coef[n] = stg_all_coef[ilay][n];
    }
  }

  return ierr;
}

