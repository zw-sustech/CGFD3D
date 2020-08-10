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
void 
fd_set_macdrp(struct fd_t *fd)
{
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

  fd->num_of_pairs = 8;
  // BBB FFB FFF BBF
  // BFB FBB FBF BFF

  // set max
  fd->fdx_max_len = 5;
  fd->fdy_max_len = 5;
  fd->fdz_max_len = 5;
  fd->fdx_nghosts = 3;
  fd->fdy_nghosts = 3;
  fd->fdz_nghosts = 3;
  fd->fdx_num_surf_lay = 0;
  fd->fdy_num_surf_lay = 0;
  fd->fdz_num_surf_lay = 3;

  // alloc
  fd->pair_fdx_all_info = (size_t ****) fdlib_mem_calloc_4l_sizet( 
                 fd->num_of_pairs, fd->num_rk_stages, fd->mac_max_half_stentil+1,FD_INFO_SIZE,
                 0, "fd_set_macdrp");
  fd->pair_fdx_all_indx = (size_t ***) fdlib_mem_calloc_3l_sizet( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");
  fd->pair_fdx_all_coef = (int ***) fdlib_mem_calloc_3l_float( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0.0, "fd_set_macdrp");

  fd->pair_fdy_all_info = (size_t ****) fdlib_mem_calloc_4l_sizet( 
                 fd->num_of_pairs, fd->num_rk_stages, fd->mac_max_half_stentil+1,FD_INFO_SIZE,
                 0, "fd_set_macdrp");
  fd->pair_fdy_all_indx = (size_t ***) fdlib_mem_calloc_3l_sizet( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");
  fd->pair_fdy_all_coef = (int ***) fdlib_mem_calloc_3l_float( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0.0, "fd_set_macdrp");

  fd->pair_fdz_all_info = (size_t ****) fdlib_mem_calloc_4l_sizet( 
                 fd->num_of_pairs, fd->num_rk_stages, fd->mac_max_half_stentil+1,FD_INFO_SIZE,
                 0, "fd_set_macdrp");
  fd->pair_fdz_all_indx = (size_t ***) fdlib_mem_calloc_3l_sizet( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0, "fd_set_macdrp");
  fd->pair_fdz_all_coef = (int ***) fdlib_mem_calloc_3l_float( 
                 fd->num_of_pairs , fd->num_rk_stages , fd->mac_all_coef_size, 0.0, "fd_set_macdrp");

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
}

//
// set grid size
//
void
fd_blk_init(struct fd_blk_t *blk,
      int number_of_x_points,
      int number_of_y_points,
      int number_of_z_points,
      int number_of_x_procs,
      int number_of_y_procs,
      int number_of_z_procs,
      char boundary_type_name[][],
      //char *abs_type_name;
      int *abs_number_of_layers,
      int fdx_nghosts,
      int fdy_nghosts,
      int fdz_nghosts,
      int number_of_levels, // depends on time scheme, for rk4 = 4
      int *myid2,
      int *neighid,
      const int myid, const int verbose)
{
  int boundary_itype[FD_NDIM_2];

  //
  // set boundary itype
  //
  //if (strcmp(abs_type_name, "cfspml")) abs_itype = FD_BOUNDARY_TYPE_CFSPML;
  //if (strcmp(abs_type_name, "ablexp")) abs_itype = FD_BOUNDARY_TYPE_ABLEXP;

  int abs_itype = FD_BOUNDARY_TYPE_NONE;

  for (int i=0; i<FD_NDIM_2; i++) {
    // default none
    int bdry_itype = FD_BOUNDARY_TYPE_NONE;
    blk->abs_number_of_layers[i] = 0;

    // set according input
    if ( (i==0 && myid2[0]==0) ||
         (i==1 && myid2[0]==number_of_x_procs -1) ||
         (i==2 && myid2[1] == 0) ||
         (i==3 && myid2[1] == number_of_y_procs -1) ||
         (i>=4))
    {
      // cfs-pml
      if (strcmp(boundary_type_name[i], "cfspml"))
      {
        bdry_itype     = FD_BOUNDARY_TYPE_CFSPML;
        blk->abs_itype = FD_BOUNDARY_TYPE_CFSPML;
        blk->abs_number_of_layers[i] = abs_number_of_layers[i];
        abs_itype = bdry_itype;
      }

      // free
      if (strcmp(boundary_type_name[i], "free"  )) {
        bdry_itype = FD_BOUNDARY_TYPE_FREE;
      }
    }

    blk->boundary_itype[i] = bdry_itype;
    blk->abs_itype = abs_itype;
  }

  // check this proc is boundary or not
  //if (neighid[0] == MPI_PROC_NULL) ;

  //
  // set grid size
  //

  // determine ni
  int nx_et = number_of_x_points;
  // double cfspml load
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML) {
    nx_et += abs_number_of_layers[0] + abs_number_of_layers[1];
  }
  int nx_avg  = nx_et / number_of_x_procs;
  int nx_left = nx_et % number_of_x_procs;
  if (nx_avg < 2 * fdx_nghosts) {
    // error
  }
  if (nx_avg<abs_number_of_layers[0] || nx_avg<abs_number_of_layers[1]) {
    // error
  }
  int ni = nx_avg;
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[0]==0) {
    ni -= abs_number_of_layers[0];
  }
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[0]==number_of_x_procs-1) {
    ni -= abs_number_of_layers[1];
  }
  // not equal divided points given to first nx_left procs
  if (myid2[0] < nx_left) {
    ni++;
  }
  // global index
  blk->gni1 = myid2[0] * nx_avg - abs_number_of_layers[0];
  if (nx_left != 0) {
    blk->gni1 += (myid2[0] < nx_left)? myid2[0] : nx_left;
  }

  // determine nj
  int ny_et = number_of_y_points;
  // double cfspml load
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML) {
    ny_et += abs_number_of_layers[2] + abs_number_of_layers[3];
  }
  int ny_avg  = ny_et / number_of_y_procs;
  int ny_left = ny_et % number_of_y_procs;
  if (ny_avg < 2 * fdy_nghosts) {
    // error
  }
  if (ny_avg<abs_number_of_layers[2] || ny_avg<abs_number_of_layers[3]) {
    // error
  }
  int nj = ny_avg;
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[1]==0) {
    nj -= abs_number_of_layers[2];
  }
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[1]==number_of_y_procs-1) {
    nj -= abs_number_of_layers[3];
  }
  // not equal divided points given to first ny_left procs
  if (myid2[1] < ny_left) {
    nj++;
  }
  // global index
  blk->gnj1 = myid2[1] * ny_avg - abs_number_of_layers[2];
  if (ny_left != 0) {
    blk->gnj1 += (myid2[1] < ny_left)? myid2[1] : ny_left;
  }

  // determine nk
  int nk = number_of_z_points;
  blk->gnk1 = 0;
  
  // add ghost points
  nx = ni + 2 * fdx_nghosts;
  ny = nj + 2 * fdy_nghosts;
  nz = nk + 2 * fdz_nghosts;

  blk->ni = ni;
  blk->nj = nj;
  blk->nk = nk;

  blk->nx = nx;
  blk->ny = ny;
  blk->nz = nz;

  blk->ni1 = fdx_nghosts;
  blk->ni2 = blk->ni1 + ni - 1;

  blk->nj1 = fdy_nghosts;
  blk->nj2 = blk->nj1 + nj - 1;

  blk->nk1 = fdz_nghosts;
  blk->nk2 = blk->nk1 + nk - 1;
  
  // x dimention varies first
  blk->siz_line   = nx; 
  blk->siz_slice  = nx * ny; 
  blk->siz_volume = nx * ny * nz;
  
  // level
  blk->w3d_num_of_vars = number_of_levels;

  return ierr;
}
