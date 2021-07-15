/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fd_t.h"
#include "fdlib_mem.h"

//
// set grid size
//
void 
fd_set_macdrp(struct fd_t *fd)
{
  //----------------------------------------------------------------------------
  // rk scheme
  //----------------------------------------------------------------------------
  fd->num_rk_stages = 4;

  fd->rk_a = (float *) fdlib_mem_malloc_1d(fd->num_rk_stages*sizeof(float),"fd_set_macdrp");
  fd->rk_b = (float *) fdlib_mem_malloc_1d(fd->num_rk_stages*sizeof(float),"fd_set_macdrp");
  fd->rk_rhs_time = (float *) fdlib_mem_malloc_1d(fd->num_rk_stages*sizeof(float),
                                                                           "fd_set_macdrp");

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
  fd->fdx_num_surf_lay = 3;
  fd->fdy_num_surf_lay = 3;
  fd->fdz_num_surf_lay = 3;

  // 1d scheme for different points to surface, or different half length
#define  m_mac_max_half_stentil   3
#define  m_mac_all_coef_size      10

  // forw/back, free at 0/1/2/3 point to free, indx/total/half/left/right
  // half, left, right number will be used to pack message
  int    mac_all_info[2][m_mac_max_half_stentil+1][FD_INFO_SIZE] =
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

  int    mac_all_indx[2][m_mac_all_coef_size] =
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

  float mac_all_coef[2][m_mac_all_coef_size] = 
  {
    {
      -1.0, 1.0,
      -7.0/6.0, 8.0/6.0, -1.0/6.0, 
      -0.30874,-0.6326 ,1.233 ,-0.3334,0.04168 
    },
    {
      -1.0, 1.0,
      1.0/6.0, -8.0/6.0, 7.0/6.0, 
      -0.04168 , 0.3334, -1.233 , 0.6326 , 0.30874
    } 
  };

  // center scheme for macdrp, which is used in metric calculation
#define m_mac_center_all_coef_size 15
  int    mac_center_all_info[m_mac_max_half_stentil+1][FD_INFO_SIZE] =
    {
      // POS, TOTAL, HALF, LEFT, RIGHT
      { 0, 0, 0, 0, 0 }, // 0 free, not used
      { 0, 3, 1, 1, 1 }, // 1 to free
      { 3, 5, 2, 2, 2 }, // 2
      { 8, 7, 3, 3, 3 }  // 3, normal op for cur scheme
    };
  float mac_center_all_indx[m_mac_center_all_coef_size] =
  {
    -1, 0, 1,
    -2,-1,0, 1, 2,
    -3,-2,-1, 0, 1, 2, 3
  };
  float mac_center_all_coef[m_mac_center_all_coef_size] =
  {
    -0.5, 0.0, 0.5,
    // 2-4 scheme
    //                   -7.0/6.0, 8.0/6.0, -1.0/6.0, 
    // 1.0/6.0, -8.0/6.0, 7.0/6.0, 
    0.08333333, -0.6666666, 0.0, 0.6666666, -0.08333333,
    // drp scheme
    //                   -0.30874,-0.6326 ,1.233 ,-0.3334,0.04168 
    // -0.04168 , 0.3334, -1.233 , 0.6326 , 0.30874
    // 
    -0.02084, 0.1667, -0.7709, 0.0, 0.7709, -0.1667, 0.02084 
  };

  // for 3d op
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
  fd->pair_fdx_all_info = (int    ****) fdlib_mem_calloc_4l_int( 
                                              fd->num_of_pairs,
                                              fd->num_rk_stages,
                                              m_mac_max_half_stentil+1,
                                              FD_INFO_SIZE,
                                              0,
                                              "fd_set_macdrp");
  fd->pair_fdx_all_indx = (int    ***) fdlib_mem_calloc_3l_int( 
                                              fd->num_of_pairs ,
                                              fd->num_rk_stages ,
                                              m_mac_all_coef_size,
                                              0,
                                              "fd_set_macdrp");
  fd->pair_fdx_all_coef = (float ***) fdlib_mem_calloc_3l_float( 
                                              fd->num_of_pairs ,
                                              fd->num_rk_stages ,
                                              m_mac_all_coef_size,
                                              0.0,
                                              "fd_set_macdrp");

  fd->pair_fdy_all_info = (int    ****) fdlib_mem_calloc_4l_int( 
                                              fd->num_of_pairs,
                                              fd->num_rk_stages,
                                              m_mac_max_half_stentil+1,
                                              FD_INFO_SIZE,
                                              0,
                                              "fd_set_macdrp");
  fd->pair_fdy_all_indx = (int    ***) fdlib_mem_calloc_3l_int( 
                                              fd->num_of_pairs ,
                                              fd->num_rk_stages ,
                                              m_mac_all_coef_size,
                                              0,
                                              "fd_set_macdrp");
  fd->pair_fdy_all_coef = (float ***) fdlib_mem_calloc_3l_float( 
                                              fd->num_of_pairs ,
                                              fd->num_rk_stages ,
                                              m_mac_all_coef_size,
                                              0.0,
                                              "fd_set_macdrp");

  fd->pair_fdz_all_info = (int    ****) fdlib_mem_calloc_4l_int( 
                                              fd->num_of_pairs,
                                              fd->num_rk_stages,
                                              m_mac_max_half_stentil+1,
                                              FD_INFO_SIZE,
                                              0,
                                              "fd_set_macdrp");
  fd->pair_fdz_all_indx = (int    ***) fdlib_mem_calloc_3l_int( 
                                              fd->num_of_pairs ,
                                              fd->num_rk_stages ,
                                              m_mac_all_coef_size,
                                              0,
                                              "fd_set_macdrp");
  fd->pair_fdz_all_coef = (float ***) fdlib_mem_calloc_3l_float( 
                                              fd->num_of_pairs ,
                                              fd->num_rk_stages ,
                                              m_mac_all_coef_size,
                                              0.0,
                                              "fd_set_macdrp");

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

      // for each point near surface
      for (int ipoint=0; ipoint <= m_mac_max_half_stentil; ipoint++)
      {
        for (int i=0; i < FD_INFO_SIZE; i++)
        {
          fd->pair_fdx_all_info[ipair][istage][ipoint][i] = mac_all_info[idir][ipoint][i];
          fd->pair_fdy_all_info[ipair][istage][ipoint][i] = mac_all_info[jdir][ipoint][i];
          fd->pair_fdz_all_info[ipair][istage][ipoint][i] = mac_all_info[kdir][ipoint][i];
        }
      }

      // indx and coef
      for (int i=0; i < m_mac_all_coef_size; i++)
      {
        fd->pair_fdx_all_indx[ipair][istage][i] = mac_all_indx[idir][i];
        fd->pair_fdy_all_indx[ipair][istage][i] = mac_all_indx[jdir][i];
        fd->pair_fdz_all_indx[ipair][istage][i] = mac_all_indx[kdir][i];

        fd->pair_fdx_all_coef[ipair][istage][i] = mac_all_coef[idir][i];
        fd->pair_fdy_all_coef[ipair][istage][i] = mac_all_coef[jdir][i];
        fd->pair_fdz_all_coef[ipair][istage][i] = mac_all_coef[kdir][i];
      }
    }
  }

  // centered op for metrics
  int fd_len   = mac_center_all_info[3][1];
  int fd_pos   = mac_center_all_info[3][0];

  fd->fd_indx = (int *) fdlib_mem_malloc_1d(fd_len*sizeof(int),"fd_set_macdrp");
  fd->fd_coef = (float  *) fdlib_mem_malloc_1d(fd_len*sizeof(float),"fd_set_macdrp");
  for (size_t i=0; i<fd_len; i++) {
    fd->fd_indx[i] = mac_center_all_indx[fd_pos+i];
    fd->fd_coef[i] = mac_center_all_coef[fd_pos+i];
  }
  fd->fd_len = fd_len;
  fd->fd_half_len = mac_center_all_info[3][2];
  fd->fd_nghosts  = fd->fd_half_len;
}

//
// set grid size
//
void
fd_blk_init(struct fd_blk_t *blk,
            int number_of_total_grid_points_x,
            int number_of_total_grid_points_y,
            int number_of_total_grid_points_z,
            int number_of_mpiprocs_x,
            int number_of_mpiprocs_y,
            char **boundary_type_name,
            int *abs_num_of_layers,
            char *output_dir,
            char *grid_export_dir,
            char *media_export_dir,
            int fdx_nghosts,
            int fdy_nghosts,
            int fdz_nghosts,
            int number_of_levels, // depends on time scheme, for rk4 = 4
            MPI_Comm comm, 
            const int myid, const int verbose)
{
  int boundary_itype[FD_NDIM_2];

  // set name
  //sprintf(blk->name, "%s", name);

  //
  // mpi topo
  //

  int pdims[2]={number_of_mpiprocs_x,number_of_mpiprocs_y};
  int periods[2] = {0,0};
  int rank; 
  // create Cartesian topology
  MPI_Cart_create(comm, 2, pdims, periods, 0, &blk->topocomm);
  // get my local x,y coordinates
  MPI_Cart_coords(blk->topocomm, myid, 2, blk->myid2);

  // neighour
  MPI_Cart_shift(blk->topocomm, 0, 1, &(blk->neighid[0]), &(blk->neighid[1]));
  MPI_Cart_shift(blk->topocomm, 1, 1, &(blk->neighid[2]), &(blk->neighid[3]));

  int myid2[]   = { blk->myid2[0], blk->myid2[1] };
  int neighid[] = { blk->neighid[0], 
                    blk->neighid[1], 
                    blk->neighid[2], 
                    blk->neighid[3] };

  // output name
  sprintf(blk->output_fname_part,"px%d_py%d", myid2[0],myid2[1]);

  //
  // set boundary itype
  //

  int abs_itype = FD_BOUNDARY_TYPE_NONE;

  for (int i=0; i<FD_NDIM_2; i++) {
    // default none
    int bdry_itype = FD_BOUNDARY_TYPE_NONE;
    blk->abs_num_of_layers[i] = 0;

    // set according input
    if ( (i==0 && myid2[0]==0) ||
         (i==1 && myid2[0]==number_of_mpiprocs_x -1) ||
         (i==2 && myid2[1] == 0) ||
         (i==3 && myid2[1] == number_of_mpiprocs_y -1) ||
         (i>=4))
    {
      // cfs-pml
      if (strcmp(boundary_type_name[i], "cfspml")==0)
      {
        bdry_itype = FD_BOUNDARY_TYPE_CFSPML;
        abs_itype  = bdry_itype;
        blk->abs_num_of_layers[i] = abs_num_of_layers[i];
      }

      // free
      if (strcmp(boundary_type_name[i], "free"  )==0) {
        bdry_itype = FD_BOUNDARY_TYPE_FREE;
      }
    }

    blk->boundary_itype[i] = bdry_itype;
    blk->abs_itype         = abs_itype;
  }

  // check this proc is boundary or not
  //if (neighid[0] == MPI_PROC_NULL) ;

  //
  // set grid size
  //

  // determine ni
  int nx_et = number_of_total_grid_points_x;
  // double cfspml load
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML) {
    nx_et += abs_num_of_layers[0] + abs_num_of_layers[1];
  }
  int nx_avg  = nx_et / number_of_mpiprocs_x;
  int nx_left = nx_et % number_of_mpiprocs_x;
  if (nx_avg < 2 * fdx_nghosts) {
    // error
  }
  if (nx_avg<abs_num_of_layers[0] || nx_avg<abs_num_of_layers[1]) {
    // error
  }
  int ni = nx_avg;
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[0]==0) {
    ni -= abs_num_of_layers[0];
  }
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[0]==number_of_mpiprocs_x-1) {
    ni -= abs_num_of_layers[1];
  }
  // not equal divided points given to first nx_left procs
  if (myid2[0] < nx_left) {
    ni++;
  }
  // global index
  if (myid2[0]==0) {
    blk->gni1 = 0;
  } else {
    blk->gni1 = myid2[0] * nx_avg - abs_num_of_layers[0];
  }
  if (nx_left != 0) {
    blk->gni1 += (myid2[0] < nx_left)? myid2[0] : nx_left;
  }

  // determine nj
  int ny_et = number_of_total_grid_points_y;
  // double cfspml load
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML) {
    ny_et += abs_num_of_layers[2] + abs_num_of_layers[3];
  }
  int ny_avg  = ny_et / number_of_mpiprocs_y;
  int ny_left = ny_et % number_of_mpiprocs_y;
  if (ny_avg < 2 * fdy_nghosts) {
    // error
  }
  if (ny_avg<abs_num_of_layers[2] || ny_avg<abs_num_of_layers[3]) {
    // error
  }
  int nj = ny_avg;
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[1]==0) {
    nj -= abs_num_of_layers[2];
  }
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML && myid2[1]==number_of_mpiprocs_y-1) {
    nj -= abs_num_of_layers[3];
  }
  // not equal divided points given to first ny_left procs
  if (myid2[1] < ny_left) {
    nj++;
  }
  // global index
  if (myid2[1]==0) {
    blk->gnj1 = 0;
  } else {
    blk->gnj1 = myid2[1] * ny_avg - abs_num_of_layers[2];
  }
  if (ny_left != 0) {
    blk->gnj1 += (myid2[1] < ny_left)? myid2[1] : ny_left;
  }

  // determine nk
  int nk = number_of_total_grid_points_z;
  blk->gnk1 = 0;
  
  // add ghost points
  int nx = ni + 2 * fdx_nghosts;
  int ny = nj + 2 * fdy_nghosts;
  int nz = nk + 2 * fdz_nghosts;

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

  // global index end
  blk->gni2 = blk->gni1 + blk->ni - 1;
  blk->gnj2 = blk->gnj1 + blk->nj - 1;
  blk->gnk2 = blk->gnk1 + blk->nk - 1;
  
  // x dimention varies first
  blk->siz_line   = nx; 
  blk->siz_slice  = nx * ny; 
  blk->siz_volume = nx * ny * nz;
  
  // level
  blk->w3d_num_of_levels = number_of_levels;

  // output
  sprintf(blk->output_dir, "%s", output_dir);
  sprintf(blk->grid_export_dir, "%s", grid_export_dir);
  sprintf(blk->media_export_dir, "%s", media_export_dir);

  // alloc pointer
  blk->sta_info = (struct fd_sta_all_t *)malloc(sizeof(struct fd_sta_all_t));

  //fprintf(stdout,"in output_dir=%s\n",output_dir);
  //fprintf(stdout,"blk output_dir=%s\n",blk->output_dir);
  //fflush(stdout);
}

void
fd_blk_set_snapshot(struct fd_blk_t *blk,
                    int  fd_nghosts,
                    int  number_of_snapshot,
                    char **snapshot_name,
                    int *snapshot_index_start,
                    int *snapshot_index_count,
                    int *snapshot_index_incre,
                    int *snapshot_time_start,
                    int *snapshot_time_incre,
                    int *snapshot_save_velocity,
                    int *snapshot_save_stress,
                    int *snapshot_save_strain)
{
  if (number_of_snapshot>0) {
    blk->snap_fname = (char **) fdlib_mem_malloc_2l_char(number_of_snapshot,FD_MAX_STRLEN,"snap_fname");
    blk->snap_info = (int *) malloc(number_of_snapshot * FD_SNAP_INFO_SIZE * sizeof(int));
  }

  // init
  blk->num_of_snap = 0;

  for (int n=0; n < number_of_snapshot; n++)
  {
    int iptr0 = n*FD_NDIM;

    // scan output k-index in this proc
    int gk1 = -1; int ngk =  0; int k_in_nc = 0;
    for (int n3=0; n3<snapshot_index_count[iptr0+2]; n3++)
    {
      int gk = snapshot_index_start[iptr0+2] + n3 * snapshot_index_incre[iptr0+2];
      if (gk >= blk->gnk1 && gk <= blk->gnk2) {
        if (gk1 == -1) {
          gk1 = gk;
          k_in_nc = n3;
        }
        ngk++;
      }
      if (gk > blk->gnk2) break; // no need to larger k
    }

    // scan output j-index in this proc
    int gj1 = -1; int ngj =  0; int j_in_nc = 0;
    for (int n2=0; n2<snapshot_index_count[iptr0+1]; n2++)
    {
      int gj = snapshot_index_start[iptr0+1] + n2 * snapshot_index_incre[iptr0+1];
      if (gj >= blk->gnj1 && gj <= blk->gnj2) {
        if (gj1 == -1) {
          gj1 = gj;
          j_in_nc = n2;
        }
        ngj++;
      }
      if (gj > blk->gnj2) break;
    }

    // scan output i-index in this proc
    int gi1 = -1; int ngi =  0; int i_in_nc = 0;
    for (int n1=0; n1<snapshot_index_count[iptr0+0]; n1++)
    {
      int gi = snapshot_index_start[iptr0+0] + n1 * snapshot_index_incre[iptr0+0];
      if (gi >= blk->gni1 && gi <= blk->gni2) {
        if (gi1 == -1) {
          gi1 = gi;
          i_in_nc = n1;
        }
        ngi++;
      }
      if (gi > blk->gni2) break;
    }

    // if in this proc
    if (ngi>0 && ngj>0 & ngk>0)
    {
      int isnap = blk->num_of_snap;
      int iptr = isnap*FD_SNAP_INFO_SIZE;

      blk->snap_info[iptr + FD_SNAP_INFO_I1]  = gi1 - blk->gni1 + fd_nghosts;
      blk->snap_info[iptr + FD_SNAP_INFO_J1]  = gj1 - blk->gnj1 + fd_nghosts;
      blk->snap_info[iptr + FD_SNAP_INFO_K1]  = gk1 - blk->gnk1 + fd_nghosts;
      blk->snap_info[iptr + FD_SNAP_INFO_NI]  = ngi;
      blk->snap_info[iptr + FD_SNAP_INFO_NJ]  = ngj;
      blk->snap_info[iptr + FD_SNAP_INFO_NK]  = ngk;
      blk->snap_info[iptr + FD_SNAP_INFO_DI]  = snapshot_index_incre[iptr0+0];
      blk->snap_info[iptr + FD_SNAP_INFO_DJ]  = snapshot_index_incre[iptr0+1];
      blk->snap_info[iptr + FD_SNAP_INFO_DK]  = snapshot_index_incre[iptr0+2];

      blk->snap_info[iptr + FD_SNAP_INFO_IT1]  = snapshot_time_start[n];
      blk->snap_info[iptr + FD_SNAP_INFO_DIT]  = snapshot_time_incre[n];

      blk->snap_info[iptr + FD_SNAP_INFO_VEL   ]  = snapshot_save_velocity[n];
      blk->snap_info[iptr + FD_SNAP_INFO_STRESS]  = snapshot_save_stress[n];
      blk->snap_info[iptr + FD_SNAP_INFO_STRAIN]  = snapshot_save_strain[n];

      blk->snap_info[iptr + FD_SNAP_INFO_GI1]  = i_in_nc;
      blk->snap_info[iptr + FD_SNAP_INFO_GJ1]  = j_in_nc;
      blk->snap_info[iptr + FD_SNAP_INFO_GK1]  = k_in_nc;

      sprintf(blk->snap_fname[isnap],"%s/%s_%s.nc",blk->output_dir,
                                                   snapshot_name[n],
                                                   blk->output_fname_part);

      blk->num_of_snap += 1;
    }
  }
}

void fd_blk_malloc_station(struct fd_blk_t *blk, int nt_total)
{
  blk->sta_seismo = (float *) malloc(blk->sta_info->total_number * 
                                     blk->w3d_num_of_vars *
                                     nt_total * sizeof(int));
}

void
fd_blk_set_slice(struct fd_blk_t *blk,
                 int  fd_nghosts,
                 int  number_of_slice_x,
                 int  number_of_slice_y,
                 int  number_of_slice_z,
                 int *slice_x_index,
                 int *slice_y_index,
                 int *slice_z_index)
{
  if (number_of_slice_x>0) {
    blk->slice_x_fname = (char **) fdlib_mem_malloc_2l_char(number_of_slice_x,
                                                            FD_MAX_STRLEN,
                                                            "slice_x_fname");
    blk->slice_x_indx = (int *) malloc(number_of_slice_x * sizeof(int));
  }
  if (number_of_slice_y>0) {
    blk->slice_y_fname = (char **) fdlib_mem_malloc_2l_char(number_of_slice_y,
                                                            FD_MAX_STRLEN,
                                                            "slice_y_fname");
    blk->slice_y_indx = (int *) malloc(number_of_slice_y * sizeof(int));
  }
  if (number_of_slice_z>0) {
    blk->slice_z_fname = (char **) fdlib_mem_malloc_2l_char(number_of_slice_z,
                                                            FD_MAX_STRLEN,
                                                            "slice_z_fname");
    blk->slice_z_indx = (int *) malloc(number_of_slice_z * sizeof(int));
  }

  // init
  blk->num_of_slice_x = 0;
  blk->num_of_slice_y = 0;
  blk->num_of_slice_z = 0;

  // x slice
  for (int n=0; n < number_of_slice_x; n++)
  {
    int gi = slice_x_index[n];
    if (gi >= blk->gni1 && gi <= blk->gni2)
    {
      int islc = blk->num_of_slice_x;

      blk->slice_x_indx[islc]  = gi - blk->gni1 + fd_nghosts;
      sprintf(blk->slice_x_fname[islc],"%s/slicex_i%d_%s.nc", blk->output_dir,gi,blk->output_fname_part);

      blk->num_of_slice_x += 1;
    }
  }

  // y slice
  for (int n=0; n < number_of_slice_y; n++)
  {
    int gj = slice_y_index[n];
    if (gj >= blk->gnj1 && gj <= blk->gnj2)
    {
      int islc = blk->num_of_slice_y;

      blk->slice_y_indx[islc]  = gj - blk->gnj1 + fd_nghosts;
      sprintf(blk->slice_y_fname[islc],"%s/slicey_j%d_%s.nc", blk->output_dir,gj,blk->output_fname_part);

      blk->num_of_slice_y += 1;
    }
  }

  // z slice
  for (int n=0; n < number_of_slice_z; n++)
  {
    int gk = slice_z_index[n];
    if (gk >= blk->gnk1 && gk <= blk->gnk2)
    {
      int islc = blk->num_of_slice_z;

      blk->slice_z_indx[islc]  = gk - blk->gnk1 + fd_nghosts;
      sprintf(blk->slice_z_fname[islc],"%s/slicez_k%d_%s.nc", blk->output_dir,gk,blk->output_fname_part);

      blk->num_of_slice_z += 1;
    }
  }
}

void
fd_blk_locate_inline(struct fd_blk_t *blk,
                     int    fd_nghosts,
                     int    nt_total,
                     int    number_of_receiver_line,
                     int   *receiver_line_index_start,
                     int   *receiver_line_index_incre,
                     int   *receiver_line_count,
                     char **receiver_line_name)
{
  // init
  blk->num_of_point  = 0;

  // first run for count
  for (int n=0; n < number_of_receiver_line; n++)
  {
    for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
    {
      int gi = receiver_line_index_start[n*FD_NDIM+0] + ipt * receiver_line_index_incre[n*FD_NDIM  ];
      int gj = receiver_line_index_start[n*FD_NDIM+1] + ipt * receiver_line_index_incre[n*FD_NDIM+1];
      int gk = receiver_line_index_start[n*FD_NDIM+2] + ipt * receiver_line_index_incre[n*FD_NDIM+2];
      if ( gi >= blk->gni1 && gi <= blk->gni2 &&
           gj >= blk->gnj1 && gj <= blk->gnj2 &&
           gk >= blk->gnk1 && gk <= blk->gnk2 )
      {
        blk->num_of_point += 1;
      }
    }
  }

  // alloc
  if (blk->num_of_point>0)
  {
    blk->point_loc_indx = (int *) malloc(blk->num_of_point * FD_NDIM * sizeof(int));
    blk->point_line_sno = (int *) malloc(blk->num_of_point *           sizeof(int));
    blk->point_line_offset = (int *) malloc(blk->num_of_point *        sizeof(int));
    blk->point_seismo = (float *) fdlib_mem_malloc_1d(
                                    blk->num_of_point * blk->w3d_num_of_vars * nt_total * sizeof(float),
                                    "point_seismo");
    blk->point_coord = (float *) fdlib_mem_malloc_1d(
                                    blk->num_of_point * FD_NDIM * sizeof(float),
                                    "point_seismo");
  }

  // second run for value
  int num_all = 0;
  float *x3d = blk->c3d + blk->c3d_pos[0];
  float *y3d = blk->c3d + blk->c3d_pos[1];
  float *z3d = blk->c3d + blk->c3d_pos[2];
  size_t   siz_line = blk->siz_line;
  size_t   siz_slice = blk->siz_slice;

  for (int n=0; n < number_of_receiver_line; n++)
  {
    for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
    {
      int gi = receiver_line_index_start[n*FD_NDIM+0] + ipt * receiver_line_index_incre[n*FD_NDIM  ];
      int gj = receiver_line_index_start[n*FD_NDIM+1] + ipt * receiver_line_index_incre[n*FD_NDIM+1];
      int gk = receiver_line_index_start[n*FD_NDIM+2] + ipt * receiver_line_index_incre[n*FD_NDIM+2];
      if ( gi >= blk->gni1 && gi <= blk->gni2 &&
           gj >= blk->gnj1 && gj <= blk->gnj2 &&
           gk >= blk->gnk1 && gk <= blk->gnk2 )
      {
        int gi = receiver_line_index_start[n*FD_NDIM+0] + ipt * receiver_line_index_incre[n*FD_NDIM+0];
        int gj = receiver_line_index_start[n*FD_NDIM+1] + ipt * receiver_line_index_incre[n*FD_NDIM+1];
        int gk = receiver_line_index_start[n*FD_NDIM+2] + ipt * receiver_line_index_incre[n*FD_NDIM+2];

        int i = gi - blk->gni1 + fd_nghosts;
        int j = gj - blk->gnj1 + fd_nghosts;
        int k = gk - blk->gnk1 + fd_nghosts;

        int iptr = i + j * siz_line + k * siz_slice;

        blk->point_line_sno[num_all] = n;
        blk->point_line_offset[num_all] = ipt;

        blk->point_loc_indx[num_all*FD_NDIM+0] = i;
        blk->point_loc_indx[num_all*FD_NDIM+1] = j;
        blk->point_loc_indx[num_all*FD_NDIM+2] = k;

        blk->point_coord[num_all*FD_NDIM+0] = x3d[iptr];
        blk->point_coord[num_all*FD_NDIM+1] = y3d[iptr];
        blk->point_coord[num_all*FD_NDIM+2] = z3d[iptr];

        num_all += 1;
      }
    }
  }
}

void
fd_blk_init_mpi_mesg(struct fd_blk_t *blk,
                     int fdx_nghosts,
                     int fdy_nghosts)
{
  int ni = blk->ni;
  int nj = blk->nj;
  int nk = blk->nk;

  // mpi mesg
  int siz_x = (nj * nk * fdx_nghosts) * blk->w3d_num_of_vars;
  int siz_y = (ni * nk * fdy_nghosts) * blk->w3d_num_of_vars;

  int siz_buff =(  ni * nk * fdy_nghosts * 2
                    + nj * nk * fdx_nghosts * 2) * blk->w3d_num_of_vars;
  blk->sbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc sbuff");
  blk->rbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc rbuff");
  for (int iptr=0; iptr < siz_buff; iptr++) {
    blk->rbuff[iptr] = 0.0;
  }

  float *sbuff_x1 = blk->sbuff;
  float *sbuff_x2 = sbuff_x1 + siz_x;
  float *sbuff_y1 = sbuff_x2 + siz_x;
  float *sbuff_y2 = sbuff_y1 + siz_y;

  int tag[4] = { 11, 12, 21, 22 };

  // send
  MPI_Send_init(sbuff_x1, siz_x, MPI_FLOAT, blk->neighid[0], tag[0], blk->topocomm, &(blk->s_reqs[0]));
  MPI_Send_init(sbuff_x2, siz_x, MPI_FLOAT, blk->neighid[1], tag[1], blk->topocomm, &(blk->s_reqs[1]));
  MPI_Send_init(sbuff_y1, siz_y, MPI_FLOAT, blk->neighid[2], tag[2], blk->topocomm, &(blk->s_reqs[2]));
  MPI_Send_init(sbuff_y2, siz_y, MPI_FLOAT, blk->neighid[3], tag[3], blk->topocomm, &(blk->s_reqs[3]));

  float *rbuff_x1 = blk->rbuff;
  float *rbuff_x2 = rbuff_x1 + siz_x;
  float *rbuff_y1 = rbuff_x2 + siz_x;
  float *rbuff_y2 = rbuff_y1 + siz_y;

  // recv
  MPI_Recv_init(rbuff_x1, siz_x, MPI_FLOAT, blk->neighid[0], tag[1], blk->topocomm, &(blk->r_reqs[0]));
  MPI_Recv_init(rbuff_x2, siz_x, MPI_FLOAT, blk->neighid[1], tag[0], blk->topocomm, &(blk->r_reqs[1]));
  MPI_Recv_init(rbuff_y1, siz_y, MPI_FLOAT, blk->neighid[2], tag[3], blk->topocomm, &(blk->r_reqs[2]));
  MPI_Recv_init(rbuff_y2, siz_y, MPI_FLOAT, blk->neighid[3], tag[2], blk->topocomm, &(blk->r_reqs[3]));
}

//- for future: consider different left/right length
//void
//fd_blk_pack_mesg(float *restrict w_cur,
//                 int num_of_vars,
//                 int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
//                 size_t siz_line, size_t siz_slice, size_t siz_volume,
//                 int   *fdx_info,
//                 int   *fdy_info,
//                 int   *fdz_info,
//                 float *restrict buff)
//{
//  int npt_to_x1 = fdx_info[FD_INFO_LENGTH_RIGTH];
//  int npt_to_x2 = fdx_info[FD_INFO_LENGTH_LEFT ];
//
//  size_t iptr_b = 0;
//
//  for (int ivar=0; ivar<num_of_vars; ivar++)
//  {
//    size_t iptr_var = ivar * siz_volume;
//
//    // x1
//    for (int k=nk1; k<=nk2; k++)
//    {
//      size_t iptr_k = iptr_var + k * siz_slice;
//
//      for (int j=nj1; j<=nj2; j++)
//      {
//        size_t iptr_j = iptr_k + j * siz_line;
//
//        for (int i=ni1; i<ni1+npt_to_x1; i++)
//        {
//          buff[iptr_b] = w_cur[iptr_j + i];
//          iptr_b++;
//        }
//      }
//    } 
//
//    // y1y2
//
//  } // ivar
//}

void
fd_blk_pack_mesg(float *restrict w_cur,float *restrict sbuff,
                 int num_of_vars,
                 int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
                 size_t siz_line, size_t siz_slice, size_t siz_volume,
                 int   fdx_nghosts,
                 int   fdy_nghosts)
{
  size_t iptr_b = 0;

  // x1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<ni1+fdx_nghosts; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // x2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni2-fdx_nghosts+1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // y1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<nj1+fdy_nghosts; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // y2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj2-fdy_nghosts+1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  } // ivar

}

void
fd_blk_unpack_mesg(float *restrict rbuff,float *restrict w_cur,
                 int num_of_vars,
                 int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
                 int nx, int ny,
                 size_t siz_line, size_t siz_slice, size_t siz_volume)
{
  size_t iptr_b = 0;

  // x1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=0; i<ni1; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  }

  // x2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni2+1; i<nx; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  }

  // y1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=0; j<nj1; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  }

  // y2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj2+1; j<ny; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  } // ivar
}

void
fd_print(struct fd_t *fd)
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

  fprintf(stdout, " fd_len = %d\n", fd->fd_len);
  fprintf(stdout, " fd_half_len = %d\n", fd->fd_half_len);
  fprintf(stdout, " fd_nghosts = %d\n", fd->fd_nghosts);
  fprintf(stdout, "  fd_indx = ");
  for (int i=0; i<fd->fd_len; i++) {
    fprintf(stdout, " %d", fd->fd_indx[i]);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "  fd_coef = ");
  for (int i=0; i<fd->fd_len; i++) {
    fprintf(stdout, " %e", fd->fd_coef[i]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, " fdx_max_len = %d\n", fd->fdx_max_len);
  // rest

}

void
fd_blk_print(struct fd_blk_t *blk)
{    

  fprintf(stdout, "\n-------------------------------------------------------\n");
  fprintf(stdout, "print blk %s:\n", blk->name);
  fprintf(stdout, "-------------------------------------------------------\n\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> ESTIMATE MEMORY INFO.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "total memory size Byte: %20.5f  B\n", PSV->total_memory_size_Byte);
  //fprintf(stdout, "total memory size KB  : %20.5f KB\n", PSV->total_memory_size_KB  );
  //fprintf(stdout, "total memory size MB  : %20.5f MB\n", PSV->total_memory_size_MB  );
  //fprintf(stdout, "total memory size GB  : %20.5f GB\n", PSV->total_memory_size_GB  );
  //fprintf(stdout, "\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> FOLDER AND FILE INFO.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "   OutFolderName: %s\n", OutFolderName);
  //fprintf(stdout, "       EventName: %s\n", OutPrefix);
  //fprintf(stdout, "     LogFilename: %s\n", LogFilename);
  //fprintf(stdout, " StationFilename: %s\n", StationFilename);
  //fprintf(stdout, "  SourceFilename: %s\n", SourceFilename);
  //fprintf(stdout, "   MediaFilename: %s\n", MediaFilename);
  //fprintf(stdout, "\n");

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " nx    = %-10d\n", blk->nx);
  fprintf(stdout, " ny    = %-10d\n", blk->ny);
  fprintf(stdout, " nz    = %-10d\n", blk->nz);
  fprintf(stdout, " ni    = %-10d\n", blk->ni);
  fprintf(stdout, " nj    = %-10d\n", blk->nj);
  fprintf(stdout, " nk    = %-10d\n", blk->nk);

  fprintf(stdout, " ni1   = %-10d\n", blk->ni1);
  fprintf(stdout, " ni2   = %-10d\n", blk->ni2);
  fprintf(stdout, " nj1   = %-10d\n", blk->nj1);
  fprintf(stdout, " nj2   = %-10d\n", blk->nj2);
  fprintf(stdout, " nk1   = %-10d\n", blk->nk1);
  fprintf(stdout, " nk2   = %-10d\n", blk->nk2);

  fprintf(stdout, " gni1   = %-10d\n", blk->gni1);
  fprintf(stdout, " gni2   = %-10d\n", blk->gni2);
  fprintf(stdout, " gnj1   = %-10d\n", blk->gnj1);
  fprintf(stdout, " gnj2   = %-10d\n", blk->gnj2);
  fprintf(stdout, " gnk1   = %-10d\n", blk->gnk1);
  fprintf(stdout, " gnk2   = %-10d\n", blk->gnk2);

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> media info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //if (blk->media_type == MEDIA_TYPE_LAYER)
  //{
  //    strcpy(str, "layer");
  //}
  //else if (blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    strcpy(str, "grid");
  //}
  //fprintf(stdout, " media_type = %s\n", str);
  //if(blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    fprintf(stdout, "\n --> the media filename is:\n");
  //    fprintf(stdout, " velp_file  = %s\n", blk->fnm_velp);
  //    fprintf(stdout, " vels_file  = %s\n", blk->fnm_vels);
  //    fprintf(stdout, " rho_file   = %s\n", blk->fnm_rho);
  //}
  //fprintf(stdout, "\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> source info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, " number_of_force  = %d\n", blk->number_of_force);
  //if(blk->number_of_force > 0)
  //{
  //    fprintf(stdout, " force_source           x           z     x_shift     z_shift           i           k:\n");
  //    for(n=0; n<blk->number_of_force; n++)
  //    {
  //        indx = 2*n;
  //        fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //                blk->force_coord[indx], blk->force_coord[indx+1],
  //                blk->force_shift[indx], blk->force_shift[indx+1],
  //                blk->force_indx [indx], blk->force_indx [indx+1]);
  //    }
  //    fprintf(stdout, "\n");
  //}

  //fprintf(stdout, "\n");
  //fprintf(stdout, " number_of_moment = %d\n", blk->number_of_moment);
  //if(blk->number_of_moment > 0)
  //{
  //    fprintf(stdout, " moment_source          x           z     x_shift     z_shift           i           k:\n");
  //    for(n=0; n<blk->number_of_moment; n++)
  //    {
  //        indx = 2*n;
  //        fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //                blk->moment_coord[indx], blk->moment_coord[indx+1],
  //                blk->moment_shift[indx], blk->moment_shift[indx+1],
  //                blk->moment_indx [indx], blk->moment_indx [indx+1]);
  //    }
  //    fprintf(stdout, "\n");
  //}

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> boundary layer information:\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //ierr = boundary_id2type(type1, blk->boundary_type[0], errorMsg);
  //ierr = boundary_id2type(type2, blk->boundary_type[1], errorMsg);
  //ierr = boundary_id2type(type3, blk->boundary_type[2], errorMsg);
  //ierr = boundary_id2type(type4, blk->boundary_type[3], errorMsg);
  //fprintf(stdout, " boundary_type         = %10s%10s%10s%10s\n", 
  //        type1, type2, type3, type4);
  //fprintf(stdout, " boundary_layer_number = %10d%10d%10d%10d\n", 
  //        blk->boundary_layer_number[0], blk->boundary_layer_number[1], 
  //        blk->boundary_layer_number[2], blk->boundary_layer_number[3]);
  //fprintf(stdout, "\n");
  //fprintf(stdout, " absorb_velocity       = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->absorb_velocity[0], blk->absorb_velocity[1], blk->absorb_velocity[2], 
  //        blk->absorb_velocity[3]);
  //fprintf(stdout, "\n");
  //fprintf(stdout, " CFS_alpha_max         = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->CFS_alpha_max[0], blk->CFS_alpha_max[1], blk->CFS_alpha_max[2], 
  //        blk->CFS_alpha_max[3]);
  //fprintf(stdout, " CFS_beta_max          = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->CFS_beta_max[0], blk->CFS_beta_max[1], blk->CFS_beta_max[2], 
  //        blk->CFS_beta_max[3]);
  
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output information:\n");
  fprintf(stdout, "-------------------------------------------------------\n");

  fprintf(stdout, "--> num_of_slice_x = %d\n", blk->num_of_slice_x);
  for (int n=0; n<blk->num_of_slice_x; n++)
  {
    fprintf(stdout, "  #%d, i=%d, fname=%s\n", n, blk->slice_x_indx[n],blk->slice_x_fname[n]);
  }
  fprintf(stdout, "--> num_of_slice_y = %d\n", blk->num_of_slice_y);
  for (int n=0; n<blk->num_of_slice_y; n++)
  {
    fprintf(stdout, "  #%d, j=%d, fname=%s\n", n, blk->slice_y_indx[n],blk->slice_y_fname[n]);
  }
  fprintf(stdout, "--> num_of_slice_z = %d\n", blk->num_of_slice_z);
  for (int n=0; n<blk->num_of_slice_z; n++)
  {
    fprintf(stdout, "  #%d, k=%d, fname=%s\n", n, blk->slice_z_indx[n],blk->slice_z_fname[n]);
  }

  fprintf(stdout, "--> num_of_snap = %d\n", blk->num_of_snap);
  fprintf(stdout, "#   i0 j0 k0 ni nj nk di dj dk it0 dit vel stress strain gi1 gj1 gk1\n");
  for (int n=0; n<blk->num_of_snap; n++)
  {
    int iptr = n * FD_SNAP_INFO_SIZE;

    fprintf(stdout, " %d", n);
    for (int i=0; i<FD_SNAP_INFO_SIZE; i++)
    {
      fprintf(stdout, " %d", blk->snap_info[iptr+i]);
    }
    fprintf(stdout, "\n");
  }

  //fprintf(stdout, "\n");
  //fprintf(stdout, "--> station information.\n");
  //fprintf(stdout, " number_of_station  = %4d\n", blk->number_of_station);
  //fprintf(stdout, " seismo_format_sac  = %4d\n", blk->seismo_format_sac );
  //fprintf(stdout, " seismo_format_segy = %4d\n", blk->seismo_format_segy);
  //fprintf(stdout, " SeismoPrefix = %s\n", SeismoPrefix);
  //fprintf(stdout, "\n");

  //if(blk->number_of_station > 0)
  //{
  //    //fprintf(stdout, " station_indx:\n");
  //    fprintf(stdout, " stations             x           z           i           k:\n");
  //}

  //for(n=0; n<blk->number_of_station; n++)
  //{
  //    indx = 2*n;
  //    fprintf(stdout, "       %04d  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //            blk->station_coord[indx], blk->station_coord[indx+1],
  //            blk->station_indx [indx], blk->station_indx [indx+1]);
  //}
  //fprintf(stdout, "\n");

  return;
}
