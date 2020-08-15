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

  fd->rk_a[0] = 0.5;
  fd->rk_a[1] = 0.5;
  fd->rk_a[2] = 1.0;

  fd->rk_b[0] = 1.0/6.0;
  fd->rk_b[1] = 1.0/3.0;
  fd->rk_b[2] = 1.0/3.0;
  fd->rk_b[3] = 1.0/6.0;

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
      1.0, -1.0,
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
            char *name,
            int number_of_total_grid_points_x,
            int number_of_total_grid_points_y,
            int number_of_total_grid_points_z,
            int number_of_mpiprocs_x,
            int number_of_mpiprocs_y,
            char **boundary_type_name,
            int *abs_num_of_layers,
            char *output_dir,
            int fdx_nghosts,
            int fdy_nghosts,
            int fdz_nghosts,
            int number_of_levels, // depends on time scheme, for rk4 = 4
            int *myid2,
            int *neighid,
            const int myid, const int verbose)
{
  int boundary_itype[FD_NDIM_2];

  // set name
  sprintf(blk->name, "%s", name);

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
  
  // x dimention varies first
  blk->siz_line   = nx; 
  blk->siz_slice  = nx * ny; 
  blk->siz_volume = nx * ny * nz;
  
  // level
  blk->w3d_num_of_levels = number_of_levels;

  // output
  sprintf(blk->output_dir, "%s", output_dir);
}

/*
void
fd_blk_print(struct fd_blk_t *blk)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
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
  fprintf(stdout, "--> GRID INFO.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " ni    = %-10d\n", PSV->ni);
  fprintf(stdout, " nk    = %-10d\n", PSV->nk);
  fprintf(stdout, " nx    = %-10d\n", PSV->nx);
  fprintf(stdout, " nz    = %-10d\n", PSV->nz);
  fprintf(stdout, " ni1   = %-10d\n", PSV->ni1);
  fprintf(stdout, " ni2   = %-10d\n", PSV->ni2);
  fprintf(stdout, " nk1   = %-10d\n", PSV->nk1);
  fprintf(stdout, " nk2   = %-10d\n", PSV->nk2);
  fprintf(stdout, " LenFD = %-10d\n", PSV->LenFD);
  fprintf(stdout, " nt    = %-10d\n", PSV->nt);
  fprintf(stdout, " stept = %10.4e\n", PSV->stept);
  fprintf(stdout, " maxdt = %10.4e\n", PSV->maxdt);
  fprintf(stdout, " steph = %10.4e\n", PSV->steph);
  fprintf(stdout, " x0    = %10.4e\n", PSV->x0);
  fprintf(stdout, " x1    = %10.4e\n", PSV->x1);
  fprintf(stdout, " z0    = %10.4e\n", PSV->z0);
  fprintf(stdout, " z1    = %10.4e\n", PSV->z1);
  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> media info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  if (PSV->media_type == MEDIA_TYPE_LAYER)
  {
      strcpy(str, "layer");
  }
  else if (PSV->media_type == MEDIA_TYPE_GRID)
  {
      strcpy(str, "grid");
  }
  fprintf(stdout, " media_type = %s\n", str);
  if(PSV->media_type == MEDIA_TYPE_GRID)
  {
      fprintf(stdout, "\n --> the media filename is:\n");
      fprintf(stdout, " velp_file  = %s\n", PSV->fnm_velp);
      fprintf(stdout, " vels_file  = %s\n", PSV->fnm_vels);
      fprintf(stdout, " rho_file   = %s\n", PSV->fnm_rho);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> source info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " number_of_force  = %d\n", PSV->number_of_force);
  if(PSV->number_of_force > 0)
  {
      fprintf(stdout, " force_source           x           z     x_shift     z_shift           i           k:\n");
      for(n=0; n<PSV->number_of_force; n++)
      {
          indx = 2*n;
          fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
                  PSV->force_coord[indx], PSV->force_coord[indx+1],
                  PSV->force_shift[indx], PSV->force_shift[indx+1],
                  PSV->force_indx [indx], PSV->force_indx [indx+1]);
      }
      fprintf(stdout, "\n");
  }

  fprintf(stdout, "\n");
  fprintf(stdout, " number_of_moment = %d\n", PSV->number_of_moment);
  if(PSV->number_of_moment > 0)
  {
      fprintf(stdout, " moment_source          x           z     x_shift     z_shift           i           k:\n");
      for(n=0; n<PSV->number_of_moment; n++)
      {
          indx = 2*n;
          fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
                  PSV->moment_coord[indx], PSV->moment_coord[indx+1],
                  PSV->moment_shift[indx], PSV->moment_shift[indx+1],
                  PSV->moment_indx [indx], PSV->moment_indx [indx+1]);
      }
      fprintf(stdout, "\n");
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> boundary layer information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  ierr = boundary_id2type(type1, PSV->boundary_type[0], errorMsg);
  ierr = boundary_id2type(type2, PSV->boundary_type[1], errorMsg);
  ierr = boundary_id2type(type3, PSV->boundary_type[2], errorMsg);
  ierr = boundary_id2type(type4, PSV->boundary_type[3], errorMsg);
  fprintf(stdout, " boundary_type         = %10s%10s%10s%10s\n", 
          type1, type2, type3, type4);
  fprintf(stdout, " boundary_layer_number = %10d%10d%10d%10d\n", 
          PSV->boundary_layer_number[0], PSV->boundary_layer_number[1], 
          PSV->boundary_layer_number[2], PSV->boundary_layer_number[3]);
  fprintf(stdout, "\n");
  fprintf(stdout, " absorb_velocity       = %10.2f%10.2f%10.2f%10.2f\n", 
          PSV->absorb_velocity[0], PSV->absorb_velocity[1], PSV->absorb_velocity[2], 
          PSV->absorb_velocity[3]);
  fprintf(stdout, "\n");
  fprintf(stdout, " CFS_alpha_max         = %10.2f%10.2f%10.2f%10.2f\n", 
          PSV->CFS_alpha_max[0], PSV->CFS_alpha_max[1], PSV->CFS_alpha_max[2], 
          PSV->CFS_alpha_max[3]);
  fprintf(stdout, " CFS_beta_max          = %10.2f%10.2f%10.2f%10.2f\n", 
          PSV->CFS_beta_max[0], PSV->CFS_beta_max[1], PSV->CFS_beta_max[2], 
          PSV->CFS_beta_max[3]);
  
  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> snapshot information.\n");
  if (PSV->number_of_snapshot > 0)
  {
      fprintf(stdout, "number_of_snapshot=%d\n", PSV->number_of_snapshot);
      fprintf(stdout, "#   x0    z0    nx    nz    dx    dz    dt     tdim_max    component\n");
      for(n=0; n<PSV->number_of_snapshot; n++)
      {
          indx = 10*n;
          componentV = ' ';
          componentT = ' ';

          if(PSV->snapshot_information[indx+8] == 1) {
              componentV = 'V';
          }

          if(PSV->snapshot_information[indx+9] == 1) {
              componentT = 'T';
          }

          for(i=0; i<7; i++)
          {
              fprintf(stdout, "%6d", PSV->snapshot_information[indx+i]);
          }

          fprintf(stdout, "%12d", PSV->snapshot_information[indx+7]);
          fprintf(stdout, "         %c%c\n", componentV, componentT);
      }
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "--> station information.\n");
  fprintf(stdout, " number_of_station  = %4d\n", PSV->number_of_station);
  fprintf(stdout, " seismo_format_sac  = %4d\n", PSV->seismo_format_sac );
  fprintf(stdout, " seismo_format_segy = %4d\n", PSV->seismo_format_segy);
  fprintf(stdout, " SeismoPrefix = %s\n", SeismoPrefix);
  fprintf(stdout, "\n");

  if(PSV->number_of_station > 0)
  {
      //fprintf(stdout, " station_indx:\n");
      fprintf(stdout, " stations             x           z           i           k:\n");
  }

  for(n=0; n<PSV->number_of_station; n++)
  {
      indx = 2*n;
      fprintf(stdout, "       %04d  %10.4e  %10.4e  %10d  %10d\n", n+1, 
              PSV->station_coord[indx], PSV->station_coord[indx+1],
              PSV->station_indx [indx], PSV->station_indx [indx+1]);
  }
  fprintf(stdout, "\n");

  return;
}
*/

void
fd_mpi_create_topo(struct fd_mpi_t *fdmpi, int myid, MPI_Comm comm, int nprocx, int nprocy)
{
  int pdims[2]={nprocx,nprocy};
  int periods[2] = {0,0};

  // create Cartesian topology
  MPI_Cart_create(comm, 2, pdims, periods, 0, &fdmpi->topocomm);

  // get my local x,y coordinates
  MPI_Cart_coords(fdmpi->topocomm, myid, 2, fdmpi->myid2);

  // neighour
  MPI_Cart_shift(fdmpi->topocomm, 0, 1, &(fdmpi->neighid[0]), &(fdmpi->neighid[1]));
  MPI_Cart_shift(fdmpi->topocomm, 1, 1, &(fdmpi->neighid[2]), &(fdmpi->neighid[3]));
}

void
fd_blk_set_snapshot(struct fd_blk_t *blk,
                    int  fd_nghosts,
                    int  number_of_snapshot,
                    int *snapshot_index_start,
                    int *snapshot_index_count,
                    int *snapshot_index_stride,
                    int *snapshot_time_start,
                    int *snapshot_time_count,
                    int *snapshot_time_stride)
{
  blk->num_of_snap = number_of_snapshot;

  blk->snap_grid_indx = (int *) malloc(number_of_snapshot * FD_NDIM*3 * sizeof(int));
  blk->snap_time_indx = (int *) malloc(number_of_snapshot * 3 * sizeof(int));

  for (int n=0; n < number_of_snapshot; n++)
  {
    // start
    blk->snap_grid_indx[n*FD_NDIM*3+0] = snapshot_index_start[n*FD_NDIM+0] + fd_nghosts;
    blk->snap_grid_indx[n*FD_NDIM*3+1] = snapshot_index_start[n*FD_NDIM+1] + fd_nghosts;
    blk->snap_grid_indx[n*FD_NDIM*3+2] = snapshot_index_start[n*FD_NDIM+2] + fd_nghosts;

    // count, convert -1 to num
    for (int i=0; i < FD_NDIM; i++)
    {
      blk->snap_grid_indx[n*FD_NDIM*3+3+i] = snapshot_index_count[n*FD_NDIM+i];
      //if (snapshot_index_count[n*FD_NDIM+i]==-1) {
      //  blk->snap_grid_indx[n*FD_NDIM*3+3] = blk->ni
      //}
    }

    // stride
    for (int i=0; i < FD_NDIM; i++)
    {
      blk->snap_grid_indx[n*FD_NDIM*3+6+i] = snapshot_index_stride[n*FD_NDIM+i];
    }
    
    // time
    blk->snap_time_indx[n*3+0] = snapshot_time_start[n];
    blk->snap_time_indx[n*3+1] = snapshot_time_count[n];
    blk->snap_time_indx[n*3+2] = snapshot_time_stride[n];
  }
}
