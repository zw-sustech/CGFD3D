#ifndef FD_T_H
#define FD_T_H

#include <mpi.h>

// consts
#define FD_NDIM   3
#define FD_NDIM_2 6 // 2 * ndim
#define FD_MAX_STRLEN 1000

// for fd_info array
#define FD_INFO_POS_OF_INDX  0
#define FD_INFO_LENGTH_TOTAL 1
#define FD_INFO_LENGTH_HALF  2
#define FD_INFO_LENGTH_LEFT  3
#define FD_INFO_LENGTH_RIGTH 4
#define FD_INFO_SIZE         5

// for boundary_ityp
#define FD_BOUNDARY_TYPE_NONE   0
#define FD_BOUNDARY_TYPE_FREE   1
#define FD_BOUNDARY_TYPE_CFSPML 2
#define FD_BOUNDARY_TYPE_ABLEXP 3
#define FD_BOUNDARY_TYPE_MPML   4
#define FD_BOUNDARY_TYPE_DPML   5

// slice info
#define FD_SLICE_INFO_ID   0
#define FD_SLICE_INFO_INDX 0

// snap info
#define FD_SNAP_INFO_I1 0
#define FD_SNAP_INFO_J1 1
#define FD_SNAP_INFO_K1 2
#define FD_SNAP_INFO_NI 3
#define FD_SNAP_INFO_NJ 4
#define FD_SNAP_INFO_NK 5
#define FD_SNAP_INFO_DI 6
#define FD_SNAP_INFO_DJ 7
#define FD_SNAP_INFO_DK 8
#define FD_SNAP_INFO_IT1 9
#define FD_SNAP_INFO_DIT 10
#define FD_SNAP_INFO_VEL 11
#define FD_SNAP_INFO_STRESS 12
#define FD_SNAP_INFO_STRAIN 13
#define FD_SNAP_INFO_GI1 14
#define FD_SNAP_INFO_GJ1 15
#define FD_SNAP_INFO_GK1 16
#define FD_SNAP_INFO_SIZE 17

// macro for fd opterators

// use siz_shift to find adjacent point of the stentil for 3d var
#define M_FD_SHIFT(deriv, var, iptr, fd_length, fd_shift, fd_coef, n) \
   deriv = 0.0; \
   for (n=0; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n]]; \
   }

#define M_FD_SHIFT_PTR(deriv, var_ptr, fd_length, fd_shift, fd_coef, n) \
  deriv = 0.0;                                                        \
  for (n = 0; n < fd_length; n++)                                     \
  {                                                                   \
    deriv += fd_coef[n] * *(var_ptr + fd_shift[n]);                    \
  }

#define M_FD_SHIFT_PTR_UNLOOP5(deriv, var_ptr, fd_length, fd_shift, fd_coef, n) \
  deriv =  fd_coef[0] * *(var_ptr + fd_shift[0])                    \
          +fd_coef[1] * *(var_ptr + fd_shift[1])                    \
          +fd_coef[2] * *(var_ptr + fd_shift[2])                    \
          +fd_coef[3] * *(var_ptr + fd_shift[3])                    \
          +fd_coef[4] * *(var_ptr + fd_shift[4]);

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

/*******************************************************************************
 * structure for different fd schemes
 ******************************************************************************/

struct fd_t{

  //----------------------------------------------------------------------------
  // Runge-Kutta time scheme
  //----------------------------------------------------------------------------

  int num_rk_stages;
  float *rk_a;
  float *rk_b;

  //----------------------------------------------------------------------------
  // single centered scheme for all dim
  //----------------------------------------------------------------------------

  int     fd_len;
  int     fd_half_len;
  int     fd_nghosts;
  int    *fd_indx;
  float  *fd_coef;

  //----------------------------------------------------------------------------
  // para for different schemes at points to boundaries for different dim
  //----------------------------------------------------------------------------

  // max total len of op
  int fdx_max_len;
  int fdy_max_len;
  int fdz_max_len;

  // max half len
  int fdx_max_half_len;
  int fdy_max_half_len;
  int fdz_max_half_len;

  // number of layers that need to use biased op near boundary
  int fdx_num_surf_lay;
  int fdy_num_surf_lay;
  int fdz_num_surf_lay;

  // ghost point required 
  int fdx_nghosts;
  int fdy_nghosts;
  int fdz_nghosts;

  //----------------------------------------------------------------------------
  // center schemes at different points to boundary
  //  fd_ is 2d pointer, the first pointer means the grid layer to free surface (from 0 to fd_len-1),
  //  the second pointer points to fd op for that layer, the size could be different, larger than
  //  inner op
  //----------------------------------------------------------------------------

  int    **fdx_all_info; // [k2free][pos, total, half, left, right]
  int    **fdx_all_indx;
  float  **fdx_all_coef;

  int    **fdy_all_info;
  int    **fdy_all_indx;
  float  **fdy_all_coef;

  int    **fdz_all_info;
  int    **fdz_all_indx;
  float  **fdz_all_coef;

  //----------------------------------------------------------------------------
  // filter schemes at different points to boundary
  //----------------------------------------------------------------------------

  int    **filtx_all_info; // [k2free][pos, total, half, left, right] 
  int    **filtx_all_indx;
  float  **filtx_all_coef;

  int    **filty_all_info;
  int    **filty_all_indx;
  float  **filty_all_coef;

  int    **filtz_all_info;
  int    **filtz_all_indx;
  float  **filtz_all_coef;

  //----------------------------------------------------------------------------
  // pairs for 3d space for MacCormack-type schemes
  //----------------------------------------------------------------------------

  int num_of_pairs;

  int    ****pair_fdx_all_info;  // [pair][stage][k2free][pos, total, half, left, right]
  int     ***pair_fdx_all_indx;  // [pair][stage][1-10 ele], not include [hlen] due to different length
  float   ***pair_fdx_all_coef;

  int    ****pair_fdy_all_info; 
  int     ***pair_fdy_all_indx; 
  float   ***pair_fdy_all_coef;

  int    ****pair_fdz_all_info;
  int     ***pair_fdz_all_indx;
  float   ***pair_fdz_all_coef;
};

/*******************************************************************************
 * block structure
 ******************************************************************************/

struct fd_blk_t
{
  // name for output file name
  char name[FD_MAX_STRLEN];

  //
  // grid index
  //
  int nx, ny, nz;
  int ni1, ni2, nj1, nj2, nk1, nk2, ni, nj, nk;
  int gni1, gnj1, gnk1; // global index
  int gni2, gnj2, gnk2; // global index
  
  // size of a single var
  size_t siz_line;
  size_t siz_slice;
  size_t siz_volume; // number of points per var
  
  //size_t siz_vars; // volume * num_of_vars, not easy for understand, may named with w3d and aux

  // mapping coord name
  char **coord_name;
  
  // coordnate: x3d, y3d, z3d
  float  *c3d; // grid 3d vars
  int     c3d_num_of_vars;
  size_t *c3d_pos; // postion 0 for each var
  char  **c3d_name;
  
  // grid metrics: jac, xi_x, etc
  float  *g3d; // grid 3d vars
  int     g3d_num_of_vars;
  size_t *g3d_pos; // postion 0 for each var
  char  **g3d_name;
  
  //
  // media: rho, lambda, mu etc
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
  //size_t  w3d_size_per_level; // can be cal by num_of_vars * siz_volume on fly
  size_t *w3d_pos; // pos of 0 elem for each vars in one level
  char  **w3d_name;
  
  // auxiliary vars, may be used in the future
  //
  float  *a3d; // auxiliary vars
  int     a3d_num_of_vars;
  int     a3d_num_of_levels;
  //size_t  a3d_size_per_level;
  size_t *a3d_pos;
  char  **a3d_name;
  
  // boundary 
  int boundary_itype[ FD_NDIM_2 ];
  
  // free surface
  float *matVx2Vz, *matVy2Vz;
  
  // source term
  int     num_of_force;
  int    *force_info; // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
  float  *force_vec_stf;
  int     num_of_moment;
  int    *moment_info; // num_of_force * 7 : si,sj,sk,start_pos_in_rate,start_it, end_it, n_ext
  float  *moment_ten_rate; // stage, it, Mij, num
  int    *moment_ext_indx;
  float  *moment_ext_coef;
  
  //
  // abs
  //
  
  // only support one type each run, may combine in the future
  //  which requires different vars in blk_t
  int  abs_itype;
  
  // only for this block, may diff with global values from input par
  int  abs_num_of_layers[ FD_NDIM_2 ];
  
  // grid index of each face
  int    abs_indx[FD_NDIM_2 * FD_NDIM_2];
  
  int      abs_coefs_facepos0[FD_NDIM_2];  // 
  float   *abs_coefs; // all coefs all faces 
  
  // may add abs_numb_of_vars later for mpml
  //   may not need it since 
  
  int     abs_vars_volsiz[FD_NDIM_2]; // size of single var in abs_blk_vars for each block
  int     abs_vars_facepos0[FD_NDIM_2]; // start pos of each blk in first level; other level similar
  int     abs_vars_size_per_level; // size of vars per level
  float   *abs_vars; //  order: vars_one_block -> all block -> one level -> 4 levels
  
  // io
  int     num_of_sta;
  int    *sta_loc_point;
  float  *sta_loc_dxyz;
  float  *sta_seismo;

  int     num_of_point;  // for saving in solver
  int    *point_loc_point;
  float  *point_seismo;
  int     num_of_inline; // for writing to file
  char  **inline_fname;
  int    *inline_pos;
  int    *inline_num;

  int     num_of_slice_x;
  int*    slice_x_indx;
  char  **slice_x_fname;

  int     num_of_slice_y;
  int*    slice_y_indx;
  char  **slice_y_fname;

  int     num_of_slice_z;
  int*    slice_z_indx;
  char  **slice_z_fname;
  
  int     num_of_snap;
  int    *snap_info;
  char  **snap_fname;
  
  // fname and dir
  char output_fname_part[FD_MAX_STRLEN];
  // wavefield output
  char output_dir[FD_MAX_STRLEN];
  // seperate grid output to save grid for repeat simulation
  char grid_export_dir[FD_MAX_STRLEN];
  // seperate medium output to save medium for repeat simulation
  char media_export_dir[FD_MAX_STRLEN];

  // mpi mesg
  int    myid2[2];
  int    neighid[4];
  float *sbuff;
  float *rbuff;
  MPI_Comm    topocomm;
  MPI_Request r_reqs[4];
  MPI_Request s_reqs[4];

  // exchange between blocks
  // point-to-point values
  int num_of_conn_points;
  int *conn_this_indx;
  int *conn_out_blk;
  int *conn_out_indx;
  // interp point
  int num_of_interp_points;
  int *interp_this_indx;
  int *interp_out_blk;
  int *interp_out_indx;
  float *interp_out_dxyz;

  // mem usage
  size_t number_of_float;
  size_t number_of_btye;
};

/*******************************************************************************
 * function prototype
 ******************************************************************************/

void 
fd_set_macdrp(struct fd_t *fd);

void
fd_blk_init(struct fd_blk_t *blk,
            int number_of_total_grid_points_x,
            int number_of_total_grid_points_y,
            int number_of_total_grid_points_z,
            int number_of_mpiprocs_x,
            int number_of_mpiprocs_y,
            char **boundary_type_name,
            int *abs_number_of_layers,
            char *output_dir,
            char *grid_export_dir,
            char *media_export_dir,
            int fdx_nghosts,
            int fdy_nghosts,
            int fdz_nghosts,
            int number_of_levels, // depends on time scheme, for rk4 = 4
            MPI_Comm comm, 
            const int myid, const int verbose);

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
                    int *snapshot_save_strain);

void
fd_blk_set_slice(struct fd_blk_t *blk,
                 int  fd_nghosts,
                 int  number_of_slice_x,
                 int  number_of_slice_y,
                 int  number_of_slice_z,
                 int *slice_x_index,
                 int *slice_y_index,
                 int *slice_z_index);

void
fd_blk_set_inline(struct fd_blk_t *blk,
                  int  fd_nghosts,
                  int  nt_total,
                  int  number_of_receiver_line,
                  int *receiver_line_index_start,
                  int *receiver_line_index_incre,
                  int *receiver_line_count,
                  char **receiver_line_name);

void
fd_blk_set_sta(struct fd_blk_t *blk);

void
fd_blk_init_mpi_mesg(struct fd_blk_t *blk,
                     int fdx_nghosts,
                     int fdy_nghosts);

void
fd_blk_pack_mesg(float *restrict w_cur,float *restrict sbuff,
                 int num_of_vars,
                 int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
                 size_t siz_line, size_t siz_slice, size_t siz_volume,
                 int   fdx_nghosts,
                 int   fdy_nghosts);

void
fd_blk_unpack_mesg(float *restrict rbuff,float *restrict w_cur,
                 int num_of_vars,
                 int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
                 int nx, int ny,
                 size_t siz_line, size_t siz_slice, size_t siz_volume);

void
fd_print(struct fd_t *fd);

void
fd_blk_print(struct fd_blk_t *blk);

#endif
