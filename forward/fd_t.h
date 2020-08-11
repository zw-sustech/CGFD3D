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

// macro for fd opterators

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

/*******************************************************************************
 * structure for different fd schemes
 ******************************************************************************/

struct fd_t{

  //----------------------------------------------------------------------------
  // Runge-Kutta time scheme
  //----------------------------------------------------------------------------
  int num_rk_stages = 4;
  float rk_a[3] = { 0.5, 0.5, 1.0 };
  float rk_b[4] = { 1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0 };

  //----------------------------------------------------------------------------
	// common para for different schemes
  //----------------------------------------------------------------------------

  int fdx_max_len, fdx_max_half_len, fdx_num_surf_lay;
  int fdy_max_len, fdy_max_half_len, fdy_num_surf_lay;
  int fdz_max_len, fdz_max_half_len, fdz_num_surf_lay;
                // number of layers that need to use biased op near boundary
  int fdx_nghosts; // ghost point required for x-dim
  int fdy_nghosts;
  int fdz_nghosts;

  //----------------------------------------------------------------------------
  // center schemes at different points to boundary
  //  fd_ is 2d pointer, the first pointer means the grid layer to free surface (from 0 to fd_len-1),
  //  the second pointer points to fd op for that layer, the size could be different, larger than
  //  inner op
  //----------------------------------------------------------------------------
  int    **fdx_all_info; // [k2free][pos, total, half, left, right]
  size_t **fdx_all_indx;
  float  **fdx_all_coef;

  int    **fdy_all_info;
  size_t **fdy_all_indx;
  float  **fdy_all_coef;

  int    **fdz_all_info;
  size_t **fdz_all_indx;
  float  **fdz_all_coef;

  //----------------------------------------------------------------------------
  // filter schemes at different points to boundary
  //----------------------------------------------------------------------------
  int    **filtx_all_info; // [k2free][pos, total, half, left, right] 
  size_t **filtx_all_indx;
  float  **filtx_all_coef;

  int    **filty_all_info;
  size_t **filty_all_indx;
  float  **filty_all_coef;

  int    **filtz_all_info;
  size_t **filtz_all_indx;
  float  **filtz_all_coef;

  //----------------------------------------------------------------------------
  // MacCormack-type scheme
  //----------------------------------------------------------------------------

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

  // center scheme for macdrp, which is used in metric calculation
  int mac_center_all_coef_size = 15;
  size_t mac_center_all_info[mac_max_half_stentil+1][FD_INFO_SIZE] =
    {
      // POS, TOTAL, HALF, LEFT, RIGHT
      { 0, 0, 0, 0, 0 }, // 0 free, not used
      { 0, 3, 1, 1, 1 }, // 1 to free
      { 2, 5, 2, 2, 2 }, // 2
      { 5, 7, 3, 3, 3 }  // 3, normal op for cur scheme
    };
  float mac_center_all_indx[mac_center_all_coef_size]
  {
    -1, 0, 1,
    -2,-1,0, 1, 2,
    -3,-2,-1, 0, 1, 2, 3
  };
  float mac_center_all_coef[10] =
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
};

/*******************************************************************************
 * block structure
 ******************************************************************************/

struct fd_blk_t
{
    //
    // grid index
    //
    int nx, ny, nz;
    int ni1, ni2, nj1, nj2, nk1, nk2, ni, nj, nk;
    int gni1, gnj1, gnk1; // global index

    // size of a single var
    size_t siz_line;
    size_t siz_slice;
    size_t siz_volume; // number of points per var

    //size_t siz_vars; // volume * num_of_vars, not easy for understand, may named with w3d and aux

    //
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
    size_t *force_info; // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
    float  *force_vec_stf;
    int     num_of_moment;
    size_t *moment_info; // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
    float  *moment_ten_rate; // stage, it, Mij, num

    //
    // abs
    //

    // only support one type each run, may combine in the future
    //  which requires different vars in blk_t
    int  abs_itype;

    // only for this block, may diff with global values from input par
    int  abs_num_of_layers[ FD_NDIM_2 ];

    // grid index of each face
    size_t abs_indx[FD_NDIM_2 * FD_NDIM_2];

    size_t   abs_coefs_facepos0[FD_NDIM_2];  // 
    float   *abs_coefs; // all coefs all faces 

    // may add abs_numb_of_vars later for mpml
    //   may not need it since 

    size_t  abs_vars_volsiz[FD_NDIM_2]; // size of single var in abs_blk_vars for each block
    size_t  abs_vars_facepos0[FD_NDIM_2]; // start pos of each blk in first level; other level similar
    size_t  abs_vars_size_per_level; // size of vars per level
    float   *abs_vars; //  order: vars_one_block -> all block -> one level -> 4 levels

    // io
    int num_of_sta;
    size_t *sta_loc_point;
    float  *sta_seismo;

    int num_of_snap;
    int *snap_grid_indx;
    int *snap_time_indx;

    // dir
    char out_dir[FD_MAX_STRLEN];

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


/*******************************************************************************
 * mpi info for each process
 ******************************************************************************/

struct fd_mpi_t
{
  int myid2[2];
  int neighid[FD_NDIM_2];
  MPI_Comm topocomm;
};

/*******************************************************************************
 * function prototype
 ******************************************************************************/

void 
fd_set_macdrp(struct fd_t *fd);

void
fd_blk_init(struct fd_blk_t *blk,
      struct fd_blk_t *blk,
      int number_of_total_grid_points_x,
      int number_of_total_grid_points_y,
      int number_of_total_grid_points_z,
      int number_of_mpiprocs_x,
      int number_of_mpiprocs_y,
      char **boundary_type_name,
      int *abs_number_of_layers,
      int fdx_nghosts,
      int fdy_nghosts,
      int fdz_nghosts,
      int number_of_levels, // depends on time scheme, for rk4 = 4
      int *myid2,
      int *neighid,
      const int myid, const int verbose);

void
fd_mpi_create_topo(struct fd_mpi_t *fdmpi, int myid, MPI_Comm comm, int nprocx, int nprocy);

#endif
