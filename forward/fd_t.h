#ifndef FD_T_H
#define FD_T_H

/*******************************************************************************
 *macro for fd opterators
 *******************************************************************************/

// use siz_shift to find adjacent point of the stentil for 3d var
#define M_FD_SHIFT(deriv, var, iptr, fd_length, fd_shift, fd_coef, n) \
   deriv = fd_coef[0] * var[iptr + fd_shift[0]]; \
   for (n=1; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n]]; \
   }

// use pointer for cur point for speedup
#define M_FD_SHIFT_PTR(deriv, var_ptr, fd_length, fd_shift, fd_coef, n) \
  deriv = fd_coef[0] * *(var_ptr + fd_shift[0]);                                                        \
  for (n = 1; n < fd_length; n++)                                     \
  {                                                                   \
    deriv += fd_coef[n] * *(var_ptr + fd_shift[n]);                    \
  }

// only valid for macdrp etc with len = 5, may be faster? 
#define M_FD_SHIFT_PTR_MACDRP(deriv, var_ptr, fd_length, fd_shift, fd_coef, n) \
  deriv =  fd_coef[0] * *(var_ptr + fd_shift[0])                    \
          +fd_coef[1] * *(var_ptr + fd_shift[1])                    \
          +fd_coef[2] * *(var_ptr + fd_shift[2])                    \
          +fd_coef[3] * *(var_ptr + fd_shift[3])                    \
          +fd_coef[4] * *(var_ptr + fd_shift[4]);

// assume var has the same size as fd_coef, ordered one by one, thus no index needed
#define M_FD_NOINDX(deriv, var, fd_length, fd_coef, n) \
   deriv = fd_coef[0] * var[0]; \
   for (n=1; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[n]; \
   }

// use indx relative to cur point as (-1,0,1), need to multiply siz_shift for 3d array
#define M_FD_INDX(deriv, var, iptr, fd_length, fd_indx, fd_coef, shift, n) \
   deriv = fd_coef[0] * var[iptr + fd_shift[0] * shift]; \
   for (n=1; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n] * shift]; \
   }

/*******************************************************************************
 * structure for different fd schemes
 ******************************************************************************/

/*
 * elementary operator
 */

typedef struct
{
  int total_len;
  int half_len;
  int left_len;
  int right_len;
  int   *indx;  // indx change to cur point as 0 for 1d
  int   *shift; // num of grid points skipped
  float *coef;
} fd_op_t; 

/*
 * collocated grid scheme
 */

typedef struct {

  //----------------------------------------------------------------------------
  // Runge-Kutta time scheme
  //----------------------------------------------------------------------------

  int num_rk_stages;
  float *rk_a;
  float *rk_b;
  float *rk_rhs_time; // relative time for rhs eval

  //----------------------------------------------------------------------------
  // central scheme
  //----------------------------------------------------------------------------

  int     fdc_len;
  int     fdc_half_len;
  int     fdc_nghosts;
  int    *fdc_indx;
  float  *fdc_coef;

  //----------------------------------------------------------------------------
  // para for different schemes at points to boundaries for different dim
  //----------------------------------------------------------------------------

  // ghost point required 
  int fdx_nghosts;
  int fdy_nghosts;
  int fdz_nghosts;

  // max total len of op
  int fdx_max_len;
  int fdy_max_len;
  int fdz_max_len;

  // max half len
  int fdx_max_half_len;
  int fdy_max_half_len;
  int fdz_max_half_len;

  ////----------------------------------------------------------------------------
  //// center schemes at different points to boundary
  ////----------------------------------------------------------------------------

  //// center scheme can be separated into 1d array as
  //int    fdx_max_nlay;
  //int    fdx_max_total;
  //int   *fdx_num_total;
  //int   *fdx_num_half;
  //int   *fdx_num_left;
  //int   *fdx_num_right;
  //int   *fdx_indx; // fdx_max_total * fdx_max_nlay
  //float *fdx_coef;

  ////  fd_ is 2d pointer, the first pointer means the grid layer to free surface (from 0 to fd_len-1),
  ////  the second pointer points to fd op for that layer, the size could be different, larger than
  ////  inner op
  //// not used yet
  //int    **fdx_all_info; // [k2free][pos, total, half, left, right]
  //int    **fdx_all_indx;
  //float  **fdx_all_coef;

  //int    **fdy_all_info;
  //int    **fdy_all_indx;
  //float  **fdy_all_coef;

  //int    **fdz_all_info;
  //int    **fdz_all_indx;
  //float  **fdz_all_coef;

  ////----------------------------------------------------------------------------
  //// filter schemes at different points to boundary
  ////----------------------------------------------------------------------------

  //// not used yet

  //int    **filtx_all_info; // [k2free][pos, total, half, left, right] 
  //int    **filtx_all_indx;
  //float  **filtx_all_coef;

  //int    **filty_all_info;
  //int    **filty_all_indx;
  //float  **filty_all_coef;

  //int    **filtz_all_info;
  //int    **filtz_all_indx;
  //float  **filtz_all_coef;

  //----------------------------------------------------------------------------
  // pairs for 3d space for MacCormack-type schemes
  //----------------------------------------------------------------------------

  int num_of_pairs;

  // number of layers that need to use biased op near boundary
  int num_of_fdx_op;
  int num_of_fdy_op;
  int num_of_fdz_op;

  fd_op_t ***pair_fdx_op; // [pair][stage][nlay]
  fd_op_t ***pair_fdy_op;
  fd_op_t ***pair_fdz_op;

} fd_t;

/*
 * staggered grid scheme
 */

typedef struct {

  // ghost point required 
  int fdx_nghosts;
  int fdy_nghosts;
  int fdz_nghosts;

  // max total len of op
  int fdx_max_len;
  int fdy_max_len;
  int fdz_max_len;

  // max half len
  int fdx_max_half_len;
  int fdy_max_half_len;
  int fdz_max_half_len;

  // number of layers that need to use biased op near boundary
  int num_of_fdx_op;
  int num_of_fdy_op;
  int num_of_fdz_op;

  fd_op_t *Vx_fdx_op; // [nlay]
  fd_op_t *Vx_fdy_op;
  fd_op_t *Vx_fdz_op;

  fd_op_t *Vy_fdx_op; // [nlay]
  fd_op_t *Vy_fdy_op;
  fd_op_t *Vy_fdz_op;

  fd_op_t *Vz_fdx_op; // [nlay]
  fd_op_t *Vz_fdy_op;
  fd_op_t *Vz_fdz_op;

  fd_op_t *Txx_fdx_op; // [nlay]
  fd_op_t *Txx_fdy_op;
  fd_op_t *Txx_fdz_op;

  fd_op_t *Tyy_fdx_op; // [nlay]
  fd_op_t *Tyy_fdy_op;
  fd_op_t *Tyy_fdz_op;

  fd_op_t *Tzz_fdx_op; // [nlay]
  fd_op_t *Tzz_fdy_op;
  fd_op_t *Tzz_fdz_op;

  fd_op_t *Tyz_fdx_op; // [nlay]
  fd_op_t *Tyz_fdy_op;
  fd_op_t *Tyz_fdz_op;

  fd_op_t *Txz_fdx_op; // [nlay]
  fd_op_t *Txz_fdy_op;
  fd_op_t *Txz_fdz_op;

  fd_op_t *Txy_fdx_op; // [nlay]
  fd_op_t *Txy_fdy_op;
  fd_op_t *Txy_fdz_op;

  fd_op_t *lay_fdx_op;
  fd_op_t *lay_fdy_op;
  fd_op_t *lay_fdz_op;

} fdstg_t;

/*******************************************************************************
 * function prototype
 ******************************************************************************/

int 
fd_set_macdrp(fd_t *fd);

void
fd_print(fd_t *fd);

int 
fd_set_stg(fdstg_t *fd);

#endif
