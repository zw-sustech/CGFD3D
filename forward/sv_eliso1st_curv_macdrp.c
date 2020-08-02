
/*********************************************************************
 This is the main program  for multi-block forward modeling: 

 AUTHOR:
     ZANG Nan
    ZHANG Wei

    elastic_iso_curv_col.c
    elastic_iso_curv_leb.c
    elastic_vti_curv_leb.c
    elastic_tti_curv_leb.c

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "netcdf.h"

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "par_funcs.h"
#include "blk_struct.h"
#include "scheme_struct.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "solver_eliso1st_curv_macdrp.h"

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

int sv_eliso1st_curv_macdrp_allstep(
    float *restrict w3d, // wavefield
    float *restrict g3d, // grid vars
    float *restrict m3d, // medium vars
    // grid size
    size_t nx, size_t ny, size_t nz,
    // time
    float dt, int nt_total,
    // scheme
    int numPair, int numStage, float *rka, float *rkb,
    int fd_max_op_len,
    size_t ***pair_fdx_all_len, size_t **pair_fdx_all_indx, float **pair_fdx_all_coef,
    size_t ***pair_fdy_all_len, size_t **pair_fdx_all_indx, float **pair_fdy_all_coef,
    size_t ***pair_fdz_all_len, size_t **pair_fdz_all_indx, float **pair_fdz_all_coef,
    // boundary type
    int *restrict boundary_itype
    // if free surface
    float *matVx2Vz, float *matVy2Vz,
    // if abs
    int *restrict abs_numbers, float *restrict abs_coefs, float *restrict abs_vars,
    // source term
    int num_of_force, int *restrict force_loc_point, float *restrict force_term,
    int num_of_moment, int *restrict moment_loc_point, float *restrict moment_term,
    // mpi
    int myid, int *myid3,
    // io
    int num_of_sta, int *restrict sta_loc_point, float *restrict sta_waveform,
    int num_of_snap, int *restrict snap_indx, int *restrict snap_time_indx,
    char *out_dir)
{
  int ierr = 0;

  int it; // loop var for time

  size_t (*a_cur)[ FD_NDIM_2 ];
  size_t (*a_pre)[ FD_NDIM_2 ];

  pos_pre = 0 * siz_vars;
  pos_rhs = 1 * siz_vars;
  pos_cur = 2 * siz_vars; // current value for rhs
  pos_end = 3 * siz_vars;

  // pml vars for rk
  for (int i=0; i< FD_NDIM_2; i++) {
    a_pre[i] = abs_vars + abs_vars_dimpos[i];
    a_end[i] = a_pre[i];
    if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML) {
      a_end[i] = a_pre[i] + abs_vars_dimsiz[i];
    }
  }

  while ( it < nt_total )
  {
    if (myid==0) fprintf(stdout,"-> it=%d, t=\%f\n", it, (it-0)*stept);

    iPair = it % num_of_pair; // mod to ipair
    if (myid==0) fprintf(stdout, " --> iPair=%d\n",iPair);

    //// init RK
    //ierr = curv_macdrp_elastic_iso_RK_start();

    //if (par->absorbing_boundary_adecfspml==1) {
    //    ierr = curv_macdrp_elastic_iso_pml_RK_start();
    //}
    //

    scheme_op_lenth = 5;
    scheme_op_lenth_left  = 3;
    scheme_op_lenth_right = 1;

    int    *fd_bdry_opx_nlay; // x1, x2
    int    *fd_bdry_opx_length; // size: fd_bdryop_nlay[0] + [1]
    int    *fd_bdry_opx_pos_in_shift;
    size_t *fd_bdry_opx_shift; // 
    float  *fd_bdry_opx_coef;

    int fd_bdry_opx_nlay;
    int fd_bdry_opx_nlay;

    w_pre = w3d+pos_pre;
    w_end = w3d+pos_end;
    w_rhs = w3d+pos_rhs;

    // loop RK stages for one step
    for (iStage=0; iStage<numStage; iStage++)
    {
      if (myid==0) fprintf(stdout, " --> iStage=%d\n",iStage);

      if (iStage==0) {
          w_cur = w3d_pre;
      } else {
          w_cur = w3d+pos_cur;
      }

      // pack message
      ierr = wf_el_1st_pack_message(W_cur,
              ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,
              pair_fdx_len[ipar][istage][fd_max_op_len],
              buff, buff_dim_pos, buff_dim_size);

      // exchange data along y
      MPI_Request r_reqs[4];

      MPI_Startall(2, reqs+0, MPI_STATUS_IGNORE);

      // wait for communication
      MPI_Waitall(2, reqs+0, MPI_STATUS_IGNORE);

      // receive message and uppack
      ierr = wf_el_1st_unpack_messag();

      // compute
      ierr = sv_eliso1st_curv_macdrp_one_stage(
              w3d+pos_pre,
              w3d+pos_rhs,
              pair_fdx_len[ipar][istage],
              fd_opx_length,
              fd_opx_shift,
              fd_opx_coef,
              );

      // RK int for w

      // cal w_cur first, then start send etc; or w_end at last stage
      if (iStage<numStage-1)
      {
          // wavefield
          for (iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
              w_cur[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
          }
          // apply exp abl
          for (int i=0; i < FD_NDIM_2; i++) {
            if (boundary_itype[i] == FD_BOUNDARY_TYPE_ABLEXP) {
              sv_eliso1st_curv_macdrp_apply_ablexp(w_cur);
            }
          }
          // pack and isend
          ierr = pack_message(w_cur,sbuff);
          MPI_Startall(sreqs);
      }
      else // w_end for send at last stage
      {
          // wavefield
          for (iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
              w_end[iptr] += coef_b * w_rhs[iptr];
          }
          // apply exp abl
          for (int i=0; i < FD_NDIM_2; i++) {
            if (boundary_itype[i] == FD_BOUNDARY_TYPE_ABLEXP) {
              sv_eliso1st_curv_macdrp_apply_ablexp(w_end);
            }
          }
          // pack and isend
          ierr = pack_message(w_end,sbuff);
          MPI_Startall(sreqs);
      }

      // cal w_end at other stages, non-communicate vars
      if (iStage==0) {
          pos_var = 0;
          for (ivar=0; ivar<num_of_vars; ivar++) {
              pos_var += size_volume;
              pos_val = pos_var;
              for (iptr=0; iptr<siz_volume; iptr++)
              {
                  pos_val += 1;
                  w_end[pos_val] = w_pre[pos_val] + coef_b * w_rhs[pos_val];
              }
          }
      } else if (iStage<numStage-1) {
          pos_var = 0;
          for (ivar=0; ivar<num_of_vars; ivar++) {
              pos_var += size_volume;
              pos_val = pos_var;
              for (iptr=0; iptr<siz_volume; iptr++)
              {
                  pos_val += 1;
                  w_end[pos_val] += coef_b * w_rhs[pos_val];
              }
          }
      } 

      // RK int for pml
      for (i=0; i<6; i++) {
        if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML) {
          a_cur = abs_vars_dimpos[i];
          // cur var
          if (istage < num_of_stage -1) {
            for (iptr=0; iptr<siz_aux_allvars; iptr++) {
                a_cur[iptr] = a_pre[iptr] + coef_a * a_rhs[iptr];
            }
          }
          // end var
          if (istage==0) {
            a_end[iptr] = a_pre[iptr] + coef_b * a_rhs[iptr];
          } else {
            a_end[iptr] += coef_b * a_rhs[iptr];
          }
        }
      }

    } // RK stages

    // swap w_pre and w_end, avoid copying
    w_tmp = w_pre; w_pre = w_end; w_end = w_tmp;
    a_tmp = a_pre; a_pre = a_end; a_end = a_tmp;

    // time step + 1
    it++;

    // QC
    if (par->qc_check_value_ornot==1 && id/par->qc_check_value_tinv==0) {
        ierr = check_value();
    }

    // save results
    ierr = io_seismo_keep();
    ierr = io_snapshot_save();

  } // time loop
}


int sv_eliso1st_curv_macdrp_onestage(float *restrict w_cur, float *restrict rhs, 
        int *restrict fdx_length, // 0: total, 1: left, 2: right
        // include different order/stentil
        int number_of_pairs, int number_of_stages,
        int fdx_max_hlen, int fdy_max_hlen, int fdz_max_hlen,
        size_t **restrict fdx_all_info, size_t *restrict fdx_all_indx,float *restrict fdx_all_coef,
        size_t **restrict fdy_all_info, size_t *restrict fdy_all_indx,float *restrict fdy_all_coef,
        size_t **restrict fdz_all_info, size_t *restrict fdz_all_indx,float *restrict fdz_all_coef,
        int myid, int *restrict myid3
        )
{
  int ierr = 0;

  size_t i,j,k;

  float *restrict Vx;
  float *restrict Vy;
  float *restrict Vz;

  float  fd_opx_coef[fdx_max_hlen];
  size_t fd_opx_shift[fdx_max_hlen];

  iptr = fdx_dhs_indx[fd_op_half_len-1];
  for (i=0; i<fd_op_total_len; i++) {
    fd_opx_coef[i] = fdx_opcoef[iptr+i];
    fd_opx_indx[i] = fdx_opindx[iptr+i];
  }

  // inner points
  i_in_opindx = fdx_all_len[fdx_max_half_len][0];
  len_of_op   = fdx_all_len[fdx_max_half_len][1];
  ierr = eliso1st_curv_macdrp_rhs_inner(w_cur, rhs,
          ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,
          fdx_op_len, fdx_all_indx+len_of_op,
          fdy_op_len, fdy_all_indx+len_of_op,
          fdz_op_len, fdz_all_indx+len_of_op,
          )

  // for boundary, free surface should be called last, which can correct other boundary affacted by free

  // boundary x1
  switch (boundary_itype[0])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 
      ierr = eliso1st_curv_macdrp_rhs_cfspml_x();
      break;

    default:
  }

  // boundary x2
  switch (boundary_itype[1])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 
      ierr = eliso1st_curv_macdrp_rhs_cfspml_x();
      break;

    default:
  }

  // boundary y1
  switch (boundary_itype[2])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 
      ierr = eliso1st_curv_macdrp_rhs_cfspml_y();
      break;

    default:
  }

  // boundary y2
  switch (boundary_itype[3])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 
      ierr = eliso1st_curv_macdrp_rhs_cfspml_y();
      break;

    default:
  }

  // boundarz y1
  switch (boundary_itype[4])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 
      ierr = eliso1st_curv_macdrp_rhs_cfspml_z();
      break;

    default:
  }

  // boundary z2
  switch (boundary_itype[5]) {
    case FD_BOUNDARY_TYPE_FREE :
      // tractiong
      ierr = sv_eliso1st_curv_macdrp_rhs_timg_z2();
      // velocity: vlow
      ierr = sv_eliso1st_curv_macdrp_rhs_vlow_z2();
      break;

    case FD_BOUNDARY_TYPE_CFSPML : 
      ierr = eliso1st_cfspml_cal_rhs_z();
      break;

    default:
  }

  // source term
  switch (source_itype) {
    case SOURCE_TYPE_POINT :
      ierr = curv_macdrp_eliso1st_src_point();
      break;

    case SOURCE_TYPE_GAUSS :
      ierr = curv_macdrp_eliso1st_src_gauss();
      break;

    case SOURCE_TYPE_SINC :
      ierr = curv_macdrp_eliso1st_src_sinc();
      break;
  }

  return ierr;
}

/*
 * calculate inner points without considering the boundaries
 */

int eliso1st_curv_macdrp_rhs_inner(
    float *restrict w_cur, float *restrict rhs, 
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    )
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // local pointer
  float *restrict Vx;
  float *restrict Vy;
  float *restrict Vz;
  float *restrict Txx;
  float *restrict Tyy;
  float *restrict Tzz;
  float *restrict Txz;
  float *restrict Tyz;
  float *restrict Txy;

  // local var
  size_t i,j,k;
  size_t n_fd; // loop var for fd

  float DxVx, // etc

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // get each vars
  Vx  = w_cur + WF_EL_1ST_SEQ_VX;
  Vy  = w_cur + WF_EL_1ST_SEQ_VY;
  Vz  = w_cur + WF_EL_1ST_SEQ_VZ;
  Txx = w_cur + WF_EL_1ST_SEQ_TXX;
  Tyy = w_cur + WF_EL_1ST_SEQ_TYY;
  Tzz = w_cur + WF_EL_1ST_SEQ_TZZ;
  Txz = w_cur + WF_EL_1ST_SEQ_TXZ;
  Tyz = w_cur + WF_EL_1ST_SEQ_TYZ;
  Txy = w_cur + WF_EL_1ST_SEQ_TXY;

  xi_x = g3d + GRID_CURV_SEQ_XIX;
  xi_y = g3d + GRID_CURV_SEQ_XIY;
  xi_z = g3d + GRID_CURV_SEQ_XIZ;

  lambda3d = m3d + MD_ELISO_SEQ_LAMBDA;

  // loop all points
  for (k=nk1; k<nk2; k++)
  {
    iptr_k = k * siz_slice;
    for (j=nj1; j<=nj2; j++)
    {
      iptr_j = iptr_k + j * siz_line;

      iptr = iptr_j + (ni1-1);

      for (i=ni1; i<ni2; i++)
      {
        iptr += 1;

        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];

        // medium
        lam = lambda3d[iptr];
        mu  = mu3d[iptr];
        rho = rho3d[iptr];

        lam2mu = lam + 2.0 * mu;
        rrho = 1.0 / rho;

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVx, Vx, iptr, fd_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzVx, Vx, iptr, fd_len, lfdz_shift, lfdz_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT(DxVy, Vy, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVy, Vy, iptr, fd_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzVy, Vy, iptr, fd_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVz, Vz, iptr, fd_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzVz, Vz, iptr, fd_len, lfdz_shift, lfdz_coef, n_fd);

        // need to add Tij


        // moment equation
        hVx[iptr] = ( DxTxx*xix + DxTxy*xiy + DxTxz*xiz
                    + DyTxx*etx + DyTxy*ety + DyTxz*etz
                    + DzTxx*ztx + DzTxy*zty + DzTxz*ztz ) * rrho;

        // Hooke's equatoin
        hTxx[iptr] = lam2mu*DxVx*xix + lam*DxVy*xiy + lam*DxVz*xiz
                   + lam2mu*DyVx*etx + lam*DyVy*ety + lam*DyVz*etz
                   + lam2mu*DzVx*ztx + lam*DzVy*zty + lam*DzVz*ztz;
      }
    }
  }
}

/*
 * implement traction image boundary plus possible ade cfs-pml
 *  any eqn affected by rhs near surface should be reprocessed in this func
 *  other eqn should be called before free surface
 */

int eliso1st_curv_macdrp_rhs_timg_z2(
    float *restrict w_cur, float *restrict rhs, 
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype, int *restrict abs_number_of_layers,
    float *restrict *abs_coefs,
    float *restrict *abs_vars
    )
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // local pointer
  float *restrict Vx;
  float *restrict Vy;
  float *restrict Vz;
  float *restrict Txx;
  float *restrict Tyy;
  float *restrict Tzz;
  float *restrict Txz;
  float *restrict Tyz;
  float *restrict Txy;

  // local var
  size_t i,j,k;
  size_t n_fd; // loop var for fd
  size_t fdz_len;

  float DxVx, // etc

  // to save traction and other two dir var
  float vecxi[fdz_len];
  float vecet[fdz_len];
  float veczt[fdz_len];

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // get each vars
  Vx  = w_cur + WF_EL_1ST_SEQ_VX;
  Vy  = w_cur + WF_EL_1ST_SEQ_VY;
  Vz  = w_cur + WF_EL_1ST_SEQ_VZ;
  Txx = w_cur + WF_EL_1ST_SEQ_TXX;
  Tyy = w_cur + WF_EL_1ST_SEQ_TYY;
  Tzz = w_cur + WF_EL_1ST_SEQ_TZZ;
  Txz = w_cur + WF_EL_1ST_SEQ_TXZ;
  Tyz = w_cur + WF_EL_1ST_SEQ_TYZ;
  Txy = w_cur + WF_EL_1ST_SEQ_TXY;

  xi_x = g3d + GRID_CURV_SEQ_XIX;
  xi_y = g3d + GRID_CURV_SEQ_XIY;
  xi_z = g3d + GRID_CURV_SEQ_XIZ;

  lambda3d = m3d + MD_ELISO_SEQ_LAMBDA;

  // get pml coefs

  int k_min = nk2 - fdz_indx[fd_len-1]; // last indx, free surface force Tx/Ty/Tz to 0 in cal
  //int k_min = nk2 + fdz_indx[0]; // first indx is negative, so + 
  //int k_max = nk2 + fdz_indx[fd_len-1]; // last indx

  // conver to op index

  // loop all points

  //for (n=0; n < fdz_len; n++)
  for (k=k_min; k <= nk2; k++)
  {
    // k corresponding to 0 index of the fd op

    // index of free surface
    size_t n_free = nk2 - k - fdz_indx[0]; // first indx is negative

    // for 1d index
    iptr_k = k * siz_slice;

    for (j=nj1; j<=nj2; j++)
    {
      iptr_j = iptr_k + j * siz_line;

      iptr = iptr_j + (ni1-1);

      for (i=ni1; i<ni2; i++)
      {
          iptr += 1;

          // metric
          xix = xi_x[iptr];
          xiy = xi_y[iptr];
          xiz = xi_z[iptr];

          // medium
          lam = lambda3d[iptr];
          mu  = mu3d[iptr];
          rho = rho3d[iptr];

          lam2mu = lam + 2.0 * mu;
          rrho = 1.0 / rho;

          //
          // for hVx
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txx[iptr4vec]
                                          + xi_y[iptr4vec] * Txy[iptr4vec]
                                          + xi_z[iptr4vec] * Txz[iptr4vec] );
          }
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txx[iptr4vec]
                                          + et_y[iptr4vec] * Txy[iptr4vec]
                                          + et_z[iptr4vec] * Txz[iptr4vec] );
          }
          //for (n=0; n<fdx_len; n++) {
          //  iptr4vec = iptr + fdx_indx[n] * siz_slice;
          //  veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
          //                                + zt_y[iptr4vec] * Txy[iptr4vec]
          //                                + zt_z[iptr4vec] * Txz[iptr4vec] );
          //}
          // blow surface, cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdx_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                          + zt_y[iptr4vec] * Txy[iptr4vec]
                                          + zt_z[iptr4vec] * Txz[iptr4vec] );
          }
          // at surface, set to 0
          vecet[n_free] = 0.0;
          // above surface, mirror
          for (n=n_free+1; n<fdz_len; n++) {
            veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // hVx 
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVx[iptr] = ( DxTx+DyTy+DzTz )*rrhojac;

          // cfspml
          if (boundary_itype[0] == FD_BOUNDARY_TYPE_CFSPML && i<ni1+abs_number_of_layers[1])
          {
            // get aux vars
            siz_aux_volume = number_of_layers * nj * nk;
            a_Vx = aux_var + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
            a_Vy = aux_var + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
            a_Vz = aux_var + WF_EL_1ST_SEQ_PMLVZ * siz_aux_volume;

            // index of each auxiliary var
            iptr_a = (i-ni1) + j * abs_ni + k * nj * abs_ni;

            d1=Dx[i]; b1=Bx[i]; a1=Ax[i];
            rb1 = 1.0 / b1;

            // 1: make corr to moment equation
            hVx[iptr] += (rb1 - 1.0) * DxTx * rrho
                          - rb1 * a_Txx[iptr_a];

            // make corr to Hooke's equatoin
            
            // 2: aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            a_hTxx[iptr_a] = d1*(lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) -a1*a_Txx[iptr_a];
          }
          // other x2, y1, y2


          //
          // for hVy, hVz
          //

      }
    }
  }
}

/*
 * implement vlow boundary plus possible ade cfs-pml
 */

int eliso1st_curv_macdrp_rhs_vlow_z2(
    float *restrict w_cur, float *restrict rhs, 
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_max_half_len, size_t fdz_all_indx_size,
    size_t* fdz_all_info, size_t **restrict fdz_all_indx, float **restrict fdz_all_coef,
    )
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];

  // local var and pointer for different hlen along k dim
  size_t  lfdz_pos; // pos in _indx
  size_t  lfdz_len;

  // point to different hlen fz near surface
  size_t *restrict lfdz_info;
  size_t *restrict lfdz_indx;
  size_t *restrict lfdz_shift;
  float  *restrict lfdz_coef;

  // local pointer
  float *restrict Vx;
  float *restrict Vy;
  float *restrict Vz;
  float *restrict Txx;
  float *restrict Tyy;
  float *restrict Tzz;
  float *restrict Txz;
  float *restrict Tyz;
  float *restrict Txy;

  // local var
  size_t i,j,k;
  size_t n_fd; // loop var for fd
  size_t fdz_len;

  float DxVx, // etc

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }

  // get each vars
  Vx  = w_cur + WF_EL_1ST_SEQ_VX;
  Vy  = w_cur + WF_EL_1ST_SEQ_VY;
  Vz  = w_cur + WF_EL_1ST_SEQ_VZ;
  Txx = w_cur + WF_EL_1ST_SEQ_TXX;
  Tyy = w_cur + WF_EL_1ST_SEQ_TYY;
  Tzz = w_cur + WF_EL_1ST_SEQ_TZZ;
  Txz = w_cur + WF_EL_1ST_SEQ_TXZ;
  Tyz = w_cur + WF_EL_1ST_SEQ_TYZ;
  Txy = w_cur + WF_EL_1ST_SEQ_TXY;

  xi_x = g3d + GRID_CURV_SEQ_XIX;
  xi_y = g3d + GRID_CURV_SEQ_XIY;
  xi_z = g3d + GRID_CURV_SEQ_XIZ;

  lambda3d = m3d + MD_ELISO_SEQ_LAMBDA;

  // index to use lower-order scheme

  lfdz_info  = fdz_all_info[fdz_max_half_len];
  lfdz_left  = lfdz_info[FD_INFO_LENGTH_LEFT];
  lfdz_right = lfdz_info[FD_INFO_LENGTH_RIGTH];

  //int k_min = nk2 - lfdz_right + 1;
  //int k_max = nk2 - 1;

  // loop all points

  //for (k=nk2 - fdz_max_half_len + 1; k <= nk2; k++)
  for (n=1; n <= fdz_max_half_len - 1; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // fdz of lower order
    lfdz_info  = fdz_all_info[n];
    lfdz_pos   = lfdz_info[FD_INFO_POS_OF_INDX];
    lfdz_len   = lfdz_info[FD_INFO_LENGTH_TOTAL];
    lfdz_indx  = fdz_all_indx[lfdz_npos];
    lfdz_coef  = fdz_all_coef[lfdz_npos];
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_shift[n_fd] = lfdz_indx[n_fd] * siz_slice;
    }

    // for 1d index
    iptr_k = k * siz_slice;

    for (j=nj1; j<=nj2; j++)
    {
      iptr_j = iptr_k + j * siz_line;

      iptr = iptr_j + (ni1-1);

      for (i=ni1; i<ni2; i++)
      {
          iptr += 1;

          // metric
          xix = xi_x[iptr];
          xiy = xi_y[iptr];
          xiz = xi_z[iptr];

          // medium
          lam = lambda3d[iptr];
          mu  = mu3d[iptr];
          rho = rho3d[iptr];

          lam2mu = lam + 2.0 * mu;
          rrho = 1.0 / rho;

          // Vx derivatives
          M_FD_SHIFT(DxVx, Vx, iptr, fd_length, lfdx_shift, lfdx_coef, n_fd);
          M_FD_SHIFT(DyVx, Vx, iptr, fd_length, lfdy_shift, lfdy_coef, n_fd);
          M_FD_SHIFT(DzVx, Vx, iptr, fd_length, lfdz_shift, lfdz_coef, n_fd);

          // Vy derivatives
          M_FD_SHIFT(DxVy, Vy, iptr, fd_length, lfdx_shift, lfdx_coef, n_fd);
          M_FD_SHIFT(DyVy, Vy, iptr, fd_length, lfdy_shift, lfdy_coef, n_fd);
          M_FD_SHIFT(DzVy, Vy, iptr, fd_length, lfdz_shift, lfdz_coef, n_fd);

          // Vz derivatives
          M_FD_SHIFT(DxVz, Vz, iptr, fd_length, lfdx_shift, lfdx_coef, n_fd);
          M_FD_SHIFT(DyVz, Vz, iptr, fd_length, lfdy_shift, lfdy_coef, n_fd);
          M_FD_SHIFT(DzVz, Vz, iptr, fd_length, lfdz_shift, lfdz_coef, n_fd);

          // need to add Tij


          // moment equation
          hVx[iptr] = ( DxTxx*xix + DxTxy*xiy + DxTxz*xiz
                      + DyTxx*etx + DyTxy*ety + DyTxz*etz
                      + DzTxx*ztx + DzTxy*zty + DzTxz*ztz ) * rrho;

          // Hooke's equatoin
          hTxx[iptr] = lam2mu*DxVx*xix + lam*DxVy*xiy + lam*DxVz*xiz
                     + lam2mu*DyVx*etx + lam*DyVy*ety + lam*DyVz*etz
                     + lam2mu*DzVx*ztx + lam*DzVy*zty + lam*DzVz*ztz;

          // cfspml
          if (boundary_itype[0] == FD_BOUNDARY_TYPE_CFSPML && i<ni1+abs_number_of_layers[1])
          {
            // get aux vars
            siz_aux_volume = number_of_layers * nj * nk;
            a_Vx = aux_var + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
            a_Vy = aux_var + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
            a_Vz = aux_var + WF_EL_1ST_SEQ_PMLVZ * siz_aux_volume;

            // index of each auxiliary var
            iptr_a = (i-ni1) + j * abs_ni + k * nj * abs_ni;

            d1=Dx[i]; b1=Bx[i]; a1=Ax[i];
            rb1 = 1.0 / b1;

            // 1: make corr to moment equation
            hVx[iptr] += (rb1 - 1.0) * DxTx * rrho
                          - rb1 * a_Txx[iptr_a];

            // make corr to Hooke's equatoin
            
            // 2: aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            a_hTxx[iptr_a] = d1*(lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) -a1*a_Txx[iptr_a];
          }
          // other x2, y1, y2

      }
    }
  }
}

/*
 * cfspml along x
 */

int eliso1st_curv_macdrp_rhs_cfspml_x(
    float *restrict w_var, float *restrict rhs,
    float *restrict aux_var, float *restrict aux_rhs,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t abs_ni1, abs_ni2, abs_ni,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    float *restrict Ax, float *restrict Bx, float *restrict Dx,
    float *restrict abs_vars
    )
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // local pointer
  float *restrict Vx;
  float *restrict Vy;
  float *restrict Vz;
  float *restrict Txx;
  float *restrict Tyy;
  float *restrict Tzz;
  float *restrict Txz;
  float *restrict Tyz;
  float *restrict Txy;

  float *restrict aVx ;
  float *restrict aTxx;

  float *restrict ahVx ;
  float *restrict ahTxx;

  // local var
  size_t i,j,k;
  size_t n_fd; // loop var for fd
  size_t iptr, iptr_a;
  size_t siz_aux_volume;

  float DxVx, // etc

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // get each vars
  Vx  = w_cur + WF_EL_1ST_SEQ_VX * siz_volume;
  Vy  = w_cur + WF_EL_1ST_SEQ_VY * siz_volume;
  Vz  = w_cur + WF_EL_1ST_SEQ_VZ * siz_volume;
  Txx = w_cur + WF_EL_1ST_SEQ_TXX * siz_volume;
  Tyy = w_cur + WF_EL_1ST_SEQ_TYY * siz_volume;
  Tzz = w_cur + WF_EL_1ST_SEQ_TZZ * siz_volume;
  Txz = w_cur + WF_EL_1ST_SEQ_TXZ * siz_volume;
  Tyz = w_cur + WF_EL_1ST_SEQ_TYZ * siz_volume;
  Txy = w_cur + WF_EL_1ST_SEQ_TXY * siz_volume;

  xi_x = g3d + GRID_CURV_SEQ_XIX * siz_volume;
  xi_y = g3d + GRID_CURV_SEQ_XIY * siz_volume;
  xi_z = g3d + GRID_CURV_SEQ_XIZ * siz_volume;

  lambda3d = m3d + MD_ELISO_SEQ_LAMBDA * siz_volume;

  // get aux vars
  siz_aux_volume = number_of_layers * nj * nk;
  a_Vx = aux_var + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
  a_Vy = aux_var + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
  a_Vz = aux_var + WF_EL_1ST_SEQ_PMLVZ * siz_aux_volume;

  // loop all points
  for (k=nk1; k<nk2; k++)
  {
    iptr_k = k * siz_slice;
    //d3=Dz[k];
    //b3=Bz[k];
    //a3=Az[k];
     
    for (j=nj1; j<=nj2; j++)
    {
      iptr_j = iptr_k + j * siz_line;

      iptr = iptr_j + (ni1-1);

      //d2=Dy[j];
      //b2=By[j];
      //a2=Ay[j];

      for (i=abs_ni1; i<=abs_ni2; i++)
      {
        iptr += 1;

        // index of each auxiliary var
        iptr_a = (i-ni1) + j * abs_ni + k * nj * abs_ni;

        d1=Dx[i]; b1=Bx[i]; a1=Ax[i];
        rb1 = 1.0 / b1;

        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];

        // medium
        lam = lambda3d[iptr];
        mu  = mu3d[iptr];
        rho = rho3d[iptr];

        lam2mu = lam + 2.0 * mu;
        rrho = 1.0 / rho;

        // xi derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DxVy, Vy, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DxVz, Vz, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);

        // need to add Tij

        // 1: make corr to moment equation
        hVx[iptr] += (rb1 - 1.0) * ( DxTxx*xix + DxTxy*xiy + DxTxz*xiz) * rrho
                      - rb1 * a_Txx[iptr_a];

        // make corr to Hooke's equatoin
        
        // 2: aux var
        //   a1 = alpha + d / beta, dealt in abs_set_cfspml
        a_hTxx[iptr_a] = d1*(lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz ) -a1*a_Txx[iptr_a];
      }
    }
  }
}

// apply on wavefield, not on derivatives
int sv_eliso1st_curv_macdrp_apply_ablexp()
{
}

