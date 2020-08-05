
/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and macdrp schem
 *
 *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "netcdf.h"

#include "fdlib_math.h"
#include "fd_t.h"
#include "par_funcs.h"
#include "blk_t.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "sv_eliso1st_curv_macdrp.h"

int sv_eliso1st_curv_macdrp_allstep(
    float *restrict w3d, int number_of_wave_vars  , size_t *restrict w3d_pos, // wavefield
    float *restrict g3d, int number_of_grid_vars  , size_t *restrict g3d_pos, // grid vars
    float *restrict m3d, int number_of_medium_vars, size_t *restrict m3d_pos, // medium vars
    // grid size
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    // time
    float dt, int nt1, int nt2, int nt_total,
    // scheme
    int num_rk_stages, float *rk_a, float *rk_b,
    int num_of_pairs, int fd_max_half_stentil,
    size_t ****pair_fdx_all_info, size_t ***pair_fdx_all_indx, float ***pair_fdx_all_coef,
    size_t ****pair_fdy_all_info, size_t ***pair_fdy_all_indx, float ***pair_fdy_all_coef,
    size_t ****pair_fdz_all_info, size_t ***pair_fdz_all_indx, float ***pair_fdz_all_coef,
    // boundary type
    int *restrict boundary_itype,
    // if free surface
    float *matVx2Vz, float *matVy2Vz,
    // if abs
    int **restrict abs_blk_indx, float **restrict abs_blk_coefs, float **restrict abs_blk_vars,
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

  //
  // time loop
  //
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
              sv_eliso1st_curv_macdrp_apply_ablexp(w_cur, number_of_wave_vars,
                  w3d_pos, nx, ny, nz, abs_blk_indx[i],abs_blk_coefs[i]);
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
              sv_eliso1st_curv_macdrp_apply_ablexp(w_end, number_of_wave_vars,
                  w3d_pos, nx, ny, nz, abs_blk_indx[i],abs_blk_coefs[i]);
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
      float *restrict pml_cur, float *restrict pml_rhs,
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

  size_t           fdx_inn_pos;
  size_t *restrict fdx_inn_indx;
  float *restrict  fdx_inn_coef;
  size_t           fdy_inn_pos;
  size_t *restrict fdy_inn_indx;
  float *restrict  fdy_inn_coef;
  size_t           fdz_inn_pos;
  size_t *restrict fdz_inn_indx;
  float *restrict  fdz_inn_coef;

  // local pointer to pml vars
  float *restrict a_Vx ;
  float *restrict a_Txx;

  float *restrict a_hVx ;
  float *restrict a_hTxx;

  fdx_inn_pos  = fdx_all_info[fdx_max_hlen][FD_INFO_POS_OF_INDX];
  fdx_inn_len  = fdx_all_info[fdx_max_hlen][FD_INFO_LENGTH_TOTAL];
  fdx_inn_indx = fdx_all_indx + fdx_inn_pos;
  fdx_inn_coef = fdx_all_coef + fdx_inn_pos;

  fdy_inn_pos  = fdy_all_info[fdy_max_hlen][FD_INFO_POS_OF_INDX];
  fdy_inn_len  = fdy_all_info[fdy_max_hlen][FD_INFO_LENGTH_TOTAL];
  fdy_inn_indx = fdy_all_indx + fdy_inn_pos;
  fdy_inn_coef = fdy_all_coef + fdy_inn_pos;

  fdz_inn_pos  = fdz_all_info[fdz_max_hlen][FD_INFO_POS_OF_INDX];
  fdz_inn_len  = fdz_all_info[fdz_max_hlen][FD_INFO_LENGTH_TOTAL];
  fdz_inn_indx = fdz_all_indx + fdz_inn_pos;
  fdz_inn_coef = fdz_all_coef + fdz_inn_pos;

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


  // inner points
  ierr = eliso1st_curv_macdrp_rhs_inner(w_cur, rhs,
          ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,
          fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
          fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
          fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
          )

  // for boundary, free surface should be called last, which can correct other boundary affacted by free

  // boundary x1
  switch (boundary_itype[0])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 

      // get aux vars and coes
      siz_aux_volume = abs_blk_vars_siz_volume[0];

      ptr_pml = pml_cur + abs_blk_vars_blkpos[0];
      a_Vx = ptr_pml + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
      a_Vy = ptr_pml + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
      // add a_Vz etc

      ptr_pml = pml_rhs + abs_blk_vars_blkpos[0];
      a_hVx = ptr_pml + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
      a_hVy = ptr_pml + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
      // add a_Vz etc

      pml_A = abs_blk_coefs[0];
      pml_B = pml_A + abs_num_of_layers[0];
      pml_D = pml_B + abs_num_of_layers[0];

      ierr = eliso1st_curv_macdrp_rhs_cfspml_x(
          xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z, lam3d, mu3d, rho3d,
          Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
          hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
          a_Vx,a_Vy,a_Vz,a_Txx,a_Tyy,a_Tzz,a_Txz,a_Tyz,a_Txy,
          a_hVx,a_hVy,a_hVz,a_hTxx,a_hTyy,a_hTzz,a_hTxz,a_hTyz,a_hTxy,
          siz_line, siz_volume,
          abs_blk_indx[0], pml_A, pml_B, pml_D,
          fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
          fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
          fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
          );

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
  switch (boundary_itype[5])
  {
    case FD_BOUNDARY_TYPE_CFSPML : 

      // get aux vars and coes
      siz_aux_volume = abs_blk_vars_siz_volume[5];

      ptr_pml = pml_cur + abs_blk_vars_blkpos[5];
      a_Vx = ptr_pml + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
      a_Vy = ptr_pml + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
      // add a_Vz etc

      ptr_pml = pml_rhs + abs_blk_vars_blkpos[5];
      a_hVx = ptr_pml + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
      a_hVy = ptr_pml + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
      // add a_Vz etc

      pml_A = abs_blk_coefs[5];
      pml_B = pml_A + abs_num_of_layers[5];
      pml_D = pml_B + abs_num_of_layers[5];

      ierr = eliso1st_curv_macdrp_rhs_cfspml_z(
          xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z, lam3d, mu3d, rho3d,
          Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
          hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
          a_Vx,a_Vy,a_Vz,a_Txx,a_Tyy,a_Tzz,a_Txz,a_Tyz,a_Txy,
          a_hVx,a_hVy,a_hVz,a_hTxx,a_hTyy,a_hTzz,a_hTxz,a_hTyz,a_hTxy,
          siz_line, siz_volume,
          abs_blk_indx[5], pml_A, pml_B, pml_D,
          fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
          fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
          fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
          );

      break;

    default:
  }

  //
  // free surface
  //

  // boundary z2
  switch (boundary_itype[5])
  {
    case FD_BOUNDARY_TYPE_FREE :

      // get aux vars and coes
      siz_aux_volume = abs_blk_vars_siz_volume[5];

      ptr_pml = pml_cur + abs_blk_vars_blkpos[5];
      a_Vx = ptr_pml + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
      a_Vy = ptr_pml + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
      // add a_Vz etc

      ptr_pml = pml_rhs + abs_blk_vars_blkpos[5];
      a_hVx = ptr_pml + WF_EL_1ST_SEQ_PMLVX * siz_aux_volume;
      a_hVy = ptr_pml + WF_EL_1ST_SEQ_PMLVY * siz_aux_volume;
      // add a_Vz etc

      pml_A = abs_blk_coefs[5];
      pml_B = pml_A + abs_num_of_layers[5];
      pml_D = pml_B + abs_num_of_layers[5];

      // tractiong
      ierr = sv_eliso1st_curv_macdrp_rhs_timg_z2(
          xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z, lam3d, mu3d, rho3d,
          Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
          a_Txx,a_Tyy,a_Tzz,a_Txz,a_Tyz,a_Txy, a_hVx,a_hVy,a_hVz,
          ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,siz_line,siz_volume,
          abs_blk_indx[0], pml_A, pml_B, pml_D,
          fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
          fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
          fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
          );

      // velocity: vlow
      ierr = sv_eliso1st_curv_macdrp_rhs_vlow_z2();
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

int sv_eliso1st_curv_macdrp_rhs_inner(
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

int sv_eliso1st_curv_macdrp_rhs_timg_z2(
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype, int *restrict abs_num_of_layers,
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
          lam = lam3d[iptr];
          mu  =  mu3d[iptr];
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
          if (boundary_itype[0] == FD_BOUNDARY_TYPE_CFSPML && i<ni1+abs_num_of_layers[1])
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

int sv_eliso1st_curv_macdrp_rhs_vlow_z2(
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
          if (boundary_itype[0] == FD_BOUNDARY_TYPE_CFSPML && i<ni1+abs_num_of_layers[1])
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

int sv_eliso1st_curv_macdrp_rhs_cfspml_x(
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict  mu3d, float *restrict rho3d,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict a_Vx  , float *restrict a_Vy  , float *restrict a_Vz  ,
    float *restrict a_Txx , float *restrict a_Tyy , float *restrict a_Tzz ,
    float *restrict a_Txz , float *restrict a_Tyz , float *restrict a_Txy ,
    float *restrict a_hVx , float *restrict a_hVy , float *restrict a_hVz ,
    float *restrict a_hTxx, float *restrict a_hTyy, float *restrict a_hTzz,
    float *restrict a_hTxz, float *restrict a_hTyz, float *restrict a_hTxy,
    size_t siz_line, siz_t siz_volue,
    size_t *restrict abs_indx, // pml range
    float  *restrict Ax, float  *restrict Bx, float  *restrict Dx,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef
    )
{
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // local var
  size_t n_fd; // loop var for fd
  size_t i,j,k;
  size_t iptr, iptr_a;

  size_t abs_ni1, size_t abs_ni2, 
  size_t abs_nj1, size_t abs_nj2, 
  size_t abs_nk1, size_t abs_nk2, 

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

  // get index into local var
  abs_ni1 = abs_indx[0];
  abs_ni2 = abs_indx[1];
  abs_nj1 = abs_indx[2];
  abs_nj2 = abs_indx[3];
  abs_nk1 = abs_indx[4];
  abs_nk2 = abs_indx[5];

  // loop all points
  iptr_a = 0;

  for (k=abs_nk1; k<=abs_nk2; k++)
  {
    iptr_k = k * siz_slice;
     
    for (j=abs_nj1; j<=abs_nj2; j++)
    {
      iptr_j = iptr_k + j * siz_line;

      iptr = iptr_j + abs_ni1;

      for (i=abs_ni1; i<=abs_ni2; i++)
      {
        // pml coefs
        d1=Dx[i-abs_ni1];
        a1=Ax[i-abs_ni1];
        b1=Bx[i-abs_ni1];
        rb1 = 1.0 / b1;

        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        rho = rho3d[iptr];

        lam2mu = lam + 2.0 * mu;
        rrho = 1.0 / rho;

        // xi derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DxVy, Vy, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DxVz, Vz, iptr, fd_len, lfdx_shift, lfdx_coef, n_fd);

        // combine for corr and aux vars
        hVx_rhs = ( xix*DxTxx + xiy*DxTxy + xiz*DxTxz ) * rrho;

        // need to add Tij

        // 1: make corr to moment equation
        hVx[iptr] += (rb1 - 1.0) * hVx_rhs - rb1 * a_Vx[iptr_a];

        // make corr to Hooke's equatoin
        
        // 2: aux var
        //   a1 = alpha + d / beta, dealt in abs_set_cfspml
        a_hVx[iptr_a] = d1 * hVx_rhs - a1 * a_Vx[iptr_a];

        // incr index
        iptr   += 1;
        iptr_a += 1;
      }
    }
  }
}

// apply on wavefield, not on derivatives
int sv_eliso1st_curv_macdrp_apply_ablexp_pointcoef(float *restrict w_cur, 
    int number_of_vars, size_t *restrict vars_pos,
    // grid size
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    // ablexp info
    size_t *restrict abs_indx, float *restrict abs_damp
    )
{
  int ierr = 0;

  size_t iptr,iptr_abs, iptr_var, iptr_k, iptr_j;
  size_t siz_line, siz_slice;
  size_t i,j,k,ivar;

  siz_line  = nx;
  siz_slice = nx * ny;

  for (ivar=0; ivar<number_of_vars; ivar++)
  {
    iptr_var = vars_pos[ivar];
    iptr_abs = 0;

    for (k=abs_indx[4]; k<=abs_indx[5]; k++)
    {
      iptr_k = iptr_var + k * siz_slice;
      for (j=abs_indx[2]; j<=abs_indx[3]; j++)
      {
        iptr_j = j * siz_line + iptr_k;
        for (i=abs_indx[0]; i<=abs_indx[2]; i++)
        {
          iptr = i + iptr_j;

          w_cur[iptr] *= abs_damp[iptr_abs];

          // next
          iptr_abs += 1;
        }
      }
    }
  }

  return ierr;
}

int sv_eliso1st_curv_macdrp_apply_ablexp(float *restrict w_cur, 
    int number_of_vars, size_t *restrict vars_pos,
    // grid size
    size_t nx, size_t ny, size_t nz,
    // ablexp info
    size_t *restrict abs_indx, float *restrict abs_damp
    )
{
  int ierr = 0;

  size_t iptr,iptr_abs, iptr_var, iptr_k, iptr_j;
  size_t siz_line, siz_slice;
  size_t i,j,k,ivar;
  float *restrict Dx, Dy, Dz;
  float dampx, dampy, dampz, dampmin;

  siz_line  = nx;
  siz_slice = nx * ny;

  Dx = abs_damp;
  Dy = Dx + nx;
  Dz = Dy + ny; 

  for (ivar=0; ivar<number_of_vars; ivar++)
  {
    iptr_var = vars_pos[ivar];

    for (k=abs_indx[4]; k<=abs_indx[5]; k++)
    {
      iptr_k = iptr_var + k * siz_slice;
      dampz = Dz[k];

      for (j=abs_indx[2]; j<=abs_indx[3]; j++)
      {
        iptr_j = j * siz_line + iptr_k;
        dampy = Dy[j];
        dampmin = (dampy < dampz) ? dampy : dampz;

        for (i=abs_indx[0]; i<=abs_indx[2]; i++)
        {
          iptr = i + iptr_j;
          dampx = Dx[i];
          dampmin = (dampx < dampy) ? dampx : dampy;

          w_cur[iptr] *= dampmin;

        }
      }
    }
  }

  return ierr;
}

/*
int eliso1st_curv_macdrp_rhs_cfspml_x(
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict  mu3d, float *restrict rho3d,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict a_Vx  , float *restrict a_Vy  , float *restrict a_Vz  ,
    float *restrict a_Txx , float *restrict a_Tyy , float *restrict a_Tzz ,
    float *restrict a_Txz , float *restrict a_Tyz , float *restrict a_Txy ,
    float *restrict a_hVx , float *restrict a_hVy , float *restrict a_hVz ,
    float *restrict a_hTxx, float *restrict a_hTyy, float *restrict a_hTzz,
    float *restrict a_hTxz, float *restrict a_hTyz, float *restrict a_hTxy,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, siz_t siz_volue,
    size_t *restrict abs_indx, // pml range
    float  *restrict Ax, float  *restrict Bx, float  *restrict Dx,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef
    )
    */
