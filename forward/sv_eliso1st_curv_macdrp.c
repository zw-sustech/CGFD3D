/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and macdrp schem
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "fd_t.h"
#include "par_t.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "abs_funcs.h"
#include "src_funcs.h"
#include "sv_eliso1st_curv_macdrp.h"

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_allstep(
    float *restrict w3d,  // wavefield
    float *restrict g3d,  // grid vars
    float *restrict m3d,  // medium vars
    // grid size
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    // boundary type
    int *restrict boundary_itype,
    // if abs
    int              abs_itype, //
    int    *restrict abs_num_of_layers, //
    size_t *restrict abs_indx, //
    size_t *restrict abs_coefs_facepos0, //
    float  *restrict abs_coefs, //
    size_t           abs_vars_size_per_level, //
    size_t *restrict abs_vars_volsiz, //
    size_t *restrict abs_vars_facepos0, //
    float  *restrict abs_vars,
    // if free surface
    float *matVx2Vz, //
    float *matVy2Vz, //
    // source term
    int num_of_force,
    int *restrict force_info, // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
    float *restrict force_vec_stf,
    int num_of_moment,
    int *restrict moment_info, // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
    float *restrict moment_ten_rate,
    // io
    int num_of_sta, int *restrict sta_loc_point, float *restrict sta_seismo,
    int num_of_snap, int *restrict snap_grid_indx, int *restrict snap_time_indx,
    // scheme
    int num_rk_stages, float *rk_a, float *rk_b,
    int num_of_pairs, 
    size_t fdx_max_half_len, size_t fdy_max_half_len,
    size_t fdz_max_len, size_t fdz_num_surf_lay,
    size_t ****pair_fdx_all_info, size_t ***pair_fdx_all_indx, float ***pair_fdx_all_coef,
    size_t ****pair_fdy_all_info, size_t ***pair_fdy_all_indx, float ***pair_fdy_all_coef,
    size_t ****pair_fdz_all_info, size_t ***pair_fdz_all_indx, float ***pair_fdz_all_coef,
    // time
    float dt, int nt_total, float t0,
    // mpi
    int myid, int *myid2, MPI_Comm comm,
    int qc_check_nan_num_of_step,
    const int verbose, // used for fprint qc
    char *out_dir)
{
  // local allocated array
  float *force_values     = NULL;  // num_of_force * 3
  float *moment_ten_value = NULL;  // num_of_moment * 6

  // local pointer
  float *restrict w_cur;
  float *restrict w_pre;
  float *restrict w_rhs;
  float *restrict w_end;
  float *restrict w_tmp;

  float *restrict abs_vars_cur;
  float *restrict abs_vars_pre;
  float *restrict abs_vars_rhs;
  float *restrict abs_vars_end;
  float *restrict abs_vars_tmp;

  int   ipair, istage;
  float t_cur;

  size_t w3d_size_per_level = w3d_num_of_vars * siz_volume;

  // alloc
  if (num_of_force > 0) {
    force_vec_value = (float *)fdlib_mem_malloc_1d_float(num_of_force*3, "alloc force stf all step");
  }
  if (num_of_moment > 0) {
    moment_ten_value = (float *)fdlib_mem_malloc_1d_float(num_of_force*6, "alloc force stf all step");
  }

  // get wavefield
  w_pre = w3d + w3d_size_per_level * 0;
  w_rhs = w3d + w3d_size_per_level * 2;
  w_end = w3d + w3d_size_per_level * 3;

  // get pml
  abs_vars_pre = abs_vars + abs_vars_size_per_level * 0;
  abs_vars_rhs = abs_vars + abs_vars_size_per_level * 2;
  abs_vars_end = abs_vars + abs_vars_size_per_level * 3;

  //
  // time loop
  //
  for (int it=0; it<nt_total; it++)
  {
    t_cur = it * dt + t0;

    if (myid==0 && verbose>10) fprintf(stdout,"-> it=%d, t=\%f\n", it, t_cur);

    ipair = it % num_of_pair; // mod to ipair
    if (myid==0 && verbose>10) fprintf(stdout, " --> ipair=%d\n",ipair);

    // loop RK stages for one step
    for (istage=0; istage<num_rk_stages; istage++)
    {
      if (myid==0 && verbose>10) fprintf(stdout, " --> istage=%d\n",istage);

      if (istage==0) {
          w_cur   = w_pre; // to aoid one dubplication 
          abs_vars_cur = abs_vars_pre; 
      } else {
          w_cur   = w3d + w3d_size_per_level * 1;
          abs_vars_cur = abs_vars + abs_vars_size_per_level * 1;
      }

      // stf value for cur stage
      src_get_stage_stf(num_of_force,
                        force_info,
                        force_vec_stf,
                        num_of_moment,
                        moment_info,
                        moment_ten_stf,
                        it,
                        istage,
                        num_rk_stages,
                        force_vec_value,
                        moment_ten_value,
                        myid, verbose);

      /*
      // pack message
      wf_el_1st_pack_message(w_cur,
              ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,
              pair_fdx_len[ipair][istage][fd_max_op_len],
              buff, buff_dim_pos, buff_dim_size);

      // exchange data along y
      MPI_Request r_reqs[4];

      MPI_Startall(2, reqs+0, MPI_STATUS_IGNORE);

      // wait for communication
      MPI_Waitall(2, reqs+0, MPI_STATUS_IGNORE);

      // receive message and uppack
      wf_el_1st_unpack_messag();
      */

      // compute
      sv_eliso1st_curv_macdrp_onestage(w_cur, w_rhs, g3d, m3d, 
                                       ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,
                                       siz_line,siz_slice,siz_volume,
                                       boundary_itype,
                                       abs_itype,
                                       abs_num_of_layers,
                                       abs_indx,
                                       abs_coefs_facepos0,
                                       abs_coefs,
                                       abs_vars_size_per_level,
                                       abs_vars_volsiz,
                                       abs_vars_facepos0,
                                       abs_vars_cur,
                                       abs_vars_rhs,
                                       matVx2Vz, matVy2Vz,
                                       num_of_force , force_loc_point , force_vec_value,
                                       num_of_moment, moment_loc_point, moment_ten_value,
                                       fdx_max_half_len, fdy_max_half_len,
                                       fdz_max_len, fdz_num_surf_lay,
                                       pair_fdx_all_info[ipair][istage],
                                       pair_fdx_all_indx[ipair][istage],
                                       pair_fdx_all_coef[ipair][istage],
                                       pair_fdy_all_info[ipair][istage],
                                       pair_fdy_all_indx[ipair][istage],
                                       pair_fdy_all_coef[ipair][istage],
                                       pair_fdz_all_info[ipair][istage],
                                       pair_fdz_all_indx[ipair][istage],
                                       pair_fdz_all_coef[ipair][istage],
                                       myid, verbose);

      // RK int

      // cal w_cur first, then start send etc; or w_end at last stage
      if (istage != num_rk_stages-1)
      {
        float coef_a = rk_a[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
            w_cur[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
        }
        // apply abs
        /*
        if (abs_has_ablexp) {
          sv_eliso1st_curv_macdrp_ablexp(w_cur, w3d_num_of_vars,
                    w3d_pos, nx, ny, nz, abs_coefs_facepos0, abs_coefs);
        }
        */
        // pack and isend
        //pack_message(w_cur,sbuff);
        //MPI_Startall(sreqs);

        // apply pml
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_cur[iptr] = abs_vars_pre[iptr] + coef_a * abs_vars_rhs[iptr];
          }
        }
      }
      else // w_end for send at last stage
      {
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }
        /*
        if (abs_has_ablexp) {
          sv_eliso1st_curv_macdrp_ablexp(w_end, w3d_num_of_vars,
                    w3d_pos, nx, ny, nz, abs_coefs_facepos0, abs_coefs);
        }
        */
        // pack and isend
        //pack_message(w_end,sbuff);
        //MPI_Startall(sreqs);

        // apply pml
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_end[iptr] += coef_b * abs_vars_rhs[iptr];
          }
        }
      }

      // cal w_end at other stages, non-communicate vars
      if (istage==0)
      {
        float coef_b = rk_b[istage] * dh;

        // wavefield
        for (size_t iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
            w_end[iptr] = w_pre[iptr] + coef_b * w_rhs[iptr];
        }

        // apply pml
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_end[iptr] = abs_vars_pre[iptr] + coef_b * abs_vars_rhs[iptr];
          }
        }
      }
      else if (istage<num_rk_stages-1)
      {
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }

        // apply pml
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_end[iptr] += coef_b * abs_vars_rhs[iptr];
          }
        }
      } 
    } // RK stages

    // QC
    if (qc_check_nan_num_of_step >0  && (it % qc_check_nan_num_of_step) == 0) {
      if (myid==0 && verbose>10) fprintf(stdout,"-> check value nan\n");
        wf_el_1st_check_value(w_end);
    }

    // swap w_pre and w_end, avoid copying
    w_tmp = w_pre; w_pre = w_end; w_end = w_tmp;
    abs_vars_tmp = abs_vars_pre; abs_vars_pre = abs_vars_end; abs_vars_end = abs_vars_tmp;

    // save results
    //io_seismo_keep();
    //io_snapshot_save();

  } // time loop

  // postproc
  if (force_vec_value ) free(force_vec_value );
  if (moment_ten_value) free(moment_ten_value);
}

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_onestage(
    float *restrict w_cur, float *restrict rhs, 
    float *restrict g3d, float *restrict m3d,
    // grid size
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t ni, size_t nj, size_t nk, size_t nx, size_t ny, size_t nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    int *restrict boundary_itype,
    int              abs_itype,
    int    *restrict abs_num_of_layers,
    size_t *restrict abs_indx,
    size_t *restrict abs_coefs_facepos0,
    float  *restrict abs_coefs,
    size_t           abs_vars_size_per_level,
    size_t *restrict abs_vars_volsiz,
    size_t *restrict abs_vars_facepos0, //
    float  *restrict abs_vars_cur,
    float  *restrict abs_vars_rhs,
    float *matVx2Vz, float *matVy2Vz, //
    // source term
    int num_of_force,
    int *restrict force_loc_point,
    float *restrict force_vec_value, // only for cur stage, size: num_of_force
    int num_of_moment,
    int *restrict moment_loc_point,
    float *restrict moment_ten_value,
    // include different order/stentil
    int fdx_max_half_len, int fdy_max_half_len,
    int fdz_max_len, int fdz_num_surf_lay,
    size_t **restrict fdx_all_info, size_t *restrict fdx_all_indx,float *restrict fdx_all_coef,
    size_t **restrict fdy_all_info, size_t *restrict fdy_all_indx,float *restrict fdy_all_coef,
    size_t **restrict fdz_all_info, size_t *restrict fdz_all_indx,float *restrict fdz_all_coef,
    const int myid, const int verbose)
{
  size_t i,j,k;

  // local pointer get each vars
  float *restrict Vx    = w_cur + WF_EL_1ST_SEQ_VX    * siz_volume;
  float *restrict Vy    = w_cur + WF_EL_1ST_SEQ_VY    * siz_volume;
  float *restrict Vz    = w_cur + WF_EL_1ST_SEQ_VZ    * siz_volume;
  float *restrict Txx   = w_cur + WF_EL_1ST_SEQ_TXX   * siz_volume;
  float *restrict Tyy   = w_cur + WF_EL_1ST_SEQ_TYY   * siz_volume;
  float *restrict Tzz   = w_cur + WF_EL_1ST_SEQ_TZZ   * siz_volume;
  float *restrict Txz   = w_cur + WF_EL_1ST_SEQ_TXZ   * siz_volume;
  float *restrict Tyz   = w_cur + WF_EL_1ST_SEQ_TYZ   * siz_volume;
  float *restrict Txy   = w_cur + WF_EL_1ST_SEQ_TXY   * siz_volume;
  float *restrict hVx   = w_rhs + WF_EL_1ST_SEQ_VX    * siz_volume;
  float *restrict hVy   = w_rhs + WF_EL_1ST_SEQ_VY    * siz_volume;
  float *restrict hVz   = w_rhs + WF_EL_1ST_SEQ_VZ    * siz_volume;
  float *restrict hTxx  = w_rhs + WF_EL_1ST_SEQ_TXX   * siz_volume;
  float *restrict hTyy  = w_rhs + WF_EL_1ST_SEQ_TYY   * siz_volume;
  float *restrict hTzz  = w_rhs + WF_EL_1ST_SEQ_TZZ   * siz_volume;
  float *restrict hTxz  = w_rhs + WF_EL_1ST_SEQ_TXZ   * siz_volume;
  float *restrict hTyz  = w_rhs + WF_EL_1ST_SEQ_TYZ   * siz_volume;
  float *restrict hTxy  = w_rhs + WF_EL_1ST_SEQ_TXY   * siz_volume;
  float *restrict xi_x  = g3d   + GD_CURV_SEQ_XIX   * siz_volume;
  float *restrict xi_y  = g3d   + GD_CURV_SEQ_XIY   * siz_volume;
  float *restrict xi_z  = g3d   + GD_CURV_SEQ_XIZ   * siz_volume;
  float *restrict et_x  = g3d   + GD_CURV_SEQ_ETX   * siz_volume;
  float *restrict et_y  = g3d   + GD_CURV_SEQ_ETY   * siz_volume;
  float *restrict et_z  = g3d   + GD_CURV_SEQ_ETZ   * siz_volume;
  float *restrict zt_x  = g3d   + GD_CURV_SEQ_ZTX   * siz_volume;
  float *restrict zt_y  = g3d   + GD_CURV_SEQ_ZTY   * siz_volume;
  float *restrict zt_z  = g3d   + GD_CURV_SEQ_ZTZ   * siz_volume;
  float *restrict jac3d = g3d   + GD_CURV_SEQ_JAC   * siz_volume;
  float *restrict lam3d = m3d   + MD_ELISO_SEQ_LAMBDA * siz_volume;
  float *restrict  mu3d = m3d   + MD_ELISO_SEQ_MU     * siz_volume;
  float *restrict slw3d = m3d   + MD_ELISO_SEQ_RHO    * siz_volume;

  // fd op
  size_t           fdx_inn_pos;
  size_t *restrict fdx_inn_indx;
  float  *restrict fdx_inn_coef;
  size_t           fdy_inn_pos;
  size_t *restrict fdy_inn_indx;
  float  *restrict fdy_inn_coef;
  size_t           fdz_inn_pos;
  size_t *restrict fdz_inn_indx;
  float  *restrict fdz_inn_coef;

  // for get a op from 1d array, currently use fdz_num_surf_lay as index
  fdx_inn_pos  = fdx_all_info[fdx_max_half_len][FD_INFO_POS_OF_INDX];
  // length, index, coef of a op
  fdx_inn_len  = fdx_all_info[fdx_max_half_len][FD_INFO_LENGTH_TOTAL];
  fdx_inn_indx = fdx_all_indx + fdx_inn_pos;
  fdx_inn_coef = fdx_all_coef + fdx_inn_pos;

  fdy_inn_pos  = fdy_all_info[fdy_max_half_len][FD_INFO_POS_OF_INDX];
  fdy_inn_len  = fdy_all_info[fdy_max_half_len][FD_INFO_LENGTH_TOTAL];
  fdy_inn_indx = fdy_all_indx + fdy_inn_pos;
  fdy_inn_coef = fdy_all_coef + fdy_inn_pos;

  fdz_inn_pos  = fdz_all_info[fdz_num_surf_lay][FD_INFO_POS_OF_INDX];
  fdz_inn_len  = fdz_all_info[fdz_num_surf_lay][FD_INFO_LENGTH_TOTAL];
  fdz_inn_indx = fdz_all_indx + fdz_inn_pos;
  fdz_inn_coef = fdz_all_coef + fdz_inn_pos;

  // inner points
  sv_eliso1st_curv_macdrp_rhs_inner(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                    hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                    lam3d, mu3d, slw3d,
                                    ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                    fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                    fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                    fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                    myid, verbose);

  // free, abs, source in turn

  // free surface at z2
  if (boundary_itype[5] == FD_BOUNDARY_TYPE_FREE)
  {
    // tractiong
    sv_eliso1st_curv_macdrp_rhs_timg_z2(Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        jac3d, slw3d,
                                        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        myid, verbose);

    // velocity: vlow
    sv_eliso1st_curv_macdrp_rhs_vlow_z2(Vx,Vy,Vz,Txx,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        lam3d, mu3d, slw3d,matVx2Vz,matVy2Vz,
                                        ni1,ni2,nj1,nj2,nk1,nk2,ni,nj,nk,nx,ny,nz,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_num_surf_lay,fdz_max_len,
                                        fdz_all_info,fdz_all_indx,fdz_all_coef,
                                        myid, verbose);
  }

  // cfs-pml, loop face inside
  if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
  {
    sv_eliso1st_curv_macdrp_rhs_cfspml(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                       xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                       lam3d, mu3d, slw3d,
                                       siz_line,siz_slice,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       boundary_itype, abs_num_of_layers, abs_indx,
                                       abs_coefs_facepos0, abs_coefs,
                                       abs_vars_volsiz, abs_vars_facepos0,
                                       abs_vars_cur, abs_vars_rhs,
                                       myid, verbose);
    
    // if free surface,  modified 
    if (boundary_itype[5] == FD_BOUNDARY_TYPE_FREE)
    {
      // cfs-pml modified for timg
      sv_eliso1st_curv_macdrp_rhs_cfspml_timg_z2(Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
                                                 xi_x, xi_y, xi_z,
                                                 et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                                 jac3d, slw3d, nk2, siz_line,siz_slice,
                                                 fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                                 fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                                 fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                                 boundary_itype, abs_num_of_layers, abs_indx,
                                                 abs_coefs_facepos0, abs_coefs,
                                                 abs_vars_volsiz, abs_vars_facepos0,
                                                 abs_vars_cur, abs_vars_rhs,
                                                 myid, verbose);
      
      // cfs-pml modified for vfree; vlow at points below surface doesn't affect pml
      sv_eliso1st_curv_macdrp_rhs_cfspml_vfree_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                                  xi_x, xi_y, xi_z,
                                                  et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                                  lam3d, mu3d, slw3d, matVx2Vz, matVy2Vz,
                                                  nk2,siz_line,siz_slice,
                                                  fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                                  fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                                  boundary_itype, abs_num_of_layers, abs_indx,
                                                  abs_coefs_facepos0, abs_coefs,
                                                  abs_vars_volsiz, abs_vars_facepos0,
                                                  abs_vars_cur, abs_vars_rhs,
                                                  myid, verbose);
    }
  }

  // source term
  /*
  switch (source_itype) {
    case SOURCE_TYPE_POINT :
      curv_macdrp_eliso1st_src_point();
      break;

    case SOURCE_TYPE_GAUSS :
      curv_macdrp_eliso1st_src_gauss();
      break;

    case SOURCE_TYPE_SINC :
      curv_macdrp_eliso1st_src_sinc();
      break;
  }
  */
  if (num_of_force>0 || num_of_moment>0)
  {
    sv_eliso1st_curv_macdrp_rhs_src(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    jac3d, slw3d, siz_line,siz_slice,
                                    num_of_force, force_loc_point, force_vec_value,
                                    num_of_moment, moment_loc_point, moment_ten_value,
                                    myid, verbose);
  }
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_rhs_inner(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // loop var for fd
  size_t n_fd; // loop var for fd

  // local var
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
#ifdef DEBUG
    if (myid==0 && verbose>90) fprintf(stdout,"----> nk1=%d,nk2=%d,k=%d\n", nk1,nk2,k);
#endif
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
#ifdef DEBUG
      if (myid==0 && verbose>91) fprintf(stdout,"-----> nj1=%d,nj2=%d,j=%d\n", nj1,nj2,j);
#endif
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
#ifdef DEBUG
        if (myid==0 && verbose>92) fprintf(stdout,"-----> ni1=%d,ni2=%d,i=%d\n", ni1,ni2,i);
#endif
        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVx, Vx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzVx, Vx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT(DxVy, Vy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVy, Vy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzVy, Vy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVz, Vz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzVz, Vz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txx derivatives
        M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyTxx, Txx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzTxx, Txx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyy derivatives
        M_FD_SHIFT(DxTyy, Tyy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyTyy, Tyy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzTyy, Tyy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tzz derivatives
        M_FD_SHIFT(DxTzz, Tzz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyTzz, Tzz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyTxz, Txz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyz derivatives
        M_FD_SHIFT(DxTyz, Tyz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyTyz, Tyz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzTyz, Tyz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txy derivatives
        M_FD_SHIFT(DxTxy, Txy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyTxy, Txy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT(DzTxy, Txy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];
        etx = et_x[iptr];
        ety = et_y[iptr];
        etz = et_z[iptr];
        ztx = zt_x[iptr];
        zty = zt_y[iptr];
        ztz = zt_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

        // moment equation
        hVx[iptr] = slw*( xix*DxTxx + xiy*DxTxy + xiz*DxTxz  
                         +etx*DyTxx + ety*DyTxy + etz*DyTxz 
                         +ztx*DzTxx + zty*DzTxy + ztz*DzTxz );
        hVy[iptr] = slw*( xix*DxTxy + xiy*DxTyy + xiz*DxTyz
                         +etx*DyTxy + ety*DyTyy + etz*DyTyz
                         +ztx*DzTxy + zty*DzTyy + ztz*DzTyz );
        hVz[iptr] = slw*( xix*DxTxz + xiy*DxTyz + xiz*DxTzz 
                         +etx*DyTxz + ety*DyTyz + etz*DyTzz
                         +ztx*DzTxz + zty*DzTyz + ztz*DzTzz );

        // Hooke's equatoin
        hTxx[iptr] =  lam2mu * ( xix*DxVx  +etx*DyVx + ztx*DzVx);
                    + lam    * ( xiy*DxVy + ety*DyVy + zty*DzVy
                                +xiz*DxVz + etz*DyVz + ztz*DzVz);

        hTyy[iptr] = lam2mu * ( xiy*DxVy + ety*DyVy + zty*DzVy)
                    +lam    * ( xix*DxVx + etx*DyVx + ztx*DzVx
                               +xiz*DxVz + etz*DyVz + ztz*DzVz);

        hTzz[iptr] = lam2mu * ( xiz*DxVz + etz*DyVz + ztz*DzVz)
                    +lam    * ( xix*DxVx  +etx*DyVx  +ztx*DzVx
                               +xiy*DxVy + ety*DyVy + zty*DzVy);

        hTxy[iptr] = mu *(
                     xiy*DxVx + xix*DxVy
                    +ety*DyVx + etx*DyVy
                    +zty*DzVx + ztx*DzVy
                    );
        hTxz[iptr] = mu *(
                     xiz*DxVx + xix*DxVz
                    +etz*DyVx + etx*DyVz
                    +ztz*DzVx + ztx*DzVz
                    );
        hTyz[iptr] = mu *(
                     xiz*DxVy + xiy*DxVz
                    +etz*DyVy + ety*DyVz
                    +ztz*DzVy + zty*DzVz
                    );

        iptr += 1;
      }
    }
  }
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_eliso1st_curv_macdrp_rhs_timg_z2(
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // loop var for fd
  size_t n_fd; // loop var for fd

  // local var
  float DxTx,DyTy,DzTz;
  float slwjac;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // to save traction and other two dir force var
  float vecxi[fdx_len];
  float vecet[fdy_len];
  float veczt[fdz_len];
  size_t n, iptr4vec;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // last indx, free surface force Tx/Ty/Tz to 0 in cal
  size_t k_min = nk2 - fdz_indx[fdz_len-1];

  // point affected by timg
  for (size_t k=k_min; k <= nk2; k++)
  {
    // k corresponding to 0 index of the fd op

    // index of free surface
    size_t n_free = nk2 - k - fdz_indx[0]; // first indx is negative

    // for 1d index
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<ni2; i++)
      {
          // metric
          xix = xi_x[iptr];
          xiy = xi_y[iptr];
          xiz = xi_z[iptr];
          etx = et_x[iptr];
          ety = et_y[iptr];
          etz = et_z[iptr];
          ztx = zt_x[iptr];
          zty = zt_y[iptr];
          ztz = zt_z[iptr];

          // slowness and jac
          slwjac = slw3d[iptr] * jac3d[iptr];

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
          for (n=0; n<fdy_len; n++) {
            iptr4vec = iptr + fdy_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txx[iptr4vec]
                                          + et_y[iptr4vec] * Txy[iptr4vec]
                                          + et_z[iptr4vec] * Txz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                          + zt_y[iptr4vec] * Txy[iptr4vec]
                                          + zt_z[iptr4vec] * Txz[iptr4vec] );
          }

          // at surface -> set to 0
          vecet[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVx[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;

          //
          // for hVy
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txy[iptr4vec]
                                          + xi_y[iptr4vec] * Tyy[iptr4vec]
                                          + xi_z[iptr4vec] * Tyz[iptr4vec] );
          }
          for (n=0; n<fdy_len; n++) {
            iptr4vec = iptr + fdy_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txy[iptr4vec]
                                          + et_y[iptr4vec] * Tyy[iptr4vec]
                                          + et_z[iptr4vec] * Tyz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txy[iptr4vec]
                                          + zt_y[iptr4vec] * Tyy[iptr4vec]
                                          + zt_z[iptr4vec] * Tyz[iptr4vec] );
          }

          // at surface -> set to 0
          vecet[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVy[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;

          //
          // for hVz
          //

          // transform to conservative vars
          for (n=0; n<fdx_len; n++) {
            iptr4vec = iptr + fdx_indx[n];
            vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txz[iptr4vec]
                                          + xi_y[iptr4vec] * Tyz[iptr4vec]
                                          + xi_z[iptr4vec] * Tzz[iptr4vec] );
          }
          for (n=0; n<fdy_len; n++) {
            iptr4vec = iptr + fdy_indx[n] * siz_line;
            vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txz[iptr4vec]
                                          + et_y[iptr4vec] * Tyz[iptr4vec]
                                          + et_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // blow surface -> cal
          for (n=0; n<n_free; n++) {
            iptr4vec = iptr + fdz_indx[n] * siz_slice;
            veczt[n] = jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                          + zt_y[iptr4vec] * Tyz[iptr4vec]
                                          + zt_z[iptr4vec] * Tzz[iptr4vec] );
          }

          // at surface -> set to 0
          vecet[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // for hVx 
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVz[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;

          // next
          iptr += 1;
      }
    }
  }
}

/*
 * implement vlow boundary
 */

void
sv_eliso1st_curv_macdrp_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_num_surf_lay, size_t fdz_max_len, size_t **restrict fdz_all_info,
    size_t *restrict fdz_all_indx, float *restrict fdz_all_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];

  // allocate max_len because fdz may have different lens
  float  lfdz_coef [fdz_max_len];
  size_t lfdz_shift[fdz_max_len];

  // local var
  size_t i,j,k;
  size_t n_fd; // loop var for fd
  size_t fdz_len;

  // local var
  float DxVx,DxVy,DxVz;
  float DyVx,DyVy,DyVz;
  float DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }

  // loop near surface layers
  for (size_t n=0; n < fdz_num_surf_lay; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // get pos and len for this point
    size_t  lfdz_pos0 = fdz_all_info[n][FD_INFO_POS_OF_INDX];
    size_t  lfdz_len  = fdz_all_info[n][FD_INFO_LENGTH_TOTAL];
    // point to indx/coef for this point
    size_t  *p_fdz_indx  = fdz_all_indx[lfdz_pos0];
    float   *p_fdz_coef  = fdz_all_coef[lfdz_pos0];
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_shift[n_fd] = p_fdz_indx[n_fd] * siz_slice;
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd];
    }

    // for index
    size_t iptr_k = k * siz_slice;

    for (j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (i=ni1; i<=ni2; i++)
      {
        // metric
        xix = xi_x[iptr];
        xiy = xi_y[iptr];
        xiz = xi_z[iptr];
        etx = et_x[iptr];
        ety = et_y[iptr];
        etz = et_z[iptr];
        ztx = zt_x[iptr];
        zty = zt_y[iptr];
        ztz = zt_z[iptr];

        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVx, Vx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT(DxVy, Vy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVy, Vy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVz, Vz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        if (k==nk2) // at surface, convert
        {
          size_t ij = i + j * siz_line;
          DzVx = matVx2Vz[ij+3*0+0] * DxVx
               + matVx2Vz[ij+3*0+1] * DxVy
               + matVx2Vz[ij+3*0+2] * DxVz
               + matVy2Vz[ij+3*0+0] * DyVx
               + matVy2Vz[ij+3*0+1] * DyVy
               + matVy2Vz[ij+3*0+2] * DyVz;

          DzVy = matVx2Vz[ij+3*1+0] * DxVx
               + matVx2Vz[ij+3*1+1] * DxVy
               + matVx2Vz[ij+3*1+2] * DxVz
               + matVy2Vz[ij+3*1+0] * DyVx
               + matVy2Vz[ij+3*1+1] * DyVy
               + matVy2Vz[ij+3*1+2] * DyVz;

          DzVz = matVx2Vz[ij+3*2+0] * DxVx
               + matVx2Vz[ij+3*2+1] * DxVy
               + matVx2Vz[ij+3*2+2] * DxVz
               + matVy2Vz[ij+3*2+0] * DyVx
               + matVy2Vz[ij+3*2+1] * DyVy
               + matVy2Vz[ij+3*2+2] * DyVz;
        }
        else // lower than surface, lower order
        {
          M_FD_SHIFT(DzVx, Vx, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVy, Vy, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVz, Vz, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
        }

        // Hooke's equatoin
        hTxx[iptr] =  lam2mu * ( xix*DxVx  +etx*DyVx + ztx*DzVx);
                    + lam    * ( xiy*DxVy + ety*DyVy + zty*DzVy
                                +xiz*DxVz + etz*DyVz + ztz*DzVz);

        hTyy[iptr] = lam2mu * ( xiy*DxVy + ety*DyVy + zty*DzVy)
                    +lam    * ( xix*DxVx + etx*DyVx + ztx*DzVx
                               +xiz*DxVz + etz*DyVz + ztz*DzVz);

        hTzz[iptr] = lam2mu * ( xiz*DxVz + etz*DyVz + ztz*DzVz)
                    +lam    * ( xix*DxVx  +etx*DyVx  +ztx*DzVx
                               +xiy*DxVy + ety*DyVy + zty*DzVy);

        hTxy[iptr] = mu *(
                     xiy*DxVx + xix*DxVy
                    +ety*DyVx + etx*DyVy
                    +zty*DzVx + ztx*DzVy
                    );
        hTxz[iptr] = mu *(
                     xiz*DxVx + xix*DxVz
                    +etz*DyVx + etx*DyVz
                    +ztz*DzVx + ztx*DzVz
                    );
        hTyz[iptr] = mu *(
                     xiz*DxVy + xiy*DxVz
                    +etz*DyVy + ety*DyVz
                    +ztz*DzVy + zty*DzVz
                    );

        iptr += 1;
      }
    }
  }
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

void
sv_eliso1st_curv_macdrp_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose)
{
  // loop var for fd
  size_t n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // val on point
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
  float hVx_rhs,hVy_rhs,hVz_rhs;
  float hTxx_rhs,hTyy_rhs,hTzz_rhs,hTxz_rhs,hTyz_rhs,hTxy_rhs

  // local
  size_t i,j,k;
  size_t iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_B_minus_1;

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

  // loop each face
  for (int iface=0; iface<FD_NDIM_2; iface++)
  {
    // skip to next face if not cfspml
    if (boundary_itype[iface] != FD_BOUNDARY_TYPE_CFSPML) continue;

    // get index into local var
    size_t abs_ni1 = abs_indx[iface+FD_NDIM_2+0];
    size_t abs_ni2 = abs_indx[iface+FD_NDIM_2+1];
    size_t abs_nj1 = abs_indx[iface+FD_NDIM_2+2];
    size_t abs_nj2 = abs_indx[iface+FD_NDIM_2+3];
    size_t abs_nk1 = abs_indx[iface+FD_NDIM_2+4];
    size_t abs_nk2 = abs_indx[iface+FD_NDIM_2+5];

    // get coef for this face
    float *restrict ptr_ceof_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    size_t pos_cur_face = abs_vars_facepos0[iface];
    size_t pml_volsiz   = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;
    float *restrict pml_Txx = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_TXX * pml_volsiz;
    float *restrict pml_Tyy = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_TYY * pml_volsiz;
    float *restrict pml_Tzz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_TZZ * pml_volsiz;
    float *restrict pml_Txz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_TXZ * pml_volsiz;
    float *restrict pml_Tyz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_TYZ * pml_volsiz;
    float *restrict pml_Txy = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_TXY * pml_volsiz;
    float *restrict pml_hVx  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_hVy  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_hVz  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;
    float *restrict pml_hTxx = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TXX * pml_volsiz;
    float *restrict pml_hTyy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TYY * pml_volsiz;
    float *restrict pml_hTzz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TZZ * pml_volsiz;
    float *restrict pml_hTxz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TXZ * pml_volsiz;
    float *restrict pml_hTyz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TYZ * pml_volsiz;
    float *restrict pml_hTxy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TXY * pml_volsiz;

    // for each dim
    if (iface < 2) // x direction
    {
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
            int abs_i = i - abs_ni1;
            coef_D = ptr_coef_D[abs_i];
            coef_A = ptr_coef_A[abs_i];
            coef_B = ptr_coef_B[abs_i];
            coef_B_minus_1 = coef_B - 1.0;

            // metric
            xix = xi_x[iptr];
            xiy = xi_y[iptr];
            xiz = xi_z[iptr];

            // medium
            lam = lam3d[iptr];
            mu  =  mu3d[iptr];
            slw = slw3d[iptr];
            lam2mu = lam + 2.0 * mu;

            // xi derivatives
            M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxVy , Vy , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxTxx, Txx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxTyy, Tyy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxTzz, Tzz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxTxz, Txz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxTyz, Tyz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
            M_FD_SHIFT(DxTxy, Txy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

            // combine for corr and aux vars
             hVx_rhs = slw * ( xix*DxTxx + xiy*DxTxy + xiz*DxTxz );
             hVy_rhs = slw * ( xix*DxTxy + xiy*DxTyy + xiz*DxTyz );
             hVz_rhs = slw * ( xix*DxTxz + xiy*DxTyz + xiz*DxTzz );
            hTxx_rhs = lam2mu*xix*DxVx + lam*xiy*DxVy + lam*xiz*DxVz;
            hTyy_rhs = lam*xix*DxVx + lam2mu*xiy*DxVy + lam*xiz*DxVz;
            hTzz_rhs = lam*xix*DxVx + lam*xiy*DxVy + lam2mu*xiz*DxVz;
            hTxy_rhs = mu*( xiy*DxVx + xix*DxVy );
            hTxz_rhs = mu*( xiz*DxVx + xix*DxVz );
            hTyz_rhs = mu*( xiz*DxVy + xiy*DxVz );

            // 1: make corr to moment equation
            hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
            hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
            hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

            // make corr to Hooke's equatoin
            hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
            hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
            hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
            hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
            hTyz[iptr] += coef_B_minus_1 * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
            hTxy[iptr] += coef_B_minus_1 * hTxy_rhs - coef_B * pml_Txy[iptr_a];
            
            // 2: aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVx[iptr_a] = d1 * hVx_rhs - a1 * pml_Vx[iptr_a];
            pml_hVy[iptr_a] = d1 * hVy_rhs - a1 * pml_Vy[iptr_a];
            pml_hVz[iptr_a] = d1 * hVz_rhs - a1 * pml_Vz[iptr_a];
            pml_hTxx[iptr_a] = d1 * hTxx_rhs - a1 * pml_Txx[iptr_a];
            pml_hTyy[iptr_a] = d1 * hTyy_rhs - a1 * pml_Tyy[iptr_a];
            pml_hTzz[iptr_a] = d1 * hTzz_rhs - a1 * pml_Tzz[iptr_a];
            pml_hTxz[iptr_a] = d1 * hTxz_rhs - a1 * pml_Txz[iptr_a];
            pml_hTyz[iptr_a] = d1 * hTyz_rhs - a1 * pml_Tyz[iptr_a];
            pml_hTxy[iptr_a] = d1 * hTxy_rhs - a1 * pml_Txy[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
    else if (iface < 4) // y direction
    {
      iptr_a = 0;
      for (k=abs_nk1; k<=abs_nk2; k++)
      {
        iptr_k = k * siz_slice;
        for (j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;

          // pml coefs
          int abs_j = j - abs_nj1;
          coef_D = ptr_coef_D[abs_j];
          coef_A = ptr_coef_A[abs_j];
          coef_B = ptr_coef_B[abs_j];
          coef_B_minus_1 = coef_B - 1.0;

          for (i=abs_ni1; i<=abs_ni2; i++)
          {
            // metric
            etx = et_x[iptr];
            ety = et_y[iptr];
            etz = et_z[iptr];

            // medium
            lam = lam3d[iptr];
            mu  =  mu3d[iptr];
            slw = slw3d[iptr];
            lam2mu = lam + 2.0 * mu;

            // et derivatives
            M_FD_SHIFT(DyVx , Vx , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyVy , Vy , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyVz , Vz , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyTxx, Txx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyTyy, Tyy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyTzz, Tzz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyTxz, Txz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyTyz, Tyz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
            M_FD_SHIFT(DyTxy, Txy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

            // combine for corr and aux vars
             hVx_rhs = slw * ( etx*DyTxx + ety*DyTxy + etz*DyTxz );
             hVy_rhs = slw * ( etx*DyTxy + ety*DyTyy + etz*DyTyz );
             hVz_rhs = slw * ( etx*DyTxz + ety*DyTyz + etz*DyTzz );
            hTxx_rhs = lam2mu*etx*DyVx + lam*ety*DyVy + lam*etz*DyVz;
            hTyy_rhs = lam*etx*DyVx + lam2mu*ety*DyVy + lam*etz*DyVz;
            hTzz_rhs = lam*etx*DyVx + lam*ety*DyVy + lam2mu*etz*DyVz;
            hTxy_rhs = mu*( ety*DyVx + etx*DyVy );
            hTxz_rhs = mu*( etz*DyVx + etx*DyVz );
            hTyz_rhs = mu*( etz*DyVy + ety*DyVz );

            // 1: make corr to moment equation
            hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
            hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
            hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

            // make corr to Hooke's equatoin
            hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
            hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
            hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
            hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
            hTyz[iptr] += coef_B_minus_1 * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
            hTxy[iptr] += coef_B_minus_1 * hTxy_rhs - coef_B * pml_Txy[iptr_a];
            
            // 2: aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVx[iptr_a] = d1 * hVx_rhs - a1 * pml_Vx[iptr_a];
            pml_hVy[iptr_a] = d1 * hVy_rhs - a1 * pml_Vy[iptr_a];
            pml_hVz[iptr_a] = d1 * hVz_rhs - a1 * pml_Vz[iptr_a];
            pml_hTxx[iptr_a] = d1 * hTxx_rhs - a1 * pml_Txx[iptr_a];
            pml_hTyy[iptr_a] = d1 * hTyy_rhs - a1 * pml_Tyy[iptr_a];
            pml_hTzz[iptr_a] = d1 * hTzz_rhs - a1 * pml_Tzz[iptr_a];
            pml_hTxz[iptr_a] = d1 * hTxz_rhs - a1 * pml_Txz[iptr_a];
            pml_hTyz[iptr_a] = d1 * hTyz_rhs - a1 * pml_Tyz[iptr_a];
            pml_hTxy[iptr_a] = d1 * hTxy_rhs - a1 * pml_Txy[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
    else // z direction
    {
      iptr_a = 0;
      for (k=abs_nk1; k<=abs_nk2; k++)
      {
        iptr_k = k * siz_slice;

        // pml coefs
        int abs_k = k - abs_nk1;
        coef_D = ptr_coef_D[abs_k];
        coef_A = ptr_coef_A[abs_k];
        coef_B = ptr_coef_B[abs_k];
        coef_B_minus_1 = coef_B - 1.0;

        for (j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;
          for (i=abs_ni1; i<=abs_ni2; i++)
          {
            // metric
            ztx = zt_x[iptr];
            zty = zt_y[iptr];
            ztz = zt_z[iptr];

            // medium
            lam = lam3d[iptr];
            mu  =  mu3d[iptr];
            slw = slw3d[iptr];
            lam2mu = lam + 2.0 * mu;

            // zt derivatives
            M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzVy , Vy , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzTxx, Txx, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzTyy, Tyy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzTzz, Tzz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzTxz, Txz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzTyz, Tyz, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
            M_FD_SHIFT(DzTxy, Txy, iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

            // combine for corr and aux vars
             hVx_rhs = slw * ( ztx*DzTxx + zty*DzTxy + ztz*DzTxz );
             hVy_rhs = slw * ( ztx*DzTxy + zty*DzTyy + ztz*DzTyz );
             hVz_rhs = slw * ( ztx*DzTxz + zty*DzTyz + ztz*DzTzz );
            hTxx_rhs = lam2mu*ztx*DzVx + lam*zty*DzVy + lam*ztz*DzVz;
            hTyy_rhs = lam*ztx*DzVx + lam2mu*zty*DzVy + lam*ztz*DzVz;
            hTzz_rhs = lam*ztx*DzVx + lam*zty*DzVy + lam2mu*ztz*DzVz;
            hTxy_rhs = mu*( zty*DzVx + ztx*DzVy );
            hTxz_rhs = mu*( ztz*DzVx + ztx*DzVz );
            hTyz_rhs = mu*( ztz*DzVy + zty*DzVz );

            // 1: make corr to moment equation
            hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
            hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
            hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

            // make corr to Hooke's equatoin
            hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
            hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
            hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
            hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
            hTyz[iptr] += coef_B_minus_1 * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
            hTxy[iptr] += coef_B_minus_1 * hTxy_rhs - coef_B * pml_Txy[iptr_a];
            
            // 2: aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVx[iptr_a] = d1 * hVx_rhs - a1 * pml_Vx[iptr_a];
            pml_hVy[iptr_a] = d1 * hVy_rhs - a1 * pml_Vy[iptr_a];
            pml_hVz[iptr_a] = d1 * hVz_rhs - a1 * pml_Vz[iptr_a];
            pml_hTxx[iptr_a] = d1 * hTxx_rhs - a1 * pml_Txx[iptr_a];
            pml_hTyy[iptr_a] = d1 * hTyy_rhs - a1 * pml_Tyy[iptr_a];
            pml_hTzz[iptr_a] = d1 * hTzz_rhs - a1 * pml_Tzz[iptr_a];
            pml_hTxz[iptr_a] = d1 * hTxz_rhs - a1 * pml_Txz[iptr_a];
            pml_hTyz[iptr_a] = d1 * hTyz_rhs - a1 * pml_Tyz[iptr_a];
            pml_hTxy[iptr_a] = d1 * hTxy_rhs - a1 * pml_Txy[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
  }
}

// considering timg free surface in cfs-pml
void
sv_eliso1st_curv_macdrp_rhs_cfspml_timg_z2(
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    size_t nk2, size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_len, size_t *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose)
{
  // loop var for fd
  size_t n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];
  
  // val on point
  float slw, slwjac;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // local
  float DxCx, DyCy;
  float coef_A, coef_B, coef_D;
  float hVx_rhs,hVy_rhs,hVz_rhs;

  // to save traction and other two dir force var
  float vecxi[fdz_len];
  float vecet[fdz_len];
  float veczt[fdz_len];
  size_t n, iptr4vec;
  size_t iptr, iptr_j, iptr_k, iptr_a;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // last indx, free surface force Tx/Ty/Tz to 0 in cal
  size_t k_min = nk2 - fdz_indx[fdz_len-1];

  // loop x/y face
  for (int iface=0; iface<4; iface++)
  {
    // skip to next face if not cfspml
    if (boundary_itype[iface] != FD_BOUNDARY_TYPE_CFSPML) continue;

    // get index into local var
    size_t abs_ni1 = abs_indx[iface+FD_NDIM_2+0];
    size_t abs_ni2 = abs_indx[iface+FD_NDIM_2+1];
    size_t abs_nj1 = abs_indx[iface+FD_NDIM_2+2];
    size_t abs_nj2 = abs_indx[iface+FD_NDIM_2+3];
    size_t abs_nk1 = abs_indx[iface+FD_NDIM_2+4];
    size_t abs_nk2 = abs_indx[iface+FD_NDIM_2+5];

    // max of two values
    abs_nk1 = (k_min > abs_nk1 ) ? k_min : abs_nk1;

    // get coef for this face
    float *restrict ptr_ceof_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    size_t pos_cur_face = abs_vars_facepos0[iface];
    size_t pml_volsiz   = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;
    float *restrict pml_hVx  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_hVy  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_hVz  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;

    // for each dim
    if (iface < 2) // x direction
    {
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
            int abs_i = i - abs_ni1;
            coef_D = ptr_coef_D[abs_i];
            coef_A = ptr_coef_A[abs_i];
            coef_B = ptr_coef_B[abs_i];

            // metric
            xix = xi_x[iptr];
            xiy = xi_y[iptr];
            xiz = xi_z[iptr];

            // slowness and jac
            slwjac = slw3d[iptr] * jac3d[iptr];

            // for hVx
            // transform to conservative vars
            for (n=0; n<fdx_len; n++) {
              iptr4vec = iptr + fdx_indx[n];
              vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txx[iptr4vec]
                                            + xi_y[iptr4vec] * Txy[iptr4vec]
                                            + xi_z[iptr4vec] * Txz[iptr4vec] );
            }
            // fd
            M_FD_NOINDX(DxCx, vecxi, fdx_len, lfdx_coef, n_fd);
            // combine for corr and aux vars
            hVx_rhs = DxCx * slwjac;
            // make corr to moment equation
            hVx[iptr] += (coef_B - 1.0) * hVx_rhs - coef_B * pml_Vx[iptr_a];
            // aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVx[iptr_a] = d1 * hVx_rhs - a1 * pml_Vx[iptr_a];

            // for hVy
            // transform to conservative vars
            for (n=0; n<fdx_len; n++) {
              iptr4vec = iptr + fdx_indx[n];
              vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txy[iptr4vec]
                                            + xi_y[iptr4vec] * Tyy[iptr4vec]
                                            + xi_z[iptr4vec] * Tyz[iptr4vec] );
            }
            // fd
            M_FD_NOINDX(DxCx, vecxi, fdx_len, lfdx_coef, n_fd);
            // combine for corr and aux vars
            hVy_rhs = DxCx * slwjac;
            // make corr to moment equation
            hVy[iptr] += (coef_B - 1.0) * hVy_rhs - coef_B * pml_Vy[iptr_a];
            // aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVy[iptr_a] = d1 * hVy_rhs - a1 * pml_Vy[iptr_a];

            //
            // for hVz
            // transform to conservative vars
            for (n=0; n<fdx_len; n++) {
              iptr4vec = iptr + fdx_indx[n];
              vecxi[n] = jac3d[iptr4vec] * (  xi_x[iptr4vec] * Txz[iptr4vec]
                                            + xi_y[iptr4vec] * Tyz[iptr4vec]
                                            + xi_z[iptr4vec] * Tzz[iptr4vec] );
            }
            // fd
            M_FD_NOINDX(DxCx, vecxi, fdx_len, lfdx_coef, n_fd);
            // combine for corr and aux vars
            hVz_rhs = DxCx * slwjac;
            // make corr to moment equation
            hVz[iptr] += (coef_B - 1.0) * hVz_rhs - coef_B * pml_Vz[iptr_a];
            // aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVz[iptr_a] = d1 * hVz_rhs - a1 * pml_Vz[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
    else // y direction
    {
      iptr_a = 0;
      for (k=abs_nk1; k<=abs_nk2; k++)
      {
        iptr_k = k * siz_slice;
        for (j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;

          // pml coefs
          int abs_j = j - abs_nj1;
          coef_D = ptr_coef_D[abs_j];
          coef_A = ptr_coef_A[abs_j];
          coef_B = ptr_coef_B[abs_j];

          for (i=abs_ni1; i<=abs_ni2; i++)
          {
            // metric
            etx = et_x[iptr];
            ety = et_y[iptr];
            etz = et_z[iptr];

            // slowness and jac
            slwjac = slw3d[iptr] * jac3d[iptr];

            // for hVx
            // transform to conservative vars
            for (n=0; n<fdy_len; n++) {
              iptr4vec = iptr + fdy_indx[n] * siz_line;
              vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txx[iptr4vec]
                                            + et_y[iptr4vec] * Txy[iptr4vec]
                                            + et_z[iptr4vec] * Txz[iptr4vec] );
            }
            // fd
            M_FD_NOINDX(DyCy, vecet, fdy_len, lfdy_coef, n_fd);
            // combine for corr and aux vars
            hVx_rhs = DyCy * slwjac;
            // make corr to moment equation
            hVx[iptr] += (coef_B - 1.0) * hVx_rhs - coef_B * pml_Vx[iptr_a];
            // aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVx[iptr_a] = d1 * hVx_rhs - a1 * pml_Vx[iptr_a];

            // for hVy
            // transform to conservative vars
            for (n=0; n<fdy_len; n++) {
              iptr4vec = iptr + fdy_indx[n] * siz_line;
              vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txy[iptr4vec]
                                            + et_y[iptr4vec] * Tyy[iptr4vec]
                                            + et_z[iptr4vec] * Tyz[iptr4vec] );
            }
            // fd
            M_FD_NOINDX(DyCy, vecet, fdy_len, lfdy_coef, n_fd);
            // combine for corr and aux vars
            hVy_rhs = DyCy * slwjac;
            // make corr to moment equation
            hVy[iptr] += (coef_B - 1.0) * hVy_rhs - coef_B * pml_Vy[iptr_a];
            // aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVy[iptr_a] = d1 * hVy_rhs - a1 * pml_Vy[iptr_a];

            // for hVz
            // transform to conservative vars
            for (n=0; n<fdy_len; n++) {
              iptr4vec = iptr + fdy_indx[n] * siz_line;
              vecet[n] = jac3d[iptr4vec] * (  et_x[iptr4vec] * Txz[iptr4vec]
                                            + et_y[iptr4vec] * Tyz[iptr4vec]
                                            + et_z[iptr4vec] * Tzz[iptr4vec] );
            }
            // fd
            M_FD_NOINDX(DyCy, vecet, fdy_len, lfdy_coef, n_fd);
            // combine for corr and aux vars
            hVz_rhs = DyCy * slwjac;
            // make corr to moment equation
            hVz[iptr] += (coef_B - 1.0) * hVz_rhs - coef_B * pml_Vz[iptr_a];
            // aux var
            //   a1 = alpha + d / beta, dealt in abs_set_cfspml
            pml_hVz[iptr_a] = d1 * hVz_rhs - a1 * pml_Vz[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    } // if x or y dim
  } // face
}

// considering velocity free surface condition in cfs-pml
void
sv_eliso1st_curv_macdrp_rhs_cfspml_vfree_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    size_t nk2, size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose)
{
  // loop var for fd
  size_t n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];

  // val on point
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // local
  size_t i,j,k;
  size_t iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D;
  float DxVx,DxVy,DxVz;
  float DyVx,DyVy,DyVz;
  float Dx_DzVx,Dy_DzVx,Dx_DzVy,Dy_DzVy,Dx_DzVz,Dy_DzVz;
  float hTxx_rhs,hTyy_rhs,hTzz_rhs,hTxz_rhs,hTyz_rhs,hTxy_rhs

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }

  // loop x/y face
  for (int iface=0; iface<4; iface++)
  {
    // skip to next face if not cfspml
    if (boundary_itype[iface] != FD_BOUNDARY_TYPE_CFSPML) continue;

    // get index into local var
    size_t abs_ni1 = abs_indx[iface+FD_NDIM_2+0];
    size_t abs_ni2 = abs_indx[iface+FD_NDIM_2+1];
    size_t abs_nj1 = abs_indx[iface+FD_NDIM_2+2];
    size_t abs_nj2 = abs_indx[iface+FD_NDIM_2+3];
    size_t abs_nk1 = abs_indx[iface+FD_NDIM_2+4];
    size_t abs_nk2 = abs_indx[iface+FD_NDIM_2+5];

    // get coef for this face
    float *restrict ptr_ceof_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    size_t pos_cur_face = abs_vars_facepos0[iface];
    size_t pml_volsiz = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;
    float *restrict pml_hTxx = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TXX * pml_volsiz;
    float *restrict pml_hTyy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TYY * pml_volsiz;
    float *restrict pml_hTzz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TZZ * pml_volsiz;
    float *restrict pml_hTxz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TXZ * pml_volsiz;
    float *restrict pml_hTyz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TYZ * pml_volsiz;
    float *restrict pml_hTxy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_TXY * pml_volsiz;

    // for each dim
    if (iface < 2) // x direction
    {
      iptr_a = 0;

      if (nk2 > abs_nk2) continue;

      // set k to free surface
      k = nk2;
      iptr_k = k * siz_slice;

      for (j=abs_nj1; j<=abs_nj2; j++)
      {
        iptr_j = iptr_k + j * siz_line;
        iptr = iptr_j + abs_ni1;
        for (i=abs_ni1; i<=abs_ni2; i++)
        {
          // pml coefs
          int abs_i = i - abs_ni1;
          coef_D = ptr_coef_D[abs_i];
          coef_A = ptr_coef_A[abs_i];
          coef_B = ptr_coef_B[abs_i];

          // metric
          xix = xi_x[iptr];
          xiy = xi_y[iptr];
          xiz = xi_z[iptr];

          // medium
          lam = lam3d[iptr];
          mu  =  mu3d[iptr];
          slw = slw3d[iptr];
          lam2mu = lam + 2.0 * mu;

          // Vx derivatives
          M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
          // Vy derivatives
          M_FD_SHIFT(DxVy, Vy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
          // Vz derivatives
          M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

          // zeta derivatives
          size_t ij = i + j * siz_line;
          Dx_DzVx = matVx2Vz[ij+3*0+0] * DxVx
                  + matVx2Vz[ij+3*0+1] * DxVy
                  + matVx2Vz[ij+3*0+2] * DxVz;

          Dx_DzVy = matVx2Vz[ij+3*1+0] * DxVx
                  + matVx2Vz[ij+3*1+1] * DxVy
                  + matVx2Vz[ij+3*1+2] * DxVz;

          Dx_DzVz = matVx2Vz[ij+3*2+0] * DxVx
                  + matVx2Vz[ij+3*2+1] * DxVy
                  + matVx2Vz[ij+3*2+2] * DxVz;

          // keep xi derivative terms, including free surface convered
          hTxx_rhs =    lam2mu * ( xix*DxVx + ztx*Dx_DzVx);
                      + lam    * ( xiy*DxVy + zty*Dx_DzVy
                                  +xiz*DxVz + ztz*Dx_DzVz);

          hTyy_rhs =   lam2mu * ( xiy*DxVy + zty*Dx_DzVy)
                      +lam    * ( xix*DxVx + ztx*Dx_DzVx
                                 +xiz*DxVz + ztz*Dx_DzVz);

          hTzz_rhs =   lam2mu * ( xiz*DxVz + ztz*Dx_DzVz)
                      +lam    * ( xix*DxVx + ztx*Dx_DzVx
                                 +xiy*DxVy + zty*Dx_DzVy);

          hTxy_rhs = mu *(
                       xiy*DxVx + xix*DxVy
                      +zty*Dx_DzVx + ztx*Dx_DzVy
                      );
          hTxz_rhs = mu *(
                       xiz*DxVx + xix*DxVz
                      +ztz*Dx_DzVx + ztx*Dx_DzVz
                      );
          hTyz_rhs = mu *(
                       xiz*DxVy + xiy*DxVz
                      +ztz*Dx_DzVy + zty*Dx_DzVz
                      );

          // make corr to Hooke's equatoin
          hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs - coef_B * pml_Txx[iptr_a];
          hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
          hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
          hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs - coef_B * pml_Txz[iptr_a];
          hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
          hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs - coef_B * pml_Txy[iptr_a];

          // aux var
          //   a1 = alpha + d / beta, dealt in abs_set_cfspml
          pml_hTxx[iptr_a] = d1 * hTxx_rhs - a1 * pml_Txx[iptr_a];
          pml_hTyy[iptr_a] = d1 * hTyy_rhs - a1 * pml_Tyy[iptr_a];
          pml_hTzz[iptr_a] = d1 * hTzz_rhs - a1 * pml_Tzz[iptr_a];
          pml_hTxz[iptr_a] = d1 * hTxz_rhs - a1 * pml_Txz[iptr_a];
          pml_hTyz[iptr_a] = d1 * hTyz_rhs - a1 * pml_Tyz[iptr_a];
          pml_hTxy[iptr_a] = d1 * hTxy_rhs - a1 * pml_Txy[iptr_a];

          // incr index
          iptr   += 1;
          iptr_a += 1;
        } // i
      } // j
    }
    else // y direction
    {
      iptr_a = 0;

      if (nk2 > abs_nk2) continue;

      // set k to free surface
      k = nk2;
      iptr_k = k * siz_slice;

      for (j=abs_nj1; j<=abs_nj2; j++)
      {
        iptr_j = iptr_k + j * siz_line;
        iptr = iptr_j + abs_ni1;

        // pml coefs
        int abs_j = j - abs_nj1;
        coef_D = ptr_coef_D[abs_j];
        coef_A = ptr_coef_A[abs_j];
        coef_B = ptr_coef_B[abs_j];

        for (i=abs_ni1; i<=abs_ni2; i++)
        {
          // metric
          etx = et_x[iptr];
          ety = et_y[iptr];
          etz = et_z[iptr];

          // medium
          lam = lam3d[iptr];
          mu  =  mu3d[iptr];
          slw = slw3d[iptr];
          lam2mu = lam + 2.0 * mu;

          // Vx derivatives
          M_FD_SHIFT(DyVx, Vx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
          // Vy derivatives
          M_FD_SHIFT(DyVy, Vy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
          // Vz derivatives
          M_FD_SHIFT(DyVz, Vz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

          // zeta derivatives
          size_t ij = i + j * siz_line;
          Dy_DzVx = matVy2Vz[ij+3*0+0] * DyVx
                  + matVy2Vz[ij+3*0+1] * DyVy
                  + matVy2Vz[ij+3*0+2] * DyVz;

          Dy_DzVy = matVy2Vz[ij+3*1+0] * DyVx
                  + matVy2Vz[ij+3*1+1] * DyVy
                  + matVy2Vz[ij+3*1+2] * DyVz;

          Dy_DzVz = matVy2Vz[ij+3*2+0] * DyVx
                  + matVy2Vz[ij+3*2+1] * DyVy
                  + matVy2Vz[ij+3*2+2] * DyVz;

          hTxx_rhs =    lam2mu * (  etx*DyVx + ztx*Dy_DzVx);
                      + lam    * (  ety*DyVy + zty*Dy_DzVy
                                  + etz*DyVz + ztz*Dy_DzVz);

          hTyy_rhs =   lam2mu * (  ety*DyVy + zty*Dy_DzVy)
                      +lam    * (  etx*DyVx + ztx*Dy_DzVx
                                 + etz*DyVz + ztz*Dy_DzVz);

          hTzz_rhs =   lam2mu * (  etz*DyVz + ztz*Dy_DzVz)
                      +lam    * (+ etx*DyVx + ztx*Dy_DzVx
                                 + ety*DyVy + zty*Dy_DzVy);

          hTxy_rhs = mu *(
                       ety*DyVx + etx*DyVy
                      +zty*Dy_DzVx + ztx*Dy_DzVy
                      );
          hTxz_rhs = mu *(
                       etz*DyVx + etx*DyVz
                      +ztz*Dy_DzVx + ztx*Dy_DzVz
                      );
          hTyz_rhs = mu *(
                       etz*DyVy + ety*DyVz
                      +ztz*Dy_DzVy + zty*Dy_DzVz
                    );

          // make corr to Hooke's equatoin
          hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs - coef_B * pml_Txx[iptr_a];
          hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
          hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
          hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs - coef_B * pml_Txz[iptr_a];
          hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
          hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs - coef_B * pml_Txy[iptr_a];

          // aux var
          //   a1 = alpha + d / beta, dealt in abs_set_cfspml
          pml_hTxx[iptr_a] = d1 * hTxx_rhs - a1 * pml_Txx[iptr_a];
          pml_hTyy[iptr_a] = d1 * hTyy_rhs - a1 * pml_Tyy[iptr_a];
          pml_hTzz[iptr_a] = d1 * hTzz_rhs - a1 * pml_Tzz[iptr_a];
          pml_hTxz[iptr_a] = d1 * hTxz_rhs - a1 * pml_Txz[iptr_a];
          pml_hTyz[iptr_a] = d1 * hTyz_rhs - a1 * pml_Tyz[iptr_a];
          pml_hTxy[iptr_a] = d1 * hTxy_rhs - a1 * pml_Txy[iptr_a];

          // incr index
          iptr   += 1;
          iptr_a += 1;
        } // i
      } // j
    } // if x or y dim
  } // face
}

// considering vlow free surface in cfs-pml, no need to deal with points below free
/*
int
sv_eliso1st_curv_macdrp_rhs_cfspml_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    size_t nk2, size_t siz_line, size_t siz_slice,
    size_t fdx_len, size_t *restrict fdx_indx, float *restrict fdx_coef,
    size_t fdy_len, size_t *restrict fdy_indx, float *restrict fdy_coef,
    size_t fdz_num_surf_lay, size_t fdz_max_len, size_t **restrict fdz_all_info,
    size_t *restrict fdz_all_indx, float *restrict fdz_all_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    size_t *restrict abs_indx, // pml range of each face
    size_t *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    size_t *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    size_t *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs)
{
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  size_t lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  size_t lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  size_t lfdz_shift[fdz_len];

  // loop var for fd
  size_t i,j,k;
  size_t n_fd; // loop var for fd
  size_t fdz_len;

  // local var
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // put fd op into local array
  for (int i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i];
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j];
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k];
    lfdz_shift[k] = fdz_indx[k] * siz_slice;
  }

  // loop x/y face
  for (int iface=0; iface<4; iface++)
  {
    // skip to next face if not cfspml
    if (boundary_itype[iface] != FD_BOUNDARY_TYPE_CFSPML) continue;

    // get index into local var
    size_t abs_ni1 = abs_indx[iface][0];
    size_t abs_ni2 = abs_indx[iface][1];
    size_t abs_nj1 = abs_indx[iface][2];
    size_t abs_nj2 = abs_indx[iface][3];
    size_t abs_nk1 = abs_indx[iface][4];
    size_t abs_nk2 = abs_indx[iface][5];

    // get coef for this face
    float *restrict ptr_ceof_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    size_t pos_cur_face = abs_vars_facepos0[iface];
    size_t pml_volsiz = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;
    float *restrict pml_hVx  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VX  * pml_volsiz;
    float *restrict pml_hVy  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VY  * pml_volsiz;
    float *restrict pml_hVz  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_VZ  * pml_volsiz;

    // for each dim
    if (iface < 2) // x direction
    {
      iptr_a = 0;
      for (size_t n=0; n < fdz_num_surf_lay; n++)
      {
        // conver to k index, from surface to inner
        k = nk2 - n;
        iptr_k = k * siz_slice;

        // get different hlen fz near surface
        size_t  lfdz_pos0 = fdz_info[n][FD_INFO_POS_OF_INDX];
        size_t  lfdz_len  = fdz_info[n][FD_INFO_LENGTH_TOTAL];
        size_t  *p_fdz_indx  = fdz_all_indx[lfdz_pos0];
        float   *p_fdz_coef  = fdz_all_coef[lfdz_pos0];
        for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
          lfdz_shift[n_fd] = p_fdz_indx[n_fd] * siz_slice;
          lfdz_coef[n_fd]  = p_fdz_coef[n_fd];
        }

        for (j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;
          for (i=abs_ni1; i<=abs_ni2; i++)
          {
            // pml coefs
            int abs_i = i - abs_ni1;
            coef_D = ptr_coef_D[abs_i];
            coef_A = ptr_coef_A[abs_i];
            coef_B = ptr_coef_B[abs_i];

            // metric
            xix = xi_x[iptr];
            xiy = xi_y[iptr];
            xiz = xi_z[iptr];

            // medium
            lam = lam3d[iptr];
            mu  =  mu3d[iptr];
            slw = slw3d[iptr];
            lam2mu = lam + 2.0 * mu;

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
    else // y direction
    {
      iptr_a = 0;
      for (k=abs_nk1; k<=abs_nk2; k++)
      {
        iptr_k = k * siz_slice;
        for (j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;

          // pml coefs
          int abs_j = j - abs_nj1;
          coef_D = ptr_coef_D[abs_j];
          coef_A = ptr_coef_A[abs_j];
          coef_B = ptr_coef_B[abs_j];

          for (i=abs_ni1; i<=abs_ni2; i++)
          {
            // metric
            etx = et_x[iptr];
            ety = et_y[iptr];
            etz = et_z[iptr];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    } // if x or y dim
  } // face
}
*/

/*******************************************************************************
 * ABL-EXP boundary
 ******************************************************************************/

/* 
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
*/

/*******************************************************************************
 * add source terms
 ******************************************************************************/

int
sv_eliso1st_curv_macdrp_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    size_t siz_line, size_t siz_slice,
    int num_of_force,
    int *restrict force_loc_point,
    float *restrict force_vec_value,
    int num_of_moment,
    int *restrict moment_loc_point,
    float *restrict moment_ten_value, // size: num_of_moment * 6
    const int myid, const int verbose)
{
  int ierr = 0;

  // local var
  size_t si,sj,sk, iptr;

  // add force
  for (int n=0; n<num_of_force; n++)
  {
    size_t *ptr_force_loc = force_loc_point + n * FD_NDIM;
    size_t *ptr_force_val = force_vec_value + n * 3;

    si = force_loc_point[0];
    sj = force_loc_point[1];
    sk = force_loc_point[2];

    iptr = si + sj * siz_line + sk * siz_slice;

    float V = slw3d[iptr] / jac3d[iptr];

    hVx[iptr] += ptr_force_val[ 0 ] * V;
    hVy[iptr] += ptr_force_val[ 1 ] * V;
    hVz[iptr] += ptr_force_val[ 2 ] * V;
  }

  // add moment source
  for (int n=0; n<num_of_moment; n++)
  {
    size_t *ptr_moment_loc = moment_loc_point + n * FD_NDIM;
    size_t *ptr_moment_val = moment_ten_value + n * 6;

    si = moment_loc_point[0];
    sj = moment_loc_point[1];
    sk = moment_loc_point[2];

    iptr = si + sj * siz_line + sk * siz_slice;

    float jac = jac3d[iptr];

    hTxx[iptr] -= ptr_moment_val[ 0 ] / jac;
    hTyy[iptr] -= ptr_moment_val[ 1 ] / jac;
    hTzz[iptr] -= ptr_moment_val[ 2 ] / jac;
    hTxz[iptr] -= ptr_moment_val[ 3 ] / jac;
    hTyz[iptr] -= ptr_moment_val[ 4 ] / jac;
    hTxy[iptr] -= ptr_moment_val[ 5 ] / jac;
  }

  return ierr;
}

/*******************************************************************************
 * related functions
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_vel_dxy2dz(
    float *restrict g3d,
    float *restrict m3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    float **restrict p_matVx2Vz,
    float **restrict p_matVy2Vz,
    const int myid, const int verbose)
{
  float e11, e12, e13, e21, e22, e23, e31, e32, e33;
  float lam2mu, lam, mu;
  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];

  (float *) matVx2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * FD_NDIM * FD_NDIM,
                                      0.0,
                                      "sv_eliso1st_curv_macdrp_vel_dxy2dz");

  (float *) matVy2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * FD_NDIM * FD_NDIM,
                                      0.0,
                                      "sv_eliso1st_curv_macdrp_vel_dxy2dz");

  float *restrict xi_x = g3d + GD_CURV_SEQ_XIX * siz_volume;
  float *restrict xi_y = g3d + GD_CURV_SEQ_XIY * siz_volume;
  float *restrict xi_z = g3d + GD_CURV_SEQ_XIZ * siz_volume;
  float *restrict et_x = g3d + GD_CURV_SEQ_ETX * siz_volume;
  float *restrict et_y = g3d + GD_CURV_SEQ_ETY * siz_volume;
  float *restrict et_z = g3d + GD_CURV_SEQ_ETZ * siz_volume;
  float *restrict zt_x = g3d + GD_CURV_SEQ_ZTX * siz_volume;
  float *restrict zt_y = g3d + GD_CURV_SEQ_ZTY * siz_volume;
  float *restrict zt_z = g3d + GD_CURV_SEQ_ZTZ * siz_volume;
  float *restrict lam3d = m3d + MD_ELISO_SEQ_LAMBDA * siz_volume;
  float *restrict  mu3d = m3d + MD_ELISO_SEQ_MU     * siz_volume;

  size_t k = nk2;

  for (size_t j = nj1; j < nj2; j++)
  {
    for (size_t i = ni1; i < ni2; i++)
    {
      size_t iptr = i + j * siz_line + k * siz_slice;

      e11 = xi_x[iptr];
      e12 = xi_y[iptr];
      e13 = xi_z[iptr];
      e21 = et_x[iptr];
      e22 = et_y[iptr];
      e23 = et_z[iptr];
      e31 = zt_x[iptr];
      e32 = zt_y[iptr];
      e33 = zt_z[iptr];

      lam    = lam3d[iptr];
      mu     =  mu3d[iptr];
      lam2mu = lam + 2.0f * mu;

      // first dim: x; sec dim: y, as Fortran code
      A[0][0] = lam2mu*e31*e31 + mu*(e32*e32+e33*e33);
      A[0][1] = lam*e31*e32 + mu*e32*e31;
      A[0][2] = lam*e31*e33 + mu*e33*e31;
      A[1][0] = lam*e32*e31 + mu*e31*e32;
      A[1][1] = lam2mu*e32*e32 + mu*(e31*e31+e33*e33);
      A[1][2] = lam*e32*e33 + mu*e33*e32;
      A[2][0] = lam*e33*e31 + mu*e31*e33;
      A[2][1] = lam*e33*e32 + mu*e32*e33;
      A[2][2] = lam2mu*e33*e33 + mu*(e31*e31+e32*e32);
      fdlib_math_invert3x3(A);

      B[0][0] = -lam2mu*e31*e11 - mu*(e32*e12+e33*e13);
      B[0][1] = -lam*e31*e12 - mu*e32*e11;
      B[0][2] = -lam*e31*e13 - mu*e33*e11;
      B[1][0] = -lam*e32*e11 - mu*e31*e12;
      B[1][1] = -lam2mu*e32*e12 - mu*(e31*e11+e33*e13);
      B[1][2] = -lam*e32*e13 - mu*e33*e12;
      B[2][0] = -lam*e33*e11 - mu*e31*e13;
      B[2][1] = -lam*e33*e12 - mu*e32*e13;
      B[2][2] = -lam2mu*e33*e13 - mu*(e31*e11+e32*e12);

      C[0][0] = -lam2mu*e31*e21 - mu*(e32*e22+e33*e23);
      C[0][1] = -lam*e31*e22 - mu*e32*e21;
      C[0][2] = -lam*e31*e23 - mu*e33*e21;
      C[1][0] = -lam*e32*e21 - mu*e31*e22;
      C[1][1] = -lam2mu*e32*e22 - mu*(e31*e21+e33*e23);
      C[1][2] = -lam*e32*e23 - mu*e33*e22;
      C[2][0] = -lam*e33*e21 - mu*e31*e23;
      C[2][1] = -lam*e33*e22 - mu*e32*e23;
      C[2][2] = -lam2mu*e33*e23 - mu*(e31*e21+e32*e22);

      fdlib_math_matmul3x3(A, B, AB);
      fdlib_math_matmul3x3(A, C, AC);

      size_t ij = (i * siz_line + j) * 9;

      // save into mat
      for(int m = 0; m < 3; m++)
        for(int n = 0; n < 3; n++){
          matVx2Vz[ij + m*3 + n] = AB[m][n];
          matVy2Vz[ij + m*3 + n] = AC[m][n];
        }
    }
  }

  *p_matVx2Vz = matVx2Vz;
  *p_matVy2Vz = matVy2Vz;
}
