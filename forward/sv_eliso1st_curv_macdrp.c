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
#include "io_funcs.h"
#include "sv_eliso1st_curv_macdrp.h"

//#define M_NCERR(ierr) {fprintf(stderr,"sv_ nc error: %s\n", nc_strerror(ierr)); exit(1);}
#ifndef M_NCERR
#define M_NCERR {fprintf(stderr,"sv_ nc error\n"); exit(1);}
#endif

//#define SV_ELISO1ST_CURV_MACDRP_DEBUG

//#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
//#include "fd_debug.h"
//#endif

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_allstep(
    float *restrict w3d,  // wavefield
    size_t *restrict w3d_pos,
    char **w3d_name,
    int   w3d_num_of_vars,
    char **coord_name,
    float *restrict g3d,  // grid vars
    float *restrict m3d,  // medium vars
    // grid size
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    int ni, int nj, int nk, int nx, int ny, int nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    // boundary type
    int *restrict boundary_itype,
    // if abs
    int              abs_itype, //
    int    *restrict abs_num_of_layers, //
    int *restrict abs_indx, //
    int *restrict abs_coefs_facepos0, //
    float  *restrict abs_coefs, //
    int           abs_vars_size_per_level, //
    int *restrict abs_vars_volsiz, //
    int *restrict abs_vars_facepos0, //
    float  *restrict abs_vars,
    // if free surface
    float *matVx2Vz, //
    float *matVy2Vz, //
    // source term
    int num_of_force,
    int *restrict force_info, // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
    float *restrict force_vec_stf,
    int   *restrict force_ext_indx,
    float *restrict force_ext_coef,
    int             num_of_moment,
    int   *restrict moment_info, // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
    float *restrict moment_ten_rate,
    int   *restrict moment_ext_indx,
    float *restrict moment_ext_coef,
    // io
    int num_of_sta, int *restrict sta_loc_indx, float *restrict sta_loc_dxyz, float *restrict sta_seismo,
    int num_of_point, int *restrict point_loc_indx, float *restrict point_seismo,
    int num_of_slice_x, int *restrict slice_x_indx, char **restrict slice_x_fname,
    int num_of_slice_y, int *restrict slice_y_indx, char **restrict slice_y_fname,
    int num_of_slice_z, int *restrict slice_z_indx, char **restrict slice_z_fname,
    int num_of_snap, int *restrict snap_info, char **snap_fname,
    char *output_fname_part,
    char *output_dir,
    // scheme
    int num_rk_stages, float *rk_a, float *rk_b,
    int num_of_pairs, 
    int fdx_max_half_len, int fdy_max_half_len,
    int fdz_max_len, int fdz_num_surf_lay,
    int ****pair_fdx_all_info, int ***pair_fdx_all_indx, float ***pair_fdx_all_coef,
    int ****pair_fdy_all_info, int ***pair_fdy_all_indx, float ***pair_fdy_all_coef,
    int ****pair_fdz_all_info, int ***pair_fdz_all_indx, float ***pair_fdz_all_coef,
    // time
    float dt, int nt_total, float t0,
    // mpi
    int myid, int *myid2, MPI_Comm comm,
    float *restrict sbuff,
    float *restrict rbuff,
    MPI_Request *s_reqs,
    MPI_Request *r_reqs,
    int qc_check_nan_num_of_step,
    const int output_all, // qc all var
    const int verbose)
{
  // local allocated array
  float *force_vec_value  = NULL;  // num_of_force * 3
  float *moment_ten_value = NULL;  // num_of_moment * 6
  char ou_file[FD_MAX_STRLEN];

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
  float t_next; // time after this loop for nc output

  size_t w3d_size_per_level = WF_EL_1ST_NVAR * siz_volume;

  // io nc
  int ncid_slx[num_of_slice_x];
  int timeid_slx;
  int varid_slx[w3d_num_of_vars];

  int ncid_sly[num_of_slice_y];
  int timeid_sly;
  int varid_sly[w3d_num_of_vars];

  int ncid_slz[num_of_slice_z];
  int timeid_slz;
  int varid_slz[w3d_num_of_vars];

  int ncid_snap[num_of_snap];
  int timeid_snap[num_of_snap];
  int varid_snap_vel[num_of_snap*FD_NDIM];
  int varid_snap_T[num_of_snap*FD_NDIM_2];
  int varid_snap_E[num_of_snap*FD_NDIM_2];
  int snap_cur_it[num_of_snap];
  for (int n=0; n<num_of_snap; n++) {
    snap_cur_it[n] = 0;
  }

  // only x/y mpi
  int num_of_r_reqs = 4;
  int num_of_s_reqs = 4;

  // alloc
  if (num_of_force > 0) {
    force_vec_value = (float *)fdlib_mem_malloc_1d(num_of_force*3*sizeof(float),
                                                   "alloc force_vec_value in all step");
  }
  if (num_of_moment > 0) {
    moment_ten_value = (float *)fdlib_mem_malloc_1d(num_of_moment*6*sizeof(float),
                                                    "alloc moment_ten_value in all step");
  }

  // get wavefield
  w_pre = w3d + w3d_size_per_level * 0; // previous level at n
  w_tmp = w3d + w3d_size_per_level * 1; // intermidate value
  w_rhs = w3d + w3d_size_per_level * 2; // for rhs
  w_end = w3d + w3d_size_per_level * 3; // end level at n+1

  // get pml
  abs_vars_pre = abs_vars + abs_vars_size_per_level * 0;
  abs_vars_tmp = abs_vars + abs_vars_size_per_level * 1;
  abs_vars_rhs = abs_vars + abs_vars_size_per_level * 2;
  abs_vars_end = abs_vars + abs_vars_size_per_level * 3;

  // create slice and snapshot nc
  for (int n=0; n<num_of_slice_x; n++)
  {
    int dimid[3];
    if (nc_create(slice_x_fname[n], NC_CLOBBER, &ncid_slx[n])) M_NCERR;
    if (nc_def_dim(ncid_slx[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid_slx[n], "k", nk      , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid_slx[n], "j" , nj      , &dimid[2])) M_NCERR;
    // time var
    if (nc_def_var(ncid_slx[n], "time", NC_FLOAT, 1, dimid+0, &timeid_slx)) M_NCERR;
    // other vars
    for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
      if (nc_def_var(ncid_slx[n], w3d_name[ivar], NC_FLOAT, 3, dimid, &varid_slx[ivar])) M_NCERR;
    }
    // attribute: index info for plot
    nc_put_att_int(ncid_slx[n],NC_GLOBAL,"i_index_with_ghosts_in_this_thread",
                   NC_INT,1,slice_x_indx+n);
    nc_put_att_int(ncid_slx[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,myid2);
    // end def
    if (nc_enddef(ncid_slx[n])) M_NCERR;
  }
  for (int n=0; n<num_of_slice_y; n++) {
    //int sli_id = slice_y_info[n*2+0];
    int dimid[3];
    if (nc_create(slice_y_fname[n], NC_CLOBBER, &ncid_sly[n])) M_NCERR;
    if (nc_def_dim(ncid_sly[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid_sly[n], "k", nk      , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid_sly[n], "i" , ni      , &dimid[2])) M_NCERR;
    // time var
    if (nc_def_var(ncid_sly[n], "time", NC_FLOAT, 1, dimid+0, &timeid_sly)) M_NCERR;
    // other vars
    for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
      if (nc_def_var(ncid_sly[n], w3d_name[ivar], NC_FLOAT, 3, dimid, &varid_sly[ivar])) M_NCERR;
    }
    // attribute: index info for plot
    nc_put_att_int(ncid_sly[n],NC_GLOBAL,"j_index_with_ghosts_in_this_thread",
                   NC_INT,1,slice_y_indx+n);
    nc_put_att_int(ncid_sly[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,myid2);
    // end def
    if (nc_enddef(ncid_sly[n])) M_NCERR;
  }
  for (int n=0; n<num_of_slice_z; n++) {
    //int sli_id = slice_z_info[n*2+0];
    int dimid[3];
    if (nc_create(slice_z_fname[n], NC_CLOBBER, &ncid_slz[n])) M_NCERR;
    if (nc_def_dim(ncid_slz[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid_slz[n], "j", nj      , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid_slz[n], "i" , ni      , &dimid[2])) M_NCERR;
    // time var
    if (nc_def_var(ncid_slz[n], "time", NC_FLOAT, 1, dimid+0, &timeid_slz)) M_NCERR;
    // other vars
    for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
      if (nc_def_var(ncid_slz[n], w3d_name[ivar], NC_FLOAT, 3, dimid, &varid_slz[ivar])) M_NCERR;
    }
    // attribute: index info for plot
    nc_put_att_int(ncid_slz[n],NC_GLOBAL,"k_index_with_ghosts_in_this_thread",
                   NC_INT,1,slice_z_indx+n);
    nc_put_att_int(ncid_slz[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,myid2);
    // end def
    if (nc_enddef(ncid_slz[n])) M_NCERR;
  }
  // snapshot
  for (int n=0; n<num_of_snap; n++)
  {
    int dimid[4];
    int *cur_snap_info = snap_info + n * FD_SNAP_INFO_SIZE;
    int snap_i1  = cur_snap_info[FD_SNAP_INFO_I1];
    int snap_j1  = cur_snap_info[FD_SNAP_INFO_J1];
    int snap_k1  = cur_snap_info[FD_SNAP_INFO_K1];
    int snap_ni  = cur_snap_info[FD_SNAP_INFO_NI];
    int snap_nj  = cur_snap_info[FD_SNAP_INFO_NJ];
    int snap_nk  = cur_snap_info[FD_SNAP_INFO_NK];
    int snap_di  = cur_snap_info[FD_SNAP_INFO_DI];
    int snap_dj  = cur_snap_info[FD_SNAP_INFO_DJ];
    int snap_dk  = cur_snap_info[FD_SNAP_INFO_DK];

    int snap_it1 = cur_snap_info[FD_SNAP_INFO_IT1];
    int snap_dit = cur_snap_info[FD_SNAP_INFO_DIT];
    int snap_nt_total = (nt_total - snap_it1) / snap_dit;
    int snap_out_V = cur_snap_info[FD_SNAP_INFO_VEL];
    int snap_out_T = cur_snap_info[FD_SNAP_INFO_STRESS];
    int snap_out_E = cur_snap_info[FD_SNAP_INFO_STRAIN];

    if (nc_create(snap_fname[n], NC_CLOBBER, &ncid_snap[n])) M_NCERR;
    if (nc_def_dim(ncid_snap[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid_snap[n], "k", snap_nk     , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid_snap[n], "j" , snap_nj     , &dimid[2])) M_NCERR;
    if (nc_def_dim(ncid_snap[n], "i" , snap_ni      , &dimid[3])) M_NCERR;
    // time var
    if (nc_def_var(ncid_snap[n], "time", NC_FLOAT, 1, dimid+0, &timeid_snap[n])) M_NCERR;
    // other vars
    if (snap_out_V==1) {
       if (nc_def_var(ncid_snap[n],"Vx",NC_FLOAT,4,dimid,&varid_snap_vel[n*FD_NDIM+0])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Vy",NC_FLOAT,4,dimid,&varid_snap_vel[n*FD_NDIM+1])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Vz",NC_FLOAT,4,dimid,&varid_snap_vel[n*FD_NDIM+2])) M_NCERR;
    }
    if (snap_out_T==1) {
       if (nc_def_var(ncid_snap[n],"Txx",NC_FLOAT,4,dimid,&varid_snap_T[n*FD_NDIM_2+0])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Tyy",NC_FLOAT,4,dimid,&varid_snap_T[n*FD_NDIM_2+1])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Tzz",NC_FLOAT,4,dimid,&varid_snap_T[n*FD_NDIM_2+2])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Txz",NC_FLOAT,4,dimid,&varid_snap_T[n*FD_NDIM_2+3])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Tyz",NC_FLOAT,4,dimid,&varid_snap_T[n*FD_NDIM_2+4])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Txy",NC_FLOAT,4,dimid,&varid_snap_T[n*FD_NDIM_2+5])) M_NCERR;
    }
    if (snap_out_E==1) {
       if (nc_def_var(ncid_snap[n],"Exx",NC_FLOAT,4,dimid,&varid_snap_E[n*FD_NDIM_2+0])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Eyy",NC_FLOAT,4,dimid,&varid_snap_E[n*FD_NDIM_2+1])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Ezz",NC_FLOAT,4,dimid,&varid_snap_E[n*FD_NDIM_2+2])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Exz",NC_FLOAT,4,dimid,&varid_snap_E[n*FD_NDIM_2+3])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Eyz",NC_FLOAT,4,dimid,&varid_snap_E[n*FD_NDIM_2+4])) M_NCERR;
       if (nc_def_var(ncid_snap[n],"Exy",NC_FLOAT,4,dimid,&varid_snap_E[n*FD_NDIM_2+5])) M_NCERR;
    }
    // attribute: index in output snapshot, index w ghost in thread
    int g_start[] = { cur_snap_info[FD_SNAP_INFO_GI1],
                      cur_snap_info[FD_SNAP_INFO_GJ1],
                      cur_snap_info[FD_SNAP_INFO_GK1] };
    nc_put_att_int(ncid_snap[n],NC_GLOBAL,"first_index_to_snapshot_output",
                   NC_INT,FD_NDIM,g_start);

    nc_put_att_int(ncid_snap[n],NC_GLOBAL,"first_index_in_this_thread_with_ghosts",
                   NC_INT,FD_NDIM,cur_snap_info+FD_SNAP_INFO_I1);

    nc_put_att_int(ncid_snap[n],NC_GLOBAL,"index_stride_in_this_thread",
                   NC_INT,FD_NDIM,cur_snap_info+FD_SNAP_INFO_DI);
    nc_put_att_int(ncid_snap[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,myid2);

    if (nc_enddef(ncid_snap[n])) M_NCERR;

    //sv_eliso1st_curv_macdrp_snap_create_nc(snap_fname[n],
    //                                       snap_nt_totle,
    //                                       snap_ni,
    //                                       snap_nj,
    //                                       snap_nk,
    //                                       snap_out_vel,
    //                                       snap_out_stress,
    //                                       snap_out_strain,
    //                                       ncid_snap+n,
    //                                       varid_snap_V+n*FD_NDIM,
    //                                       varid_snap_T+n*FD_NDIM_2,
    //                                       varid_snap_E+n*FD_NDIM_2);
  }

  //
  // time loop
  //
  for (int it=0; it<nt_total; it++)
  {
    t_cur = it * dt + t0;
    t_next = t_cur +dt;

    if (myid==0 && verbose>10) fprintf(stdout,"-> it=%d, t=\%f\n", it, t_cur);

    ipair = it % num_of_pairs; // mod to ipair
    if (myid==0 && verbose>10) fprintf(stdout, " --> ipair=%d\n",ipair);

    // loop RK stages for one step
    for (istage=0; istage<num_rk_stages; istage++)
    {
      if (myid==0 && verbose>10) fprintf(stdout, " --> istage=%d\n",istage);

      // use pointer to avoid 1 copy for previous level value
      if (istage!=0) {
        w_cur        = w_tmp;
        abs_vars_cur = abs_vars_tmp;
      } else {
        w_cur        = w_pre;
        abs_vars_cur = abs_vars_pre;
      }

      //// recv mesg
      //if (it>0 || istage>0) {
      //  MPI_Waitall(num_of_r_reqs, r_reqs, MPI_STATUS_IGNORE);

      //  fd_blk_unpack_mesg(rbuff, w_cur, w3d_num_of_vars,
      //                     ni1,ni2,nj1,nj2,nk1,nk2,nx,ny,
      //                     siz_line,siz_slice,siz_volume);

      //  MPI_Startall(num_of_r_reqs, r_reqs);

      //  MPI_Waitall(num_of_s_reqs, s_reqs, MPI_STATUS_IGNORE);
      //}
      // stf value for cur stage
      src_get_stage_stf(num_of_force,
                        force_info,
                        force_vec_stf,
                        num_of_moment,
                        moment_info,
                        moment_ten_rate,
                        it,
                        istage,
                        num_rk_stages,
                        force_vec_value,
                        moment_ten_value,
                        myid, verbose);
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
                                       num_of_force , force_info , force_vec_value,
                                       force_ext_indx,force_ext_coef,
                                       num_of_moment, moment_info, moment_ten_value,
                                       moment_ext_indx,moment_ext_coef,
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

      // cal var for send first
      if (istage < num_rk_stages-1)
      {
        float coef_a = rk_a[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr<WF_EL_1ST_NVAR*siz_volume; iptr++) {
            w_tmp[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
        }
        // apply abs
        /*
        if (abs_has_ablexp) {
          sv_eliso1st_curv_macdrp_ablexp(w_tmp, w3d_num_of_vars,
                    w3d_pos, nx, ny, nz, abs_coefs_facepos0, abs_coefs);
        }
        */

        // pack and isend
        fd_blk_pack_mesg(w_tmp, sbuff, w3d_num_of_vars,
                         ni1,ni2,nj1,nj2,nk1,nk2,
                         siz_line,siz_slice,siz_volume,
                         fdx_max_half_len, fdy_max_half_len);

        MPI_Startall(num_of_s_reqs, s_reqs);
      }
      else // w_end for send at last stage
      {
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr<WF_EL_1ST_NVAR*siz_volume; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }
        /*
        if (abs_has_ablexp) {
          sv_eliso1st_curv_macdrp_ablexp(w_end, w3d_num_of_vars,
                    w3d_pos, nx, ny, nz, abs_coefs_facepos0, abs_coefs);
        }
        */
        // pack and isend
        fd_blk_pack_mesg(w_end, sbuff, w3d_num_of_vars,
                         ni1,ni2,nj1,nj2,nk1,nk2,
                         siz_line,siz_slice,siz_volume,
                         fdx_max_half_len, fdy_max_half_len);

        MPI_Startall(num_of_s_reqs, s_reqs);
      }

      // cal rest
      if (istage==0)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // w_tmp is done
        // pml_tmp
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_tmp[iptr] = abs_vars_pre[iptr] + coef_a * abs_vars_rhs[iptr];
          }
        }

        // w_end
        for (size_t iptr=0; iptr<WF_EL_1ST_NVAR*siz_volume; iptr++) {
            w_end[iptr] = w_pre[iptr] + coef_b * w_rhs[iptr];
        }
        // pml_end
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_end[iptr] = abs_vars_pre[iptr] + coef_b * abs_vars_rhs[iptr];
          }
        }
      }
      else if (istage<num_rk_stages-1)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // w_tmp is done
        // pml_tmp
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_tmp[iptr] = abs_vars_pre[iptr] + coef_a * abs_vars_rhs[iptr];
          }
        }

        // w_end
        for (size_t iptr=0; iptr<WF_EL_1ST_NVAR*siz_volume; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }
        // pml_end
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_end[iptr] += coef_b * abs_vars_rhs[iptr];
          }
        }
      }
      else // last stage
      {
        float coef_b = rk_b[istage] * dt;

        // w_end is done
        // pml_end
        if (abs_itype == FD_BOUNDARY_TYPE_CFSPML)
        {
          for (size_t iptr=0; iptr<abs_vars_size_per_level; iptr++) {
            abs_vars_end[iptr] += coef_b * abs_vars_rhs[iptr];
          }
        }
      }

      // recv mesg
      MPI_Waitall(num_of_r_reqs, r_reqs, MPI_STATUS_IGNORE);
      if (istage != num_rk_stages-1) {
        fd_blk_unpack_mesg(rbuff, w_tmp, w3d_num_of_vars,
                           ni1,ni2,nj1,nj2,nk1,nk2,nx,ny,
                           siz_line,siz_slice,siz_volume);
      } else {
        fd_blk_unpack_mesg(rbuff, w_end, w3d_num_of_vars,
                           ni1,ni2,nj1,nj2,nk1,nk2,nx,ny,
                           siz_line,siz_slice,siz_volume);
      }

      MPI_Startall(num_of_r_reqs, r_reqs);
      MPI_Waitall(num_of_s_reqs, s_reqs, MPI_STATUS_IGNORE);
    } // RK stages

    // QC
    if (qc_check_nan_num_of_step >0  && (it % qc_check_nan_num_of_step) == 0) {
      if (myid==0 && verbose>10) fprintf(stdout,"-> check value nan\n");
        //wf_el_1st_check_value(w_end);
    }

    // save results
    //-- sta by interp
    for (int n=0; n<num_of_sta; n++)
    {
      int iptr = sta_loc_indx[0] + sta_loc_indx[1] * siz_line + sta_loc_indx[2] * siz_slice;
      // need to implement interp, now just take value
      for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
        int iptr_sta = (n * w3d_num_of_vars + ivar) * nt_total + it;
        sta_seismo[iptr_sta] = w_end[ivar*siz_volume + iptr];
      }
    }
    //-- point values
    for (int n=0; n<num_of_point; n++)
    {
      int iptr = point_loc_indx[0] + point_loc_indx[1] * siz_line + point_loc_indx[2] * siz_slice;
      for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
        int iptr_point = (n * w3d_num_of_vars + ivar) * nt_total + it;
        point_seismo[iptr_point] = w_end[ivar*siz_volume + iptr];
      }
    }
    //-- slice x, use w_rhs as buff
    int max_used_rhs = 0;
    for (int n=0; n<num_of_slice_x; n++) {
      size_t startp[] = { it, 0, 0 };
      size_t countp[] = { 1, nk, nj};
      size_t start_tdim = it;
      for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
        int iptr_var = w3d_pos[ivar];
        int iptr_slice = 0;
        for (int k=nk1; k<=nk2; k++) {
          for (int j=nj1; j<=nj2; j++) {
            int i = slice_x_indx[n];
            int iptr = i + j * siz_line + k * siz_slice + iptr_var;
            w_rhs[iptr_slice] = w_end[iptr];
            iptr_slice++;
          }
        }
        nc_put_var1_float(ncid_slx[n],timeid_slx,&start_tdim,&t_next);
        nc_put_vara_float(ncid_slx[n],varid_slx[ivar],startp,countp,w_rhs);
        max_used_rhs = (iptr_slice > max_used_rhs) ? iptr_slice : max_used_rhs;
      }
    }
    // slice y
    for (int n=0; n<num_of_slice_y; n++) {
      size_t startp[] = { it, 0, 0 };
      size_t countp[] = { 1, nk, ni};
      size_t start_tdim = it;
      for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
        int iptr_var = w3d_pos[ivar];
        int iptr_slice = 0;
        for (int k=nk1; k<=nk2; k++) {
          int j = slice_y_indx[n];
          for (int i=ni1; i<=ni2; i++) {
            int iptr = i + j * siz_line + k * siz_slice + iptr_var;
            w_rhs[iptr_slice] = w_end[iptr];
            iptr_slice++;
          }
        }
        nc_put_var1_float(ncid_sly[n],timeid_sly,&start_tdim,&t_next);
        nc_put_vara_float(ncid_sly[n],varid_sly[ivar],startp,countp,w_rhs);
        max_used_rhs = (iptr_slice > max_used_rhs) ? iptr_slice : max_used_rhs;
      }
    }
    // slice z
    for (int n=0; n<num_of_slice_z; n++) {
      size_t startp[] = { it, 0, 0 };
      size_t countp[] = { 1, nj, ni};
      size_t start_tdim = it;
      for (int ivar=0; ivar<w3d_num_of_vars; ivar++) {
        int iptr_var = w3d_pos[ivar];
        int k = slice_z_indx[n];
        int iptr_slice = 0;
        for (int j=nj1; j<=nj2; j++) {
          for (int i=ni1; i<=ni2; i++) {
            int iptr = i + j * siz_line + k * siz_slice + iptr_var;
            w_rhs[iptr_slice] = w_end[iptr];
            iptr_slice++;
          }
        }
        nc_put_var1_float(ncid_slz[n],timeid_slz,&start_tdim,&t_next);
        nc_put_vara_float(ncid_slz[n],varid_slz[ivar],startp,countp,w_rhs);
        max_used_rhs = (iptr_slice > max_used_rhs) ? iptr_slice : max_used_rhs;
      }
    }
    // snapshot
    for (int n=0; n<num_of_snap; n++)
    {
      int *cur_snap_info = snap_info + n * FD_SNAP_INFO_SIZE;
      int snap_i1  = cur_snap_info[FD_SNAP_INFO_I1];
      int snap_j1  = cur_snap_info[FD_SNAP_INFO_J1];
      int snap_k1  = cur_snap_info[FD_SNAP_INFO_K1];
      int snap_ni  = cur_snap_info[FD_SNAP_INFO_NI];
      int snap_nj  = cur_snap_info[FD_SNAP_INFO_NJ];
      int snap_nk  = cur_snap_info[FD_SNAP_INFO_NK];
      int snap_di  = cur_snap_info[FD_SNAP_INFO_DI];
      int snap_dj  = cur_snap_info[FD_SNAP_INFO_DJ];
      int snap_dk  = cur_snap_info[FD_SNAP_INFO_DK];

      int snap_it1 = cur_snap_info[FD_SNAP_INFO_IT1];
      int snap_dit = cur_snap_info[FD_SNAP_INFO_DIT];
      int snap_it_mod = (it - snap_it1) % snap_dit;
      int snap_it_num = (it - snap_it1) / snap_dit;
      int snap_nt_total = (nt_total - snap_it1) / snap_dit;
      int snap_out_V = cur_snap_info[FD_SNAP_INFO_VEL];
      int snap_out_T = cur_snap_info[FD_SNAP_INFO_STRESS];
      int snap_out_E = cur_snap_info[FD_SNAP_INFO_STRAIN];

      int snap_max_num = snap_ni * snap_nj * snap_nk;

      if (it>=snap_it1 && snap_it_num<=snap_nt_total && snap_it_mod==0)
      {
        size_t startp[] = { snap_cur_it[n], 0, 0, 0 };
        size_t countp[] = { 1, snap_nk, snap_nj, snap_ni };
        size_t start_tdim = snap_cur_it[n];

        nc_put_var1_float(ncid_snap[n],timeid_snap[n],&start_tdim,&t_next);

        // vel
        if (snap_out_V==1)
        {
          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Vx*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_vel[n*FD_NDIM+0],startp,countp,w_rhs);

          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Vy*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_vel[n*FD_NDIM+1],startp,countp,w_rhs);

          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Vz*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_vel[n*FD_NDIM+2],startp,countp,w_rhs);
        }
        if (snap_out_T==1)
        {
          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Txx*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_T[n*FD_NDIM_2+0],startp,countp,w_rhs);

          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Tyy*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_T[n*FD_NDIM_2+1],startp,countp,w_rhs);
          
          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Tzz*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_T[n*FD_NDIM_2+2],startp,countp,w_rhs);
          
          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Txz*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_T[n*FD_NDIM_2+3],startp,countp,w_rhs);

          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Tyz*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_T[n*FD_NDIM_2+4],startp,countp,w_rhs);

          sv_eliso1st_curv_macdrp_snap_buff(w_end+WF_EL_1ST_SEQ_Txy*siz_volume,
                siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                snap_dj,snap_k1,snap_nk,snap_dk,w_rhs);
          nc_put_vara_float(ncid_snap[n],varid_snap_T[n*FD_NDIM_2+5],startp,countp,w_rhs);
        }
        if (snap_out_E==1)
        {
          // need to implement
        }
        //for (int ivar=0; ivar<w3d_num_of_vars; ivar++)
        //{
        //  if (snap_save_cmp[n*w3d_num_of_vars+ivar]==1)
        //  {
        //    int iptr_var = w3d_pos[ivar];
        //    int iptr_snap=0;
        //    for (int n3=0; n3<snap_nk; n3++)
        //    {
        //      int k = cur_snap_info[2] + n3 * cur_snap_info[8];
        //      for (int n2=0; n2<snap_nj; n2++)
        //      {
        //        int j = cur_snap_info[1] + n2 * cur_snap_info[7];
        //        for (int n1=0; n1<snap_ni; n1++)
        //        {
        //          int i = cur_snap_info[0] + n1 * cur_snap_info[6];
        //          int iptr = i + j * siz_line + k * siz_slice + iptr_var;
        //          w_rhs[iptr_snap] = w_end[iptr];
        //          iptr_snap++;
        //        }
        //      }
        //    }
        //    size_t startp[] = { snap_cur_it[n], 0, 0, 0 };
        //    size_t countp[] = { 1, snap_nk, snap_nj, snap_ni };
        //    nc_put_vara_float(ncid_snap[n],varid_snap[n*w3d_num_of_vars+ivar],startp,countp,w_rhs);
        //  }
        //}
        max_used_rhs = (snap_max_num > max_used_rhs) ? snap_max_num : max_used_rhs;
        snap_cur_it[n] += 1;
      }
    }
    // zero temp used w_rsh
    for (int i=0; i<max_used_rhs; i++) w_rhs[i]=0.0;

    // debug output
    if (output_all==1) {
        io_build_fname_time(output_dir,"w3d",".nc",myid2,it,ou_file);
        io_var3d_export_nc(ou_file,
                           w_end,
                           w3d_pos,
                           w3d_name,
                           w3d_num_of_vars,
                           coord_name,
                           nx,
                           ny,
                           nz);
    }

    // swap w_pre and w_end, avoid copying
    w_cur = w_pre; w_pre = w_end; w_end = w_cur;
    abs_vars_cur = abs_vars_pre; abs_vars_pre = abs_vars_end; abs_vars_end = abs_vars_cur;

  } // time loop

  // postproc
  if (force_vec_value ) free(force_vec_value );
  if (moment_ten_value) free(moment_ten_value);

  // close nc
  for (int n=0; n<num_of_slice_x; n++) {
    nc_close(ncid_slx[n]);
  }
  for (int n=0; n<num_of_slice_y; n++) {
    nc_close(ncid_sly[n]);
  }
  for (int n=0; n<num_of_slice_z; n++) {
    nc_close(ncid_slz[n]);
  }
  for (int n=0; n<num_of_snap; n++)
  {
    nc_close(ncid_snap[n]);
  }
}

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_onestage(
    float *restrict w_cur, float *restrict rhs, 
    float *restrict g3d, float *restrict m3d,
    // grid size
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    int ni, int nj, int nk, int nx, int ny, int nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    int *restrict boundary_itype,
    int              abs_itype,
    int    *restrict abs_num_of_layers,
    int *restrict abs_indx,
    int *restrict abs_coefs_facepos0,
    float  *restrict abs_coefs,
    int           abs_vars_size_per_level,
    int *restrict abs_vars_volsiz,
    int *restrict abs_vars_facepos0, //
    float  *restrict abs_vars_cur,
    float  *restrict abs_vars_rhs,
    float *matVx2Vz, float *matVy2Vz, //
    // source term
    int num_of_force,
    int *restrict force_info,
    float *restrict force_vec_value, // only for cur stage, size: num_of_force
    int   *restrict force_ext_indx,
    float *restrict force_ext_coef,
    int             num_of_moment,
    int   *restrict moment_info, // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
    float *restrict moment_ten_value,
    int   *restrict moment_ext_indx,
    float *restrict moment_ext_coef,
    // include different order/stentil
    int fdx_max_half_len, int fdy_max_half_len,
    int fdz_max_len, int fdz_num_surf_lay,
    int **restrict fdx_all_info, int *restrict fdx_all_indx,float *restrict fdx_all_coef,
    int **restrict fdy_all_info, int *restrict fdy_all_indx,float *restrict fdy_all_coef,
    int **restrict fdz_all_info, int *restrict fdz_all_indx,float *restrict fdz_all_coef,
    const int myid, const int verbose)
{
  int i,j,k;

  // local pointer get each vars
  float *restrict Vx    = w_cur + WF_EL_1ST_SEQ_Vx    * siz_volume;
  float *restrict Vy    = w_cur + WF_EL_1ST_SEQ_Vy    * siz_volume;
  float *restrict Vz    = w_cur + WF_EL_1ST_SEQ_Vz    * siz_volume;
  float *restrict Txx   = w_cur + WF_EL_1ST_SEQ_Txx   * siz_volume;
  float *restrict Tyy   = w_cur + WF_EL_1ST_SEQ_Tyy   * siz_volume;
  float *restrict Tzz   = w_cur + WF_EL_1ST_SEQ_Tzz   * siz_volume;
  float *restrict Txz   = w_cur + WF_EL_1ST_SEQ_Txz   * siz_volume;
  float *restrict Tyz   = w_cur + WF_EL_1ST_SEQ_Tyz   * siz_volume;
  float *restrict Txy   = w_cur + WF_EL_1ST_SEQ_Txy   * siz_volume;
  float *restrict hVx   = rhs + WF_EL_1ST_SEQ_Vx    * siz_volume;
  float *restrict hVy   = rhs + WF_EL_1ST_SEQ_Vy    * siz_volume;
  float *restrict hVz   = rhs + WF_EL_1ST_SEQ_Vz    * siz_volume;
  float *restrict hTxx  = rhs + WF_EL_1ST_SEQ_Txx   * siz_volume;
  float *restrict hTyy  = rhs + WF_EL_1ST_SEQ_Tyy   * siz_volume;
  float *restrict hTzz  = rhs + WF_EL_1ST_SEQ_Tzz   * siz_volume;
  float *restrict hTxz  = rhs + WF_EL_1ST_SEQ_Txz   * siz_volume;
  float *restrict hTyz  = rhs + WF_EL_1ST_SEQ_Tyz   * siz_volume;
  float *restrict hTxy  = rhs + WF_EL_1ST_SEQ_Txy   * siz_volume;
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
  float *restrict lam3d = m3d   + MD_EL_ISO_SEQ_LAMBDA * siz_volume;
  float *restrict  mu3d = m3d   + MD_EL_ISO_SEQ_MU     * siz_volume;
  float *restrict slw3d = m3d   + MD_EL_ISO_SEQ_RHO    * siz_volume;

  // fd op
  int           fdx_inn_len;
  int           fdx_inn_pos;
  int *restrict fdx_inn_indx;
  float  *restrict fdx_inn_coef;
  int           fdy_inn_len;
  int           fdy_inn_pos;
  int *restrict fdy_inn_indx;
  float  *restrict fdy_inn_coef;
  int           fdz_inn_len;
  int           fdz_inn_pos;
  int *restrict fdz_inn_indx;
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
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        myid, verbose);

    // velocity: vlow
    sv_eliso1st_curv_macdrp_rhs_vlow_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        lam3d, mu3d, slw3d,matVx2Vz,matVy2Vz,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
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
      // cfs-pml modified for timg: can't used, double corr with main cfspml
      //sv_eliso1st_curv_macdrp_rhs_cfspml_timg_z2(Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
      //                                           xi_x, xi_y, xi_z,
      //                                           et_x, et_y, et_z, zt_x, zt_y, zt_z,
      //                                           jac3d, slw3d, nk2, siz_line,siz_slice,
      //                                           fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
      //                                           fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
      //                                           fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
      //                                           boundary_itype, abs_num_of_layers, abs_indx,
      //                                           abs_coefs_facepos0, abs_coefs,
      //                                           abs_vars_volsiz, abs_vars_facepos0,
      //                                           abs_vars_cur, abs_vars_rhs,
      //                                           myid, verbose);
      
      // cfs-pml modified for vfree; vlow at points below surface doesn't affect pml
      sv_eliso1st_curv_macdrp_rhs_cfspml_vfree_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                                  zt_x, zt_y, zt_z,
                                                  lam3d, mu3d, matVx2Vz, matVy2Vz,
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
                                    num_of_force, force_info, force_vec_value,
                                    force_ext_indx,force_ext_coef,
                                    num_of_moment, moment_info, moment_ten_value,
                                    moment_ext_indx,moment_ext_coef,
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
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  float *restrict Vx_ptr;
  float *restrict Vy_ptr;
  float *restrict Vz_ptr;
  float *restrict Txx_ptr;
  float *restrict Txy_ptr;
  float *restrict Txz_ptr;
  float *restrict Tyy_ptr;
  float *restrict Tzz_ptr;
  float *restrict Tyz_ptr;

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

#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
  /*
  for (int i=0; i < fdx_len; i++) {
    fprintf(stdout," %d", fdx_indx [i]);
  }
  fprintf(stdout,"\n");
  for (int i=0; i < fdx_len; i++) {
    fprintf(stdout," %f", fdx_coef [i]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);

  for (int j=0; j < fdy_len; j++) {
    fprintf(stdout," %d", fdy_indx [j]);
  }
  fprintf(stdout,"\n");
  for (int j=0; j < fdy_len; j++) {
    fprintf(stdout," %f", fdy_coef [j]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);

  for (int k=0; k < fdz_len; k++) {
    fprintf(stdout," %d", fdz_indx [k]);
  }
  fprintf(stdout,"\n");
  for (int k=0; k < fdz_len; k++) {
    fprintf(stdout," %f", fdz_coef [k]);
  }
  fprintf(stdout,"\n");
  fflush(stdout);
  */
#endif

  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
    if (myid==0 && verbose>900) fprintf(stdout,"----> nk1=%d,nk2=%d,k=%d\n", nk1,nk2,k);
#endif
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
      if (myid==0 && verbose>910) fprintf(stdout,"-----> nj1=%d,nj2=%d,j=%d\n", nj1,nj2,j);
#endif
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
        if (myid==0 && verbose>920) fprintf(stdout,"-----> ni1=%d,ni2=%d,i=%d\n", ni1,ni2,i);
#endif
        Vx_ptr = Vx + iptr;
        Vy_ptr = Vy + iptr;
        Vz_ptr = Vz + iptr;
        Txx_ptr = Txx + iptr;
        Tyy_ptr = Tyy + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;
        Tyz_ptr = Tyz + iptr;
        Txy_ptr = Txy + iptr;

        // Vx derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyVx, Vx_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxVy, Vy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyVy, Vy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzVy, Vy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyVz, Vz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txx derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxTxx, Txx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyTxx, Txx_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzTxx, Txx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyy derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxTyy, Tyy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyTyy, Tyy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzTyy, Tyy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tzz derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxTzz, Tzz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyTzz, Tzz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyTxz, Txz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyz derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxTyz, Tyz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyTyz, Tyz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzTyz, Tyz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txy derivatives
        M_FD_SHIFT_PTR_UNLOOP5(DxTxy, Txy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DyTxy, Txy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_UNLOOP5(DzTxy, Txy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

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
        hTxx[iptr] =  lam2mu * ( xix*DxVx  +etx*DyVx + ztx*DzVx)
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
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // loop var for fd
  int n_fd; // loop var for fd

  // local var
  float DxTx,DyTy,DzTz;
  float slwjac;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  // to save traction and other two dir force var
  float vecxi[fdx_len];
  float vecet[fdy_len];
  float veczt[fdz_len];
  int n, iptr4vec;

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
  int k_min = nk2 - fdz_indx[fdz_len-1];

  // point affected by timg
  for (size_t k=k_min; k <= nk2; k++)
  {
    // k corresponding to 0 index of the fd op

    // index of free surface
    int n_free = nk2 - k - fdz_indx[0]; // first indx is negative

    // for 1d index
    size_t iptr_k = k * siz_slice;

    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
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
          slwjac = slw3d[iptr] / jac3d[iptr];

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
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++)
          {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            //int n_img = n_free-(n-n_free);
            iptr4vec = iptr + n_img * siz_slice;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txx[iptr4vec]
                                           + zt_y[iptr4vec] * Txy[iptr4vec]
                                           + zt_z[iptr4vec] * Txz[iptr4vec] );
            //veczt[n] = -veczt[n_free-(n-n_free)];
          }

          // deri
          M_FD_NOINDX(DxTx, vecxi, fdx_len, lfdx_coef, n_fd);
          M_FD_NOINDX(DyTy, vecet, fdy_len, lfdy_coef, n_fd);
          M_FD_NOINDX(DzTz, veczt, fdz_len, lfdz_coef, n_fd);

          hVx[iptr] = ( DxTx+DyTy+DzTz ) * slwjac;
#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
            /*
            if (hVx[iptr] > 1.0e-25) {
              fprintf(stderr, "WARNING: nonzero value for zero input:\n");
              fprintf(stderr, "  i=%d,j=%d,k=%d: DxTx=%f, DyTy=%f, DzTz=%f\n",
                      i,j,k,DxTx,DyTy,DzTz);
              fprintf(stderr, "  %e,%e,%e,%e,%e\n",
                                  veczt[0],
                                  veczt[1],
                                  veczt[2],
                                  veczt[3],
                                  veczt[4]);
              fprintf(stderr, "  n_free=%d,siz_slice=%d,fdz_len=%d,fdz_indx=%d,%d,%d,%d,%d\n",
                          n_free,siz_slice,fdz_len,
                          fdz_indx[0],
                          fdz_indx[1],
                          fdz_indx[2],
                          fdz_indx[3],
                          fdz_indx[4]);

              for (n=0; n<n_free; n++) {
                iptr4vec = iptr + fdz_indx[n] * siz_slice;
                fprintf(stderr,"   n=%d,jac=%e,zt_x=%e,Txx=%e,zt_y=%e,Txy=%e,zt_z=%e,Txz=%e\n",
                           jac3d[iptr4vec] ,    zt_x[iptr4vec] , Txx[iptr4vec],
                                                zt_y[iptr4vec] , Txy[iptr4vec],
                                                zt_z[iptr4vec] , Txz[iptr4vec]);
              }

              fflush(stdout);
              exit(2);
            }
            */
#endif

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
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            iptr4vec = iptr + n_img * siz_slice;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txy[iptr4vec]
                                           + zt_y[iptr4vec] * Tyy[iptr4vec]
                                           + zt_z[iptr4vec] * Tyz[iptr4vec] );
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
          veczt[n_free] = 0.0;

          // above surface -> mirror
          for (n=n_free+1; n<fdz_len; n++) {
            int n_img = fdz_indx[n] - 2*(n-n_free);
            iptr4vec = iptr + n_img * siz_slice;
            veczt[n] = -jac3d[iptr4vec] * (  zt_x[iptr4vec] * Txz[iptr4vec]
                                           + zt_y[iptr4vec] * Tyz[iptr4vec]
                                           + zt_z[iptr4vec] * Tzz[iptr4vec] );
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
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_num_surf_lay, int fdz_max_len,
    int  **restrict fdz_all_info,
    int   *restrict fdz_all_indx,
    float *restrict fdz_all_coef,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];

  // allocate max_len because fdz may have different lens
  float  lfdz_coef [fdz_max_len];
  int lfdz_shift[fdz_max_len];

  // local var
  int i,j,k;
  int n_fd; // loop var for fd
  int fdz_len;

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
  //for (size_t n=0; n < 1; n++)
  for (size_t n=0; n < fdz_num_surf_lay; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // get pos and len for this point
    int  lfdz_pos0 = fdz_all_info[n][FD_INFO_POS_OF_INDX];
    int  lfdz_len  = fdz_all_info[n][FD_INFO_LENGTH_TOTAL];
    // point to indx/coef for this point
    int   *p_fdz_indx  = fdz_all_indx+lfdz_pos0;
    float *p_fdz_coef  = fdz_all_coef+lfdz_pos0;
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
          size_t ij = (i + j * siz_line)*9;
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
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    int *restrict abs_indx, // pml range of each face
    int *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    int *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    int *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose)
{
  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // val on point
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
  float hVx_rhs,hVy_rhs,hVz_rhs;
  float hTxx_rhs,hTyy_rhs,hTzz_rhs,hTxz_rhs,hTyz_rhs,hTxy_rhs;

  // local
  int i,j,k;
  int iptr, iptr_j, iptr_k, iptr_a;
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
    int abs_ni1 = abs_indx[iface*FD_NDIM_2+0];
    int abs_ni2 = abs_indx[iface*FD_NDIM_2+1];
    int abs_nj1 = abs_indx[iface*FD_NDIM_2+2];
    int abs_nj2 = abs_indx[iface*FD_NDIM_2+3];
    int abs_nk1 = abs_indx[iface*FD_NDIM_2+4];
    int abs_nk2 = abs_indx[iface*FD_NDIM_2+5];

#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
    //fprintf(stdout," iface=%d,ni1=%d,ni2=%d,nj1=%d,nj2=%d,nk1=%d,nk2=%d\n",
    //        iface,abs_ni1,abs_ni2,abs_nj1,abs_nj2,abs_nk1,abs_nk2);
    //fflush(stdout);
#endif

    // get coef for this face
    float *restrict ptr_coef_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    int pos_cur_face = abs_vars_facepos0[iface];
    int pml_volsiz   = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vx  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vy  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vz  * pml_volsiz;
    float *restrict pml_Txx = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Txx * pml_volsiz;
    float *restrict pml_Tyy = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Tyy * pml_volsiz;
    float *restrict pml_Tzz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Tzz * pml_volsiz;
    float *restrict pml_Txz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Txz * pml_volsiz;
    float *restrict pml_Tyz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Tyz * pml_volsiz;
    float *restrict pml_Txy = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Txy * pml_volsiz;
    float *restrict pml_hVx  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vx  * pml_volsiz;
    float *restrict pml_hVy  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vy  * pml_volsiz;
    float *restrict pml_hVz  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vz  * pml_volsiz;
    float *restrict pml_hTxx = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Txx * pml_volsiz;
    float *restrict pml_hTyy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Tyy * pml_volsiz;
    float *restrict pml_hTzz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Tzz * pml_volsiz;
    float *restrict pml_hTxz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Txz * pml_volsiz;
    float *restrict pml_hTyz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Tyz * pml_volsiz;
    float *restrict pml_hTxy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Txy * pml_volsiz;

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
            pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
            pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
            pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
            pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
            pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
            pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
            pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
            pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
            pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

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
            pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
            pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
            pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
            pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
            pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
            pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
            pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
            pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
            pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

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
            pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
            pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
            pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
            pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
            pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
            pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
            pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
            pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
            pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
  }
}

// considering timg free surface in cfs-pml: conflict with main cfspml, need to revise
//  in the future if required
/*
void
sv_eliso1st_curv_macdrp_rhs_cfspml_timg_z2(
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict jac3d, float *restrict slw3d,
    int nk2, size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    int *restrict abs_indx, // pml range of each face
    int *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    int *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    int *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose)
{
  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];
  
  // val on point
  float slw, slwjac;

  // local
  float DxCx, DyCy;
  float coef_A, coef_B, coef_D;
  float hVx_rhs,hVy_rhs,hVz_rhs;

  // to save traction and other two dir force var
  float vecxi[fdz_len];
  float vecet[fdz_len];
  float veczt[fdz_len];
  int n, iptr4vec;
  int iptr, iptr_j, iptr_k, iptr_a;

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
  int k_min = nk2 - fdz_indx[fdz_len-1];

  // loop x/y face
  //for (int iface=0; iface<2; iface++)
  //for (int iface=3; iface<4; iface++)
  for (int iface=0; iface<4; iface++)
  {
    // skip to next face if not cfspml
    if (boundary_itype[iface] != FD_BOUNDARY_TYPE_CFSPML) continue;

    // get index into local var
    int abs_ni1 = abs_indx[iface*FD_NDIM_2+0];
    int abs_ni2 = abs_indx[iface*FD_NDIM_2+1];
    int abs_nj1 = abs_indx[iface*FD_NDIM_2+2];
    int abs_nj2 = abs_indx[iface*FD_NDIM_2+3];
    int abs_nk1 = abs_indx[iface*FD_NDIM_2+4];
    int abs_nk2 = abs_indx[iface*FD_NDIM_2+5];

    // max of two values
    int abs_nk1_kmin = (k_min > abs_nk1 ) ? k_min : abs_nk1;

    // get coef for this face
    float *restrict ptr_coef_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    int pos_cur_face = abs_vars_facepos0[iface];
    int pml_volsiz   = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vx  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vy  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vz  * pml_volsiz;
    float *restrict pml_hVx  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vx  * pml_volsiz;
    float *restrict pml_hVy  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vy  * pml_volsiz;
    float *restrict pml_hVz  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vz  * pml_volsiz;

    // for each dim
    if (iface < 2) // x direction
    {
      iptr_a = (abs_nk1_kmin-abs_nk1) * (abs_nj2-abs_nj1+1) * (abs_ni2-abs_ni1+1);
      for (int k=abs_nk1_kmin; k<=abs_nk2; k++)
      {
        iptr_k = k * siz_slice;
        for (int j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;
          for (int i=abs_ni1; i<=abs_ni2; i++)
          {
            // pml coefs
            int abs_i = i - abs_ni1;
            coef_D = ptr_coef_D[abs_i];
            coef_A = ptr_coef_A[abs_i];
            coef_B = ptr_coef_B[abs_i];

            // slowness and jac
            slwjac = slw3d[iptr] / jac3d[iptr];

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
            pml_hVx[iptr_a] = coef_D * hVx_rhs - coef_A * pml_Vx[iptr_a];

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
            pml_hVy[iptr_a] = coef_D * hVy_rhs - coef_A * pml_Vy[iptr_a];

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
            pml_hVz[iptr_a] = coef_D * hVz_rhs - coef_A * pml_Vz[iptr_a];

            // incr index
            iptr   += 1;
            iptr_a += 1;
          } // i
        } // j
      } // k
    }
    else // y direction
    {
      iptr_a = (abs_nk1_kmin-abs_nk1) * (abs_nj2-abs_nj1+1) * (abs_ni2-abs_ni1+1);
      for (int k=abs_nk1_kmin; k<=abs_nk2; k++)
      {
        iptr_k = k * siz_slice;
        for (int j=abs_nj1; j<=abs_nj2; j++)
        {
          iptr_j = iptr_k + j * siz_line;
          iptr = iptr_j + abs_ni1;

          // pml coefs
          int abs_j = j - abs_nj1;
          coef_D = ptr_coef_D[abs_j];
          coef_A = ptr_coef_A[abs_j];
          coef_B = ptr_coef_B[abs_j];

          for (int i=abs_ni1; i<=abs_ni2; i++)
          {
            // slowness and jac
            slwjac = slw3d[iptr] / jac3d[iptr];

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
            pml_hVx[iptr_a] = coef_D * hVx_rhs - coef_A * pml_Vx[iptr_a];

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
            pml_hVy[iptr_a] = coef_D * hVy_rhs - coef_A * pml_Vy[iptr_a];

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
            pml_hVz[iptr_a] = coef_D * hVz_rhs - coef_A * pml_Vz[iptr_a];

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

// considering velocity free surface condition in cfs-pml
void
sv_eliso1st_curv_macdrp_rhs_cfspml_vfree_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    int nk2, size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    int *restrict abs_indx, // pml range of each face
    int *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    int *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    int *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs,
    const int myid, const int verbose)
{
  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];

  // val on point
  float lam,mu,lam2mu;

  // local
  int i,j,k;
  int iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D;
  float DxVx,DxVy,DxVz;
  float DyVx,DyVy,DyVz;
  float Dx_DzVx,Dy_DzVx,Dx_DzVy,Dy_DzVy,Dx_DzVz,Dy_DzVz;
  float hTxx_rhs,hTyy_rhs,hTzz_rhs,hTxz_rhs,hTyz_rhs,hTxy_rhs;

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
    int abs_ni1 = abs_indx[iface*FD_NDIM_2+0];
    int abs_ni2 = abs_indx[iface*FD_NDIM_2+1];
    int abs_nj1 = abs_indx[iface*FD_NDIM_2+2];
    int abs_nj2 = abs_indx[iface*FD_NDIM_2+3];
    int abs_nk1 = abs_indx[iface*FD_NDIM_2+4];
    int abs_nk2 = abs_indx[iface*FD_NDIM_2+5];

    // get coef for this face
    float *restrict ptr_coef_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    int pos_cur_face = abs_vars_facepos0[iface];
    int pml_volsiz = abs_vars_volsiz[iface];
    float *restrict pml_Txx = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Txx * pml_volsiz;
    float *restrict pml_Tyy = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Tyy * pml_volsiz;
    float *restrict pml_Tzz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Tzz * pml_volsiz;
    float *restrict pml_Txz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Txz * pml_volsiz;
    float *restrict pml_Tyz = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Tyz * pml_volsiz;
    float *restrict pml_Txy = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Txy * pml_volsiz;
    float *restrict pml_hTxx = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Txx * pml_volsiz;
    float *restrict pml_hTyy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Tyy * pml_volsiz;
    float *restrict pml_hTzz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Tzz * pml_volsiz;
    float *restrict pml_hTxz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Txz * pml_volsiz;
    float *restrict pml_hTyz = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Tyz * pml_volsiz;
    float *restrict pml_hTxy = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Txy * pml_volsiz;

    // for each dim
    if (iface < 2) // x direction
    {
      if (nk2 > abs_nk2) continue;

      // set k to free surface
      k = nk2;
      iptr_k = k * siz_slice;
      iptr_a = (k-abs_nk1) * (abs_nj2-abs_nj1+1) * (abs_ni2-abs_ni1+1);

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
          float ztx = zt_x[iptr];
          float zty = zt_y[iptr];
          float ztz = zt_z[iptr];

          // medium
          lam = lam3d[iptr];
          mu  =  mu3d[iptr];
          lam2mu = lam + 2.0 * mu;

          // Vx derivatives
          M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
          // Vy derivatives
          M_FD_SHIFT(DxVy, Vy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
          // Vz derivatives
          M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

          // zeta derivatives
          int ij = (i + j * siz_line)*9;
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
          hTxx_rhs =    lam2mu * (            ztx*Dx_DzVx)
                      + lam    * (            zty*Dx_DzVy
                                  +           ztz*Dx_DzVz);

          hTyy_rhs =   lam2mu * (            zty*Dx_DzVy)
                      +lam    * (            ztx*Dx_DzVx
                                            +ztz*Dx_DzVz);

          hTzz_rhs =   lam2mu * (            ztz*Dx_DzVz)
                      +lam    * (            ztx*Dx_DzVx
                                            +zty*Dx_DzVy);

          hTxy_rhs = mu *(
                       zty*Dx_DzVx + ztx*Dx_DzVy
                      );
          hTxz_rhs = mu *(
                       ztz*Dx_DzVx + ztx*Dx_DzVz
                      );
          hTyz_rhs = mu *(
                       ztz*Dx_DzVy + zty*Dx_DzVz
                      );

          //// make corr to Hooke's equatoin
          //hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs - coef_B * pml_Txx[iptr_a];
          //hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
          //hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
          //hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs - coef_B * pml_Txz[iptr_a];
          //hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
          //hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs - coef_B * pml_Txy[iptr_a];

          //// aux var
          ////   a1 = alpha + d / beta, dealt in abs_set_cfspml
          //pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
          //pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
          //pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
          //pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
          //pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
          //pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

          // make corr to Hooke's equatoin
          hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
          hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs;
          hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
          hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;
          hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs;
          hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs;

          // aux var
          //   a1 = alpha + d / beta, dealt in abs_set_cfspml
          pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
          pml_hTyy[iptr_a] += coef_D * hTyy_rhs;
          pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
          pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
          pml_hTyz[iptr_a] += coef_D * hTyz_rhs;
          pml_hTxy[iptr_a] += coef_D * hTxy_rhs;

          // incr index
          iptr   += 1;
          iptr_a += 1;
        } // i
      } // j
    }
    else // y direction
    {
      if (nk2 > abs_nk2) continue;

      // set k to free surface
      k = nk2;
      iptr_k = k * siz_slice;
      iptr_a = (k-abs_nk1) * (abs_nj2-abs_nj1+1) * (abs_ni2-abs_ni1+1);

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
          float ztx = zt_x[iptr];
          float zty = zt_y[iptr];
          float ztz = zt_z[iptr];

          // medium
          lam = lam3d[iptr];
          mu  =  mu3d[iptr];
          lam2mu = lam + 2.0 * mu;

          // Vx derivatives
          M_FD_SHIFT(DyVx, Vx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
          // Vy derivatives
          M_FD_SHIFT(DyVy, Vy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
          // Vz derivatives
          M_FD_SHIFT(DyVz, Vz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

          // zeta derivatives
          int ij = (i + j * siz_line)*9;
          Dy_DzVx = matVy2Vz[ij+3*0+0] * DyVx
                  + matVy2Vz[ij+3*0+1] * DyVy
                  + matVy2Vz[ij+3*0+2] * DyVz;

          Dy_DzVy = matVy2Vz[ij+3*1+0] * DyVx
                  + matVy2Vz[ij+3*1+1] * DyVy
                  + matVy2Vz[ij+3*1+2] * DyVz;

          Dy_DzVz = matVy2Vz[ij+3*2+0] * DyVx
                  + matVy2Vz[ij+3*2+1] * DyVy
                  + matVy2Vz[ij+3*2+2] * DyVz;

          hTxx_rhs =    lam2mu * (             ztx*Dy_DzVx);
                      + lam    * (             zty*Dy_DzVy
                                              +ztz*Dy_DzVz);

          hTyy_rhs =   lam2mu * (             zty*Dy_DzVy)
                      +lam    * (             ztx*Dy_DzVx
                                             +ztz*Dy_DzVz);

          hTzz_rhs =   lam2mu * (             ztz*Dy_DzVz)
                      +lam    * (             ztx*Dy_DzVx
                                             +zty*Dy_DzVy);

          hTxy_rhs = mu *(
                       zty*Dy_DzVx + ztx*Dy_DzVy
                      );
          hTxz_rhs = mu *(
                       ztz*Dy_DzVx + ztx*Dy_DzVz
                      );
          hTyz_rhs = mu *(
                       ztz*Dy_DzVy + zty*Dy_DzVz
                    );

          //// make corr to Hooke's equatoin
          //hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs - coef_B * pml_Txx[iptr_a];
          //hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
          //hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
          //hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs - coef_B * pml_Txz[iptr_a];
          //hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs - coef_B * pml_Tyz[iptr_a];
          //hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs - coef_B * pml_Txy[iptr_a];

          //// aux var
          ////   a1 = alpha + d / beta, dealt in abs_set_cfspml
          //pml_hTxx[iptr_a] = coef_D * hTxx_rhs - coef_A * pml_Txx[iptr_a];
          //pml_hTyy[iptr_a] = coef_D * hTyy_rhs - coef_A * pml_Tyy[iptr_a];
          //pml_hTzz[iptr_a] = coef_D * hTzz_rhs - coef_A * pml_Tzz[iptr_a];
          //pml_hTxz[iptr_a] = coef_D * hTxz_rhs - coef_A * pml_Txz[iptr_a];
          //pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
          //pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

          // make corr to Hooke's equatoin
          hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
          hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs;
          hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
          hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;
          hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs;
          hTxy[iptr] += (coef_B - 1.0) * hTxy_rhs;

          // aux var
          //   a1 = alpha + d / beta, dealt in abs_set_cfspml
          pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
          pml_hTyy[iptr_a] += coef_D * hTyy_rhs;
          pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
          pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
          pml_hTyz[iptr_a] += coef_D * hTyz_rhs;
          pml_hTxy[iptr_a] += coef_D * hTxy_rhs;

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
    int nk2, size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_num_surf_lay, int fdz_max_len, int **restrict fdz_all_info,
    int *restrict fdz_all_indx, float *restrict fdz_all_coef,
    int *restrict boundary_itype,
    int *restrict abs_num_of_layers,
    int *restrict abs_indx, // pml range of each face
    int *restrict abs_coefs_facepos0, // pos of 0 elem of each face
    float  *restrict abs_coefs,
    int *restrict abs_vars_volsiz, // size of single var in abs_blk_vars for each block
    int *restrict abs_vars_facepos0, // start pos of each blk in first level; other level similar
    float  *restrict abs_vars_cur, //
    float  *restrict abs_vars_rhs)
{
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int lfdz_shift[fdz_len];

  // loop var for fd
  int i,j,k;
  int n_fd; // loop var for fd
  int fdz_len;

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
    int abs_ni1 = abs_indx[iface][0];
    int abs_ni2 = abs_indx[iface][1];
    int abs_nj1 = abs_indx[iface][2];
    int abs_nj2 = abs_indx[iface][3];
    int abs_nk1 = abs_indx[iface][4];
    int abs_nk2 = abs_indx[iface][5];

    // get coef for this face
    float *restrict ptr_coef_A = abs_coefs + abs_coefs_facepos0[iface];
    float *restrict ptr_coef_B = ptr_coef_A + abs_num_of_layers[iface];
    float *restrict ptr_coef_D = ptr_coef_B + abs_num_of_layers[iface];

    // get pml vars
    int pos_cur_face = abs_vars_facepos0[iface];
    int pml_volsiz = abs_vars_volsiz[iface];
    float *restrict pml_Vx  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vx  * pml_volsiz;
    float *restrict pml_Vy  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vy  * pml_volsiz;
    float *restrict pml_Vz  = abs_vars_cur + pos_cur_face + WF_EL_1ST_SEQ_Vz  * pml_volsiz;
    float *restrict pml_hVx  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vx  * pml_volsiz;
    float *restrict pml_hVy  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vy  * pml_volsiz;
    float *restrict pml_hVz  = abs_vars_rhs + pos_cur_face + WF_EL_1ST_SEQ_Vz  * pml_volsiz;

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
        int  lfdz_pos0 = fdz_info[n][FD_INFO_POS_OF_INDX];
        int  lfdz_len  = fdz_info[n][FD_INFO_LENGTH_TOTAL];
        int  *p_fdz_indx  = fdz_all_indx[lfdz_pos0];
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
    int number_of_vars, int *restrict vars_pos,
    // grid size
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    int ni, int nj, int nk, int nx, int ny, int nz,
    // ablexp info
    int *restrict abs_indx, float *restrict abs_damp
    )
{
  int ierr = 0;

  int iptr,iptr_abs, iptr_var, iptr_k, iptr_j;
  size_t siz_line, siz_slice;
  int i,j,k,ivar;

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
    int number_of_vars, int *restrict vars_pos,
    // grid size
    int nx, int ny, int nz,
    // ablexp info
    int *restrict abs_indx, float *restrict abs_damp
    )
{
  int ierr = 0;

  int iptr,iptr_abs, iptr_var, iptr_k, iptr_j;
  size_t siz_line, siz_slice;
  int i,j,k,ivar;
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

void
sv_eliso1st_curv_macdrp_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    size_t siz_line, size_t siz_slice,
    int num_of_force,
    int *restrict force_info,
    float *restrict force_vec_value,
    int   *restrict force_ext_indx,
    float *restrict force_ext_coef,
    int             num_of_moment,
    int   *restrict moment_info,
    float *restrict moment_ten_value, // size: num_of_moment * 6
    int   *restrict moment_ext_indx,
    float *restrict moment_ext_coef,
    const int myid, const int verbose)
{
  // local var
  int si,sj,sk, iptr;

  // add force
  for (int n=0; n<num_of_force; n++)
  {
    int    *ptr_force_info = force_info + n * M_SRC_INFO_NVAL;
    // float  *ptr_force_stf = force_vec_value + ptr_force_info[M_SRC_INFO_SEQ_POS];
    float  *ptr_force_stf = force_vec_value + n*3;
    //-- add at single point
    //si = ptr_force_info[0];
    //sj = ptr_force_info[1];
    //sk = ptr_force_info[2];

    //iptr = si + sj * siz_line + sk * siz_slice;

    //float V = slw3d[iptr] / jac3d[iptr];

    //hVx[iptr] += ptr_force_val[ 0 ] * V;
    //hVy[iptr] += ptr_force_val[ 1 ] * V;
    //hVz[iptr] += ptr_force_val[ 2 ] * V;

    //-- add with smoothing
    int n_ext_max  = ptr_force_info[M_SRC_INFO_SEQ_NEXT_MAX];
    int n_ext_this = ptr_force_info[M_SRC_INFO_SEQ_NEXT_THIS];

    int   *ptr_force_ext_indx = force_ext_indx + n * n_ext_max;
    float *ptr_force_ext_coef = force_ext_coef + n * n_ext_max;
    for (int i_ext=0; i_ext<n_ext_this; i_ext++)
    {
      iptr = ptr_force_ext_indx[i_ext];

      float V = ptr_force_ext_coef[i_ext] * slw3d[iptr] / jac3d[iptr];
      hVx[iptr] += ptr_force_stf[ 0 ] * V;
      hVy[iptr] += ptr_force_stf[ 1 ] * V;
      hVz[iptr] += ptr_force_stf[ 2 ] * V;
    }
  }

  // add moment source
  for (int n=0; n<num_of_moment; n++)
  {
    int    *ptr_moment_info = moment_info + n * M_SRC_INFO_NVAL;
    //float  *ptr_moment_mrf  = moment_ten_value + ptr_moment_info[M_SRC_INFO_SEQ_POS];
    float  *ptr_moment_mrf  = moment_ten_value + n*6;

    //si = ptr_moment_info[0];
    //sj = ptr_moment_info[1];
    //sk = ptr_moment_info[2];
    //iptr = si + sj * siz_line + sk * siz_slice;

    int n_ext_max  = ptr_moment_info[M_SRC_INFO_SEQ_NEXT_MAX];
    int n_ext_this = ptr_moment_info[M_SRC_INFO_SEQ_NEXT_THIS];

    int   *ptr_moment_ext_indx = moment_ext_indx + n * n_ext_max;
    float *ptr_moment_ext_coef = moment_ext_coef + n * n_ext_max;
    for (int i_ext=0; i_ext<n_ext_this; i_ext++)
    {
      iptr = ptr_moment_ext_indx[i_ext];

      float rjac = ptr_moment_ext_coef[i_ext] / jac3d[iptr];
      hTxx[iptr] -= ptr_moment_mrf[ 0 ] * rjac;
      hTyy[iptr] -= ptr_moment_mrf[ 1 ] * rjac;
      hTzz[iptr] -= ptr_moment_mrf[ 2 ] * rjac;
      hTxz[iptr] -= ptr_moment_mrf[ 3 ] * rjac;
      hTyz[iptr] -= ptr_moment_mrf[ 4 ] * rjac;
      hTxy[iptr] -= ptr_moment_mrf[ 5 ] * rjac;
    }
  }
}

/*******************************************************************************
 * related functions
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_vel_dxy2dz(
    float *restrict g3d,
    float *restrict m3d,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    float **restrict p_matVx2Vz,
    float **restrict p_matVy2Vz,
    const int myid, const int verbose)
{
  float e11, e12, e13, e21, e22, e23, e31, e32, e33;
  float lam2mu, lam, mu;
  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];

  float *matVx2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * FD_NDIM * FD_NDIM,
                                      0.0,
                                      "sv_eliso1st_curv_macdrp_vel_dxy2dz");

  float *matVy2Vz = (float *)fdlib_mem_calloc_1d_float(
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
  float *restrict lam3d = m3d + MD_EL_ISO_SEQ_LAMBDA * siz_volume;
  float *restrict  mu3d = m3d + MD_EL_ISO_SEQ_MU     * siz_volume;

  int k = nk2;

  for (size_t j = nj1; j <= nj2; j++)
  {
    for (size_t i = ni1; i <= ni2; i++)
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

      // first dim: irow; sec dim: jcol, as Fortran code
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

      size_t ij = (j * siz_line + i) * 9;

      // save into mat
      for(int irow = 0; irow < 3; irow++)
        for(int jcol = 0; jcol < 3; jcol++){
          matVx2Vz[ij + irow*3 + jcol] = AB[irow][jcol];
          matVy2Vz[ij + irow*3 + jcol] = AC[irow][jcol];
        }
    }
  }

  *p_matVx2Vz = matVx2Vz;
  *p_matVy2Vz = matVy2Vz;
}

void
sv_eliso1st_curv_macdrp_snap_buff(float *restrict var,
                                  size_t siz_line,
                                  size_t siz_slice,
                                  int starti,
                                  int counti,
                                  int increi,
                                  int startj,
                                  int countj,
                                  int increj,
                                  int startk,
                                  int countk,
                                  int increk,
                                  float *restrict buff)
{
  int iptr_snap=0;
  for (int n3=0; n3<countk; n3++)
  {
    int k = startk + n3 * increk;
    for (int n2=0; n2<countj; n2++)
    {
      int j = startj + n2 * increj;
      for (int n1=0; n1<counti; n1++)
      {
        int i = starti + n1 * increi;
        int iptr = i + j * siz_line + k * siz_slice;
        buff[iptr_snap] = var[iptr];
        iptr_snap++;
      }
    }
  }
}

void
sv_eliso1st_curv_macdrp_snap_buff_strain(float *restrict var,
                                  size_t siz_line,
                                  size_t siz_slice,
                                  int starti,
                                  int counti,
                                  int increi,
                                  int startj,
                                  int countj,
                                  int increj,
                                  int startk,
                                  int countk,
                                  int increk,
                                  float *restrict buff)
{
  int iptr_snap=0;
  for (int n3=0; n3<countk; n3++)
  {
    int k = startk + n3 * increk;
    for (int n2=0; n2<countj; n2++)
    {
      int j = startj + n2 * increj;
      for (int n1=0; n1<counti; n1++)
      {
        int i = starti + n1 * increi;
        int iptr = i + j * siz_line + k * siz_slice;
        buff[iptr_snap] = var[iptr];
        iptr_snap++;
      }
    }
  }
}

//void
//sv_eliso1st_curv_macdrp_snap_create_nc(char *ou_file,
//                                       int   snap_nt_total,
//                                       int   snap_ni,
//                                       int   snap_nj,
//                                       int   snap_nk,
//                                       int   snap_out_vel,
//                                       int   snap_out_stress,
//                                       int   snap_out_strain,
//                                       int  *ncid,
//                                       int  *varid_V,
//                                       int  *varid_T,
//                                       int  *varid_E)
//{
//    int dimid[4];
//    //sprintf(ou_file,"%s/%s_%s.nc", output_dir,snap_name[n],output_fname_part);
//
//    int snap_ni  = cur_snap_info[FD_NDIM  ];
//    int snap_nj  = cur_snap_info[FD_NDIM+1];
//    int snap_nk  = cur_snap_info[FD_NDIM+2];
//    int snap_it1 = cur_snap_info[FD_NDIM*3+0];
//    int snap_dit = cur_snap_info[FD_NDIM*3+1];
//    int snap_nt_total = (nt_total - snap_it1) / snap_it;
//    int snap_out_vel    = cur_snap_info[FD_NDIM*3+2+0];
//    int snap_out_stress = cur_snap_info[FD_NDIM*3+2+1];
//    int snap_out_strain = cur_snap_info[FD_NDIM*3+2+2];
//
//    //if (strcmp(strrchr(snap_fname[n],'.'),".nc")==0) 
//
//}
