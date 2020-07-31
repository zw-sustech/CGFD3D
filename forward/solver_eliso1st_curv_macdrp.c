
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

#define CAL_FD_VAL(deriv, var, iptr, fd_length, fd_shift, fd_coef, n) \
   deriv = 0.0; \
   for (n=0; n<fd_length; n++) { \
       deriv += fd_coef[n] * var[iptr + fd_shift[n]]; \
   }

int elastic_iso_curv_macdrp_rkint(
    float *restrict w3d, // wavefield
    float *restrict g3d, // grid vars
    float *restrict m3d, // medium vars
    float *restrict aux, // aux vars, eg. PML
    // grid size
    size_t nx, size_t ny, size_t nz,
    // time
    float dt, int nt_total,
    // scheme
    int numPair, int numStage, float *rka, float *rkb,
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


  pos_pre = 0 * siz_vars;
  pos_rhs = 1 * siz_vars;
  pos_cur = 2 * siz_vars; // current value for rhs
  pos_end = 3 * siz_vars;

  while ( it < nt_total )
  {
      if (myid==0) fprintf(stdout,"-> it=%d, t=\%f\n", it, (it-0)*stept);

      for (iPair=0; iPair<numPair; iPair++)
      {
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
              ierr = elastic_iso_curv_pack_message(W_cur,buff,
                      fdy_length[1],fdy_length[2]);

              // exchange data along y
              MPI_Request r_reqs[4];

              MPI_Startall(2, reqs+0, MPI_STATUS_IGNORE);

              // wait for communication
              MPI_Waitall(2, reqs+0, MPI_STATUS_IGNORE);

              // receive message and uppack
              ierr = curv_macdrp_elastic_iso_unpack_message();

              // compute
              ierr = elastic_iso_curv_macdrp_cal_rhs(
                      w3d+pos_pre,
                      w3d+pos_rhs,
                      fd_opx_length,
                      fd_opx_shift,
                      fd_opx_coef,
                      fd_bdry_opx_nlay,

                      scheme_op_lenth
                      scheme_op_indx[iPair][iStage],
                      scheme_op_coef[iPair][iStage],
                      );

              // RK int, write here for hiding communication
              if (iStage<numStage-1) {
                  // wavefield
                  for (iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
                      w_cur[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
                  }
                  // aux
                  for (iptr=0; iptr<siz_aux_allvars; iptr++) {
                      a_cur[iptr] = a_pre[iptr] + coef_a * a_rhs[iptr];
                  }
                  // pack and isend
                  ierr = pack_message(w_cur,a_cur,sbuff);
                  MPI_Startall(sreqs);
              } else {
                  // wavefield
                  for (iptr=0; iptr<num_of_vars*siz_volume; iptr++) {
                      w_end[iptr] += coef_b * w_rhs[iptr];
                  }
                  // aux
                  for (iptr=0; iptr<siz_aux_allvars; iptr++) {
                      a_end[iptr] += coef_b * a_rhs[iptr];
                  }
                  // exp abs

                  // pack and isend
                  ierr = pack_message(w_end,a_end,sbuff);
                  MPI_Startall(sreqs);
              }

              // rk final value
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
          } // RK stages

          // swap w_pre and w_end, avoid copying
          w_tmp = w_pre; w_pre = w_end; w_end = w_tmp;

          // time step + 1
          it++;

          // QC
          if (par->qc_check_value_ornot==1 && id/par->qc_check_value_tinv==0) {
              ierr = check_value();
          }

          // save results
          ierr = io_seismo_keep();
          ierr = io_snapshot_save();

      } // forw/back pair
  } // time loop
}

//int elastic_iso_curv_macdrp_cal_rhs(float *restrict w_cur, float *restrict rhs, 
//        int *restrict fdx_length, size_t *restrict fdx_shift, float *restrict fdx_coef,
//        int *restrict fdx_bdry_nlay, int *restrict fdx_bdary_nlay_length, int *restrict fdx_bdary_nlay_pos,
//        int *restrict fdx_bdary_nlay_shift, float *restrict fdx_bdary_nlay_coef,
//        int *restrict fdy_length, size_t *restrict fdy_shift, float *restrict fdy_coef,
//        int *restrict fdy_bdry_nlay, int *restrict fdy_bdary_nlay_length, int *restrict fdy_bdary_nlay_pos,
//        int *restrict fdy_bdary_nlay_shift, float *restrict fdy_bdary_nlay_coef,
//        int *restrict fdz_length, size_t *restrict fdz_shift, float *restrict fdz_coef,
//        int *restrict fdz_bdry_nlay, int *restrict fdz_bdary_nlay_length, int *restrict fdz_bdary_nlay_pos,
//        int *restrict fdz_bdary_nlay_shift, float *restrict fdz_bdary_nlay_coef,
//        int myid, int *restrict myid3
//        )

int elastic_iso_curv_macdrp_cal_rhs(float *restrict w_cur, float *restrict rhs, 
        int *restrict fdx_length, // 0: total, 1: left, 2: right
        size_t *restrict fdx_shift, float *restrict fdx_coef,
        int *restrict fdy_length, size_t *restrict fdy_shift, float *restrict fdy_coef,
        int *restrict fdz_length, size_t *restrict fdz_shift, float *restrict fdz_coef,
        int *restrict fdz_bdry_nlay, int *restrict fdz_bdary_nlay_length, int *restrict fdz_bdary_nlay_pos,
        int *restrict fdz_bdary_nlay_shift, float *restrict fdz_bdary_nlay_coef,
        int myid, int *restrict myid3
        )
{
  int ierr = 0;

  size_t i,j,k;

  float *restrict Vx;
  float *restrict Vy;
  float *restrict Vz;
  float  fd_opx_coef[FD_MAX_OP_LEN];
  size_t fd_opx_possft[FD_MAX_OP_LEN];

  iptr = fdx_dhs_indx[fd_op_half_len-1];
  for (i=0; i<fd_op_total_len; i++) {
    fd_opx_coef[i] = fdx_opcoef[iptr+i];
    fd_opx_indx[i] = fdx_opindx[iptr+i];
  }

  Vx = w_cur + ELASTIC_ISO_CURV_MACDRP_SEQ_VX;
  Vy = w_cur + ELASTIC_ISO_CURV_MACDRP_SEQ_VY;
  Vz = w_cur + ELASTIC_ISO_CURV_MACDRP_SEQ_VZ;

  xi_x = g3d + GRID_CURV_SEQ_XIX;
  xi_y = g3d + GRID_CURV_SEQ_XIY;
  xi_z = g3d + GRID_CURV_SEQ_XIZ;

  // loop all points
  for (k=nk1; k<nk2; k++)
  {
      iptr_k = k * siz_slice;
      for (j=indx[2]; j<indx[3]; j++)
      {
          iptr_j = iptr_k + j * siz_line;

          iptr = iptr_j + (ni1-1);

          for (i=indx[0]; i<indx[1]; i++)
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

              // derivatives
              DxVx = 0.0;
              for (n=0; n<fd_opx_length; n++) {
                  DxVx += fd_opx_coef[n] * Vx[iptr + fd_opx_shift[n]];
              }

              DxVy = 0.0;
              for (n=0; n<fd_opx_length; n++) {
                  DxVy += fd_opx_coef[n] * Vy[iptr + fd_opx_shift[n]];
              }

              DxVz = 0.0;
              for (n=0; n<fd_opx_length; n++) {
                  DxVz += fd_opx_coef[n] * Vz[iptr + fd_opx_shift[n]];
              }

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

  // boundary x1
  switch (boundary_itype[0]) {
    case BOUNDARY_TYPE_FREE :
      //ierr = elastic_iso_free_x1();
      fprintf(stderr, "not implemented yet!\n"); fflush(stderr);
      break;

    case BOUNDARY_TYPE_CFSPML : 
      ierr = elastic_iso_cfspml_cal_rhs_x1();
      break;

    case BOUNDARY_TYPE_ABSEXP : 
      ierr = elastic_iso_absexp_apply_x1();
      break;

    default:
  }

  // boundary x2
  switch (boundary_itype[1]) {
    case BOUNDARY_TYPE_FREE :
      //ierr = elastic_iso_free_x2();
      fprintf(stderr, "not implemented yet!\n"); fflush(stderr);
      break;

    case BOUNDARY_TYPE_CFSPML : 
      ierr = elastic_iso_cfspml_cal_rhs_x2();
      break;

    default:
  }

  // y1, y2 ,z1, z2

  // boundary z2
  switch (boundary_itype[5]) {
    case BOUNDARY_TYPE_FREE :
      // tractiong
      ierr = curv_macdrp_elastic_iso_free_timg_z2();
      // velocity: vlow
      ierr = curv_macdrp_elastic_iso_free_vlow_z2();
      break;

    case BOUNDARY_TYPE_CFSPML : 
      ierr = elastic_iso_cfspml_cal_rhs_z2();
      break;

    default:
  }

  // source term
  switch (source_itype) {
    case SOURCE_TYPE_POINT :
      ierr = curv_macdrp_elastic_iso_src_point();
      break;

    case SOURCE_TYPE_GAUSS :
      ierr = curv_macdrp_elastic_iso_src_gauss();
      break;

    case SOURCE_TYPE_SINC :
      ierr = curv_macdrp_elastic_iso_src_sinc();
      break;
  }

  return ierr;
}

/*
int curv_macdrp_elastic_iso_rk_first()
{
    pos_var = 0;
    for (ivar=0; ivar<num_of_vars; ivar++) {
        // start index of ivar
        pos_var += size_volume;

        // start index of a point
        pos_val = pos_var;
        for (iptr=0; iptr<siz_volume; iptr++)
        {
            pos_val += 1;
            w_cur[pos_val] = w_pre[pos_val] + coef_a * w_rhs[pos_val];
            w_end[pos_val] = w_pre[pos_val] + coef_b * w_rhs[pos_val];
        }
    }
}

int curv_macdrp_elastic_iso_rk_inn()
{
    pos_var = 0;
    for (ivar=0; ivar<num_of_vars; ivar++) {
        // start index of ivar
        pos_var += size_volume;

        // start index of a point
        pos_val = pos_var;
        for (iptr=0; iptr<siz_volume; iptr++)
        {
            pos_val += 1;
            w_cur[pos_val] = w_pre[pos_val] + coef_a * w_rhs[pos_val];
            w_end[pos_val] += coef_b * w_rhs[pos_val];
        }
    }
}

int curv_macdrp_elastic_iso_rk_end()
{
    size_t pos_var, pos_val, iptr;
    int ivar;

    pos_var = 0;

    for (ivar=0; ivar<num_of_vars; ivar++)
    {
        // start index of ivar
        pos_var += size_volume;

        // start index of a point
        pos_val = pos_var;
        for (iptr=0; iptr<siz_volume; iptr++)
        {
            pos_val += 1;
            w_end[pos_val] += coef_b * w_rhs[pos_val];
        }
    }
}
*/
