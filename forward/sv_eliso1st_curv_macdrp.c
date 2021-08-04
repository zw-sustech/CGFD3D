/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and macdrp schem
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"
#include "sv_eliso1st_curv_macdrp.h"

//#define SV_ELISO1ST_CURV_MACDRP_DEBUG

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 *  simple MPI exchange without computing-communication overlapping
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_allstep(
  fd_t            *fd,
  gdinfo_t        *gdinfo,
  gdcurv_metric_t *metric,
  mdeliso_t      *mdeliso,
  src_t      *src,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  wfel1st_t  *wfel1st,
  mympi_t    *mympi,
  iorecv_t   *iorecv,
  ioline_t   *ioline,
  ioslice_t  *ioslice,
  iosnap_t   *iosnap,
  // time
  float dt, int nt_total, float t0,
  char *output_fname_part,
  char *output_dir,
  int qc_check_nan_num_of_step,
  const int output_all, // qc all var
  const int verbose)
{
  // retrieve from struct
  int num_rk_stages = fd->num_rk_stages;
  float *rk_a = fd->rk_a;
  float *rk_b = fd->rk_b;

  int num_of_pairs     = fd->num_of_pairs;
  int fdx_max_half_len = fd->fdx_max_half_len;
  int fdy_max_half_len = fd->fdy_max_half_len;
  int fdz_max_len      = fd->fdz_max_len;
  int num_of_fdz_op    = fd->num_of_fdz_op;

  // mpi
  int myid = mympi->myid;
  int *topoid = mympi->topoid;
  MPI_Comm comm = mympi->comm;
  float *restrict sbuff = mympi->sbuff;
  float *restrict rbuff = mympi->rbuff;
  MPI_Request *s_reqs = mympi->s_reqs;
  MPI_Request *r_reqs = mympi->r_reqs;

  // local allocated array
  char ou_file[CONST_MAX_STRLEN];

  // local pointer
  float *restrict w_cur;
  float *restrict w_pre;
  float *restrict w_rhs;
  float *restrict w_end;
  float *restrict w_tmp;

  int   ipair, istage;
  float t_cur;
  float t_end; // time after this loop for nc output

  // create slice nc output files
  if (myid==0 && verbose>0) fprintf(stdout,"prepare slice nc output ...\n"); 
  ioslice_nc_t ioslice_nc;
  io_slice_nc_create(ioslice, wfel1st->ncmp, wfel1st->cmp_name,
                     gdinfo->ni, gdinfo->nj, gdinfo->nk, topoid,
                     &ioslice_nc);

  // create snapshot nc output files
  if (myid==0 && verbose>0) fprintf(stdout,"prepare snap nc output ...\n"); 
  iosnap_nc_t  iosnap_nc;
  io_snap_nc_create(iosnap, &iosnap_nc, topoid);

  // only x/y mpi
  int num_of_r_reqs = 4;
  int num_of_s_reqs = 4;

  // get wavefield
  w_pre = wfel1st->v5d + wfel1st->siz_ilevel * 0; // previous level at n
  w_tmp = wfel1st->v5d + wfel1st->siz_ilevel * 1; // intermidate value
  w_rhs = wfel1st->v5d + wfel1st->siz_ilevel * 2; // for rhs
  w_end = wfel1st->v5d + wfel1st->siz_ilevel * 3; // end level at n+1

  // set pml for rk
  for (int idim=0; idim<CONST_NDIM; idim++) {
    for (int iside=0; iside<2; iside++) {
      if (bdrypml->is_at_sides[idim][iside]==1) {
        bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
        auxvar->pre = auxvar->var + auxvar->siz_ilevel * 0;
        auxvar->tmp = auxvar->var + auxvar->siz_ilevel * 1;
        auxvar->rhs = auxvar->var + auxvar->siz_ilevel * 2;
        auxvar->end = auxvar->var + auxvar->siz_ilevel * 3;
      }
    }
  }

  //--------------------------------------------------------
  // time loop
  //--------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"start time loop ...\n"); 

  for (int it=0; it<nt_total; it++)
  {
    t_cur  = it * dt + t0;
    t_end = t_cur +dt;

    if (myid==0 && verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // mod to get ipair
    ipair = it % num_of_pairs;
    if (myid==0 && verbose>10) fprintf(stdout, " --> ipair=%d\n",ipair);

    // loop RK stages for one step
    for (istage=0; istage<num_rk_stages; istage++)
    {
      if (myid==0 && verbose>10) fprintf(stdout, " --> istage=%d\n",istage);

      // use pointer to avoid 1 copy for previous level value
      if (istage==0) {
        w_cur = w_pre;
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            bdrypml->auxvar[idim][iside].cur = bdrypml->auxvar[idim][iside].pre;
          }
        }
      }
      else
      {
        w_cur = w_tmp;
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            bdrypml->auxvar[idim][iside].cur = bdrypml->auxvar[idim][iside].tmp;
          }
        }
      }

      // set src_t time
      src_set_time(src, it, istage);

      // compute rhs
      sv_eliso1st_curv_macdrp_onestage(
          w_cur,w_rhs,wfel1st,
          gdinfo, metric, mdeliso, bdryfree, bdrypml, src,
          fd->num_of_fdx_op, fd->pair_fdx_op[ipair][istage],
          fd->num_of_fdy_op, fd->pair_fdy_op[ipair][istage],
          fd->num_of_fdz_op, fd->pair_fdz_op[ipair][istage],
          fd->fdz_max_len,
          myid, verbose);

      // recv mesg
      MPI_Startall(num_of_r_reqs, r_reqs);

      // rk start
      if (istage==0)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr < wfel1st->siz_ilevel; iptr++) {
            w_tmp[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
        }
        // pack and isend
        blk_pack_mesg(w_tmp, sbuff, wfel1st->ncmp, gdinfo,
                         fdx_max_half_len, fdy_max_half_len);

        MPI_Startall(num_of_s_reqs, s_reqs);

        // pml_tmp
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            if (bdrypml->is_at_sides[idim][iside]==1) {
              bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
              for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                auxvar->tmp[iptr] = auxvar->pre[iptr] + coef_a * auxvar->rhs[iptr];
              }
            }
          }
        }

        // w_end
        for (size_t iptr=0; iptr < wfel1st->siz_ilevel; iptr++) {
            w_end[iptr] = w_pre[iptr] + coef_b * w_rhs[iptr];
        }
        // pml_end
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            if (bdrypml->is_at_sides[idim][iside]==1) {
              bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
              for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                auxvar->end[iptr] = auxvar->pre[iptr] + coef_b * auxvar->rhs[iptr];
              }
            }
          }
        }
      }
      else if (istage<num_rk_stages-1)
      {
        float coef_a = rk_a[istage] * dt;
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr < wfel1st->siz_ilevel; iptr++) {
            w_tmp[iptr] = w_pre[iptr] + coef_a * w_rhs[iptr];
        }
        // pack and isend
        blk_pack_mesg(w_tmp, sbuff, wfel1st->ncmp, gdinfo,
                         fdx_max_half_len, fdy_max_half_len);

        MPI_Startall(num_of_s_reqs, s_reqs);

        // pml_tmp
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            if (bdrypml->is_at_sides[idim][iside]==1) {
              bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
              for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                auxvar->tmp[iptr] = auxvar->pre[iptr] + coef_a * auxvar->rhs[iptr];
              }
            }
          }
        }

        // w_end
        for (size_t iptr=0; iptr < wfel1st->siz_ilevel; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }
        // pml_end
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            if (bdrypml->is_at_sides[idim][iside]==1) {
              bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
              for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                auxvar->end[iptr] += coef_b * auxvar->rhs[iptr];
              }
            }
          }
        }
      }
      else // last stage
      {
        float coef_b = rk_b[istage] * dt;

        // wavefield
        for (size_t iptr=0; iptr < wfel1st->siz_ilevel; iptr++) {
            w_end[iptr] += coef_b * w_rhs[iptr];
        }
        // pack and isend
        blk_pack_mesg(w_end, sbuff, wfel1st->ncmp, gdinfo,
                         fdx_max_half_len, fdy_max_half_len);

        MPI_Startall(num_of_s_reqs, s_reqs);

        // pml_end
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            if (bdrypml->is_at_sides[idim][iside]==1) {
              bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
              for (size_t iptr=0; iptr < auxvar->siz_ilevel; iptr++) {
                auxvar->end[iptr] += coef_b * auxvar->rhs[iptr];
              }
            }
          }
        }
      }

      MPI_Waitall(num_of_s_reqs, s_reqs, MPI_STATUS_IGNORE);
      MPI_Waitall(num_of_r_reqs, r_reqs, MPI_STATUS_IGNORE);

      if (istage != num_rk_stages-1) {
        blk_unpack_mesg(rbuff, w_tmp,  wfel1st->ncmp, gdinfo,
                         fdx_max_half_len, fdy_max_half_len);
     } else {
        blk_unpack_mesg(rbuff, w_end,  wfel1st->ncmp, gdinfo,
                         fdx_max_half_len, fdy_max_half_len);
     }
    } // RK stages

    //--------------------------------------------
    // QC
    //--------------------------------------------

    if (qc_check_nan_num_of_step >0  && (it % qc_check_nan_num_of_step) == 0) {
      if (myid==0 && verbose>10) fprintf(stdout,"-> check value nan\n");
        //wf_el_1st_check_value(w_end);
    }

    //--------------------------------------------
    // save results
    //--------------------------------------------

    //-- recv by interp
    io_recv_keep(iorecv, w_end, it, wfel1st->ncmp, wfel1st->siz_icmp);

    //-- line values
    io_line_keep(ioline, w_end, it, wfel1st->ncmp, wfel1st->siz_icmp);

    // write slice, use w_rhs as buff
    io_slice_nc_put(ioslice,&ioslice_nc,gdinfo,w_end,w_rhs,it,t_end,wfel1st->ncmp);

    // snapshot
    io_snap_nc_put(iosnap, &iosnap_nc, gdinfo, wfel1st, 
                   w_end, w_rhs, nt_total, it, t_end);

    // zero temp used w_rsh
    sv_eliso1st_curv_macdrp_zero_edge(gdinfo, wfel1st, w_rhs);

    // debug output
    if (output_all==1)
    {
        io_build_fname_time(output_dir,"w3d",".nc",topoid,it,ou_file);
        io_var3d_export_nc(ou_file,
                           w_end,
                           wfel1st->cmp_pos,
                           wfel1st->cmp_name,
                           wfel1st->ncmp,
                           gdinfo->index_name,
                           gdinfo->nx,
                           gdinfo->ny,
                           gdinfo->nz);
    }

    // swap w_pre and w_end, avoid copying
    w_cur = w_pre; w_pre = w_end; w_end = w_cur;

    for (int idim=0; idim<CONST_NDIM; idim++) {
      for (int iside=0; iside<2; iside++) {
        bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
        auxvar->cur = auxvar->pre;
        auxvar->pre = auxvar->end;
        auxvar->end = auxvar->cur;
      }
    }

  } // time loop

  // postproc

  // close nc
  io_slice_nc_close(&ioslice_nc);
  io_snap_nc_close(&iosnap_nc);

}

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eliso1st_curv_macdrp_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wfel1st_t  *wfel1st,
  gdinfo_t   *gdinfo,
  gdcurv_metric_t  *metric,
  mdeliso_t *mdeliso,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  // include different order/stentil
  int num_of_fdx_op, fd_op_t *fdx_op,
  int num_of_fdy_op, fd_op_t *fdy_op,
  int num_of_fdz_op, fd_op_t *fdz_op,
  int fdz_max_len, 
  const int myid, const int verbose)
{
  // local pointer get each vars
  float *restrict Vx    = w_cur + wfel1st->Vx_pos ;
  float *restrict Vy    = w_cur + wfel1st->Vy_pos ;
  float *restrict Vz    = w_cur + wfel1st->Vz_pos ;
  float *restrict Txx   = w_cur + wfel1st->Txx_pos;
  float *restrict Tyy   = w_cur + wfel1st->Tyy_pos;
  float *restrict Tzz   = w_cur + wfel1st->Tzz_pos;
  float *restrict Txz   = w_cur + wfel1st->Tyz_pos;
  float *restrict Tyz   = w_cur + wfel1st->Txz_pos;
  float *restrict Txy   = w_cur + wfel1st->Txy_pos;
  float *restrict hVx   = rhs   + wfel1st->Vx_pos ; 
  float *restrict hVy   = rhs   + wfel1st->Vy_pos ; 
  float *restrict hVz   = rhs   + wfel1st->Vz_pos ; 
  float *restrict hTxx  = rhs   + wfel1st->Txx_pos; 
  float *restrict hTyy  = rhs   + wfel1st->Tyy_pos; 
  float *restrict hTzz  = rhs   + wfel1st->Tzz_pos; 
  float *restrict hTxz  = rhs   + wfel1st->Tyz_pos; 
  float *restrict hTyz  = rhs   + wfel1st->Txz_pos; 
  float *restrict hTxy  = rhs   + wfel1st->Txy_pos; 

  float *restrict xi_x  = metric->xi_x;
  float *restrict xi_y  = metric->xi_y;
  float *restrict xi_z  = metric->xi_z;
  float *restrict et_x  = metric->eta_x;
  float *restrict et_y  = metric->eta_x;
  float *restrict et_z  = metric->eta_x;
  float *restrict zt_x  = metric->zeta_x;
  float *restrict zt_y  = metric->zeta_x;
  float *restrict zt_z  = metric->zeta_x;
  float *restrict jac3d = metric->jac;

  float *restrict lam3d = mdeliso->lambda;
  float *restrict  mu3d = mdeliso->mu;
  float *restrict slw3d = mdeliso->rho;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;

  int ni  = gdinfo->ni;
  int nj  = gdinfo->nj;
  int nk  = gdinfo->nk;
  int nx  = gdinfo->nx;
  int ny  = gdinfo->ny;
  int nz  = gdinfo->nz;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;
  size_t siz_volume = gdinfo->siz_volume;

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;

  // local fd op
  int              fdx_inn_len;
  int    *restrict fdx_inn_indx;
  float *restrict fdx_inn_coef;
  int              fdy_inn_len;
  int    *restrict fdy_inn_indx;
  float  *restrict fdy_inn_coef;
  int              fdz_inn_len;
  int    *restrict fdz_inn_indx;
  float  *restrict fdz_inn_coef;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op
  fdx_inn_len  = fdx_op[num_of_fdx_op-1].total_len;
  fdx_inn_indx = fdx_op[num_of_fdx_op-1].indx;
  fdx_inn_coef = fdx_op[num_of_fdx_op-1].coef;

  fdy_inn_len  = fdy_op[num_of_fdy_op-1].total_len;
  fdy_inn_indx = fdy_op[num_of_fdy_op-1].indx;
  fdy_inn_coef = fdy_op[num_of_fdy_op-1].coef;

  fdz_inn_len  = fdz_op[num_of_fdz_op-1].total_len;
  fdz_inn_indx = fdz_op[num_of_fdz_op-1].indx;
  fdz_inn_coef = fdz_op[num_of_fdz_op-1].coef;

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
  if (bdryfree->is_at_sides[2][1] == 1)
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
                                        lam3d, mu3d, slw3d,
                                        matVx2Vz,matVy2Vz,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        num_of_fdz_op,fdz_op,fdz_max_len,
                                        myid, verbose);
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eliso1st_curv_macdrp_rhs_cfspml(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                       xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                       lam3d, mu3d, slw3d,
                                       nk2, siz_line,siz_slice,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       bdrypml, bdryfree,
                                       myid, verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_eliso1st_curv_macdrp_rhs_src(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    jac3d, slw3d, 
                                    src,
                                    myid, verbose);
  }
  // end func
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
  int    lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int    lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

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
    int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
    const int myid, const int verbose)
{
  // use local stack array for speedup
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int    lfdy_shift[fdy_len];

  // allocate max_len because fdz may have different lens
  float  lfdz_coef [fdz_max_len];
  int    lfdz_shift[fdz_max_len];

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
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    k = nk2 - n;

    // get pos and len for this point
    int  lfdz_len  = fdz_op[n].total_len;
    // point to indx/coef for this point
    int   *p_fdz_indx  = fdz_op[n].indx;
    float *p_fdz_coef  = fdz_op[n].coef;
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
    int nk2, size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    bdrypml_t *bdrypml, bdryfree_t *bdryfree,
    const int myid, const int verbose)
{

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;

  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef [fdx_len];
  int    lfdx_shift[fdx_len];
  float  lfdy_coef [fdy_len];
  int    lfdy_shift[fdy_len];
  float  lfdz_coef [fdz_len];
  int    lfdz_shift[fdz_len];

  // val on point
  float DxTxx,DxTyy,DxTzz,DxTxy,DxTxz,DxTyz,DxVx,DxVy,DxVz;
  float DyTxx,DyTyy,DyTzz,DyTxy,DyTxz,DyTyz,DyVx,DyVy,DyVz;
  float DzTxx,DzTyy,DzTzz,DzTxy,DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,lam2mu,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
  float hVx_rhs,hVy_rhs,hVz_rhs;
  float hTxx_rhs,hTyy_rhs,hTzz_rhs,hTxz_rhs,hTyz_rhs,hTxy_rhs;
  // for free surface
  float Dx_DzVx,Dy_DzVx,Dx_DzVy,Dy_DzVy,Dx_DzVz,Dy_DzVz;

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

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      // skip to next face if not cfspml
      if (bdrypml->is_at_sides[idim][iside] == 0) continue;

      // get index into local var
      int abs_ni1 = bdrypml->ni1[idim][iside];
      int abs_ni2 = bdrypml->ni2[idim][iside];
      int abs_nj1 = bdrypml->nj1[idim][iside];
      int abs_nj2 = bdrypml->nj2[idim][iside];
      int abs_nk1 = bdrypml->nk1[idim][iside];
      int abs_nk2 = bdrypml->nk2[idim][iside];

#ifdef SV_ELISO1ST_CURV_MACDRP_DEBUG
    //fprintf(stdout," iface=%d,ni1=%d,ni2=%d,nj1=%d,nj2=%d,nk1=%d,nk2=%d\n",
    //        iface,abs_ni1,abs_ni2,abs_nj1,abs_nj2,abs_nk1,abs_nk2);
    //fflush(stdout);
#endif

      // get coef for this face
      float *restrict ptr_coef_A = bdrypml->A[idim][iside];
      float *restrict ptr_coef_B = bdrypml->B[idim][iside];
      float *restrict ptr_coef_D = bdrypml->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);

      // get pml vars
      float *restrict abs_vars_cur = auxvar->cur;
      float *restrict abs_vars_rhs = auxvar->rhs;

      float *restrict pml_Vx   = abs_vars_cur + auxvar->Vx_pos;
      float *restrict pml_Vy   = abs_vars_cur + auxvar->Vy_pos;
      float *restrict pml_Vz   = abs_vars_cur + auxvar->Vz_pos;
      float *restrict pml_Txx  = abs_vars_cur + auxvar->Txx_pos;
      float *restrict pml_Tyy  = abs_vars_cur + auxvar->Tyy_pos;
      float *restrict pml_Tzz  = abs_vars_cur + auxvar->Tzz_pos;
      float *restrict pml_Txz  = abs_vars_cur + auxvar->Txz_pos;
      float *restrict pml_Tyz  = abs_vars_cur + auxvar->Tyz_pos;
      float *restrict pml_Txy  = abs_vars_cur + auxvar->Txy_pos;
      float *restrict pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
      float *restrict pml_hVy  = abs_vars_rhs + auxvar->Vy_pos;
      float *restrict pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
      float *restrict pml_hTxx = abs_vars_rhs + auxvar->Txx_pos;
      float *restrict pml_hTyy = abs_vars_rhs + auxvar->Tyy_pos;
      float *restrict pml_hTzz = abs_vars_rhs + auxvar->Tzz_pos;
      float *restrict pml_hTxz = abs_vars_rhs + auxvar->Txz_pos;
      float *restrict pml_hTyz = abs_vars_rhs + auxvar->Tyz_pos;
      float *restrict pml_hTxy = abs_vars_rhs + auxvar->Txy_pos;

      // for each dim
      if (idim == 0 ) // x direction
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

              // add contributions from free surface condition
              //  not consider timg because conflict with main cfspml,
              //     need to revise in the future if required
              if (bdryfree->is_at_sides[idim][iside]==1 && k==nk2)
              {
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
              }

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else if (idim == 1) // y direction
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

              // add contributions from free surface condition
              if (bdryfree->is_at_sides[idim][iside]==1 && k==nk2)
              {
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

                hTxx_rhs =    lam2mu * (             ztx*Dy_DzVx)
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
              }

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
      } // if which dim
    } // iside
  } // idim

  return;
}

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

int
sv_eliso1st_curv_macdrp_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    src_t *src, // short nation for reference member
    const int myid, const int verbose)
{
  int ierr = 0;

  // local var
  int si,sj,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  int it     = src->it;
  int istage = src->istage;

  // add src; is is a commont iterater var
  for (int is=0; is < src->total_number; is++)
  {
    int   it_start = src->it_begin[is];
    int   it_end   = src->it_end  [is];

    if (it >= it_start && it <= it_end)
    {
      int   *ptr_ext_indx = src->ext_indx + is * max_ext;
      float *ptr_ext_coef = src->ext_coef + is * max_ext;
      int it_to_it_start = it - it_start;
      int iptr_cur_stage =   is * src->max_nt * src->max_stage // skip other src
                           + it_to_it_start * src->max_stage // skip other time step
                           + istage;
      float fx  = src->Fx [iptr_cur_stage];
      float fy  = src->Fy [iptr_cur_stage];
      float fz  = src->Fz [iptr_cur_stage];
      float Mxx = src->Mxx[iptr_cur_stage];
      float Myy = src->Myy[iptr_cur_stage];
      float Mzz = src->Mzz[iptr_cur_stage];
      float Mxz = src->Mxz[iptr_cur_stage];
      float Myz = src->Myz[iptr_cur_stage];
      float Mxy = src->Mxy[iptr_cur_stage];
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];

        float V = coef * slw3d[iptr] / jac3d[iptr];
        hVx[iptr] += fx * V;
        hVy[iptr] += fy * V;
        hVz[iptr] += fz * V;

        float rjac = coef / jac3d[iptr];
        hTxx[iptr] -= Mxx * rjac;
        hTyy[iptr] -= Myy * rjac;
        hTzz[iptr] -= Mzz * rjac;
        hTxz[iptr] -= Mxz * rjac;
        hTyz[iptr] -= Myz * rjac;
        hTxy[iptr] -= Mxy * rjac;
      } // i_ext

    } // it
  } // is

  return ierr;
}

/*******************************************************************************
 * related functions
 ******************************************************************************/

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

int
sv_eliso1st_curv_macdrp_zero_edge(gdinfo_t *gdinfo, wfel1st_t *wfel1st,
                                  float *restrict w4d)
{
  int ierr = 0;

  for (int icmp=0; icmp < wfel1st->ncmp; icmp++)
  {
    float *restrict var = w4d + wfel1st->cmp_pos[icmp];

    // z1
    for (int k=0; k < gdinfo->nk1; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j=0; j < gdinfo->ny; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // z2
    for (int k=gdinfo->nk2+1; k < gdinfo->nz; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j=0; j < gdinfo->ny; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // y1
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j=0; j < gdinfo->nj1; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // y2
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j = gdinfo->nj2+1; j < gdinfo->ny; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // x1
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j = gdinfo->nj1; j <= gdinfo->nj2; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->ni1; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    } 

    // x2
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j = gdinfo->nj1; j <= gdinfo->nj2; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i = gdinfo->ni2+1; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    } 

  } // icmp

  return ierr;
}

// code spool
/*
        for (int idim=0; idim<CONST_NDIM; idim++) {
          for (int iside=0; iside<2; iside++) {
            if (bdrypml->is_at_sides[idim][iside]==1) {
              bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);
              auxvar->cur = auxvar->tmp;
            }
          }
        }
*/
