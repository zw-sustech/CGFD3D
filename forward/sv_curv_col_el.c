/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"

#include "fd_t.h"
#include "mympi_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"

/*******************************************************************************
 * free surface
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_curv_col_el_rhs_timg_z2(
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

/*******************************************************************************
 * add source terms
 ******************************************************************************/

int
sv_curv_col_el_rhs_src(
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

  // get fi / mij
  float fx, fy, fz;
  float Mxx,Myy,Mzz,Mxz,Myz,Mxy;

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
      if (src->force_actived == 1) {
        fx  = src->Fx [iptr_cur_stage];
        fy  = src->Fy [iptr_cur_stage];
        fz  = src->Fz [iptr_cur_stage];
      }
      if (src->moment_actived == 1) {
        Mxx = src->Mxx[iptr_cur_stage];
        Myy = src->Myy[iptr_cur_stage];
        Mzz = src->Mzz[iptr_cur_stage];
        Mxz = src->Mxz[iptr_cur_stage];
        Myz = src->Myz[iptr_cur_stage];
        Mxy = src->Mxy[iptr_cur_stage];
      }
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];

        if (src->force_actived == 1) {
          float V = coef * slw3d[iptr] / jac3d[iptr];
          hVx[iptr] += fx * V;
          hVy[iptr] += fy * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / jac3d[iptr];
          hTxx[iptr] -= Mxx * rjac;
          hTyy[iptr] -= Myy * rjac;
          hTzz[iptr] -= Mzz * rjac;
          hTxz[iptr] -= Mxz * rjac;
          hTyz[iptr] -= Myz * rjac;
          hTxy[iptr] -= Mxy * rjac;
        }
      } // i_ext

    } // it
  } // is

  return ierr;
}

/*
 * add group of sources per time step, naming as inject to differentiate with previous
 * function, not strictly means injection method
 */

int
sv_curv_col_el_rhs_srcdd(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict jac3d, float *restrict slw3d,
    src_t *src, // short nation for reference member
    const int myid, const int verbose)
{
  int ierr = 0;

  // get info from src
  int it_here= src->dd_it_here;
  int istage = src->istage;

  if (verbose>999) {
    fprintf(stdout,"-- add srcdd: myid=%d,it_here=%d,istage=%d\n",myid,it_here,istage);
    fflush(stdout);
  }

  // return if pass time range: compare out of func for efficiency
  //if (it >= src->dd_max_nt) {
  //  return 0;
  //}

  // add injection src; is is a commont iterater var
  if (src->dd_is_add_at_point == 1)
  {
    // vi
    if (src->dd_vi_actived == 1) 
    {
      size_t iptr0_stf = it_here * src->max_stage * src->dd_total_number * CONST_NDIM
                                         + istage * src->dd_total_number * CONST_NDIM;
      float *dd_vi_pt = src->dd_vi + iptr0_stf;

      for (int is=0; is < src->dd_total_number; is++)
      {
        size_t iptr      = src->dd_indx[is];
        float V = slw3d[iptr] / jac3d[iptr];

        hVx[iptr] += dd_vi_pt[0] * V;
        hVy[iptr] += dd_vi_pt[1] * V;
        hVz[iptr] += dd_vi_pt[2] * V;

        // next source
        dd_vi_pt += CONST_NDIM;
      }
    }

    // mij
    if (src->dd_mij_actived == 1) 
    {
      size_t iptr0_stf = it_here * src->max_stage * src->dd_total_number * CONST_NDIM_2
                                         + istage * src->dd_total_number * CONST_NDIM_2;
      float *dd_mij_pt = src->dd_mij + iptr0_stf;
#ifdef DEBUG_SV_EQ1ST_CURV_COLGRD_ISO
      fprintf(stdout,"-- add srcdd mij: it_here=%d,istage=%d,max_stage=%d,dd_total_number=%d\n",
                          it_here,istage,src->max_stage,src->dd_total_number);
      fprintf(stdout,"-- add srcdd mij: iptr0_stf=%zu\n",
                          iptr0_stf);
      fflush(stdout);
#endif

      for (int is=0; is < src->dd_total_number; is++)
      {
        size_t iptr      = src->dd_indx[is];
        float rjac = 1.0 / jac3d[iptr];
#ifdef DEBUG_SV_EQ1ST_CURV_COLGRD_ISO
        fprintf(stdout,"-- add srcdd mij: is=%d,iptr=%zu,dd_mpi_pt=%p\n",
                            is,iptr,dd_mij_pt);
        fprintf(stdout,"-- add srcdd mij: mxx=%g,myy=%g,mzz=%g,myz=%g,mxz=%g,mxy=%g\n",
                  dd_mij_pt[0], dd_mij_pt[1], dd_mij_pt[2], dd_mij_pt[3], dd_mij_pt[4], dd_mij_pt[5]);
        fflush(stdout);
#endif

        hTxx[iptr] -= dd_mij_pt[0] * rjac;
        hTyy[iptr] -= dd_mij_pt[1] * rjac;
        hTzz[iptr] -= dd_mij_pt[2] * rjac;
        hTxz[iptr] -= dd_mij_pt[3] * rjac;
        hTyz[iptr] -= dd_mij_pt[4] * rjac;
        hTxy[iptr] -= dd_mij_pt[5] * rjac;

        dd_mij_pt += CONST_NDIM_2;
      }
    }
  }
  // add with spatial smooth
  else
  {
    fprintf(stderr,"ERROR: ddsrc with smoooth not implemented\n");
    fflush(stderr);
    MPI_Abort(MPI_COMM_WORLD, 1);

    //float wid_gauss = src->dd_smo_hlen / 2.0;
    //int   smo_heln  = src->dd_smo_hlen;

    //// mij
    //if (src->dd_mij_actived == 1) 
    //{
    //  size_t iptr0_stf = it_here * src->max_stage * src->dd_total_number * CONST_NDIM_2
    //                                     + istage * src->dd_total_number * CONST_NDIM_2;
    //  float *dd_mij_pt = src->dd_mij + iptr0_stf;

    //  for (int is=0; is < src->dd_total_number; is++)
    //  {
    //    size_t iptr      = src->dd_indx[is];
    //    float rjac = 1.0 / jac3d[iptr];

    //    int   si = src->dd_si[is];
    //    int   sj = src->dd_sj[is];
    //    int   sk = src->dd_sk[is];
    //    float si_inc = src->dd_si_inc[is];
    //    float sj_inc = src->dd_si_inc[is];
    //    float sk_inc = src->dd_si_inc[is];

    //    hTxx[iptr] -= dd_mij_pt[0] * rjac;
    //    hTyy[iptr] -= dd_mij_pt[1] * rjac;
    //    hTzz[iptr] -= dd_mij_pt[2] * rjac;
    //    hTxz[iptr] -= dd_mij_pt[3] * rjac;
    //    hTyz[iptr] -= dd_mij_pt[4] * rjac;
    //    hTxy[iptr] -= dd_mij_pt[5] * rjac;

    //    dd_mij_pt += CONST_NDIM_2;

    //    for (int ismo=-smo_heln; ismo<smo_heln; ismo++)
    //    {
    //      for (int jsmo=-smo_heln; jsmo<smo_heln; jsmo++)
    //      {
    //        for (int ksmo=-smo_heln; ksmo<smo_heln; ksmo++)
    //        {
    //          int i = si + ismo;
    //          int j = sj + jsmo;
    //          int k = sk + ksmo;
    //          if (    i<=ni2 && i>=ni1 
    //               && j<=nj2 && j>=nj1 
    //               && k<=nk2 && k>=nk1)
    //          {
    //          }
    //        }
    //      }
    //    }
    //  }
    //}
  }

  return ierr;
}

/*******************************************************************************
 * add attenuation
 ******************************************************************************/

int
sv_curv_col_el_graves_Qs(float *w, int ncmp, float dt, gd_t *gdinfo, md_t *md)
{
  int ierr = 0;

  float coef = - PI * md->visco_Qs_freq * dt;

  for (int icmp=0; icmp<ncmp; icmp++)
  {
    float *restrict var = w + icmp * gdinfo->siz_icmp;

    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      for (int j = gdinfo->nj1; j <= gdinfo->nj2; j++)
      {
        for (int i = gdinfo->ni1; i <= gdinfo->ni2; i++)
        {
          size_t iptr = i + j * gdinfo->siz_iy + k * gdinfo->siz_iz;

          float Qatt = expf( coef / md->Qs[iptr] );

          var[iptr] *= Qatt;
        }
      }
    }
  }

  return ierr;
}

