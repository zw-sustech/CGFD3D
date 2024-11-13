/*******************************************************************************
 * solver of isotropic visco-elastic 1st-order eqn using curv grid and collocated scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "sv_curv_col_el.h"
#include "sv_curv_col_vis_iso.h"
#include "sv_curv_col_el_iso.h"

//#define DEBUG_SV_EQ1ST_CURV_COLGRD_ISO

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_curv_col_vis_iso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gd_t   *gdinfo,
  gdcurv_metric_t  *metric,
  md_t *md,
  bdry_t  *bdry,
  src_t *src,
  // include different order/stentil
  int num_of_fdx_op, fd_op_t *fdx_op,
  int num_of_fdy_op, fd_op_t *fdy_op,
  int num_of_fdz_op, fd_op_t *fdz_op,
  int fdz_max_len, 
  const int myid, const int verbose)
{
    int nmaxwell = md->nmaxwell;

  // local pointer get each vars
  float *restrict Vx    = w_cur + wav->Vx_pos ;
  float *restrict Vy    = w_cur + wav->Vy_pos ;
  float *restrict Vz    = w_cur + wav->Vz_pos ;
  float *restrict Txx   = w_cur + wav->Txx_pos;
  float *restrict Tyy   = w_cur + wav->Tyy_pos;
  float *restrict Tzz   = w_cur + wav->Tzz_pos;
  float *restrict Txz   = w_cur + wav->Txz_pos;
  float *restrict Tyz   = w_cur + wav->Tyz_pos;
  float *restrict Txy   = w_cur + wav->Txy_pos;
  float *restrict hVx   = rhs   + wav->Vx_pos ; 
  float *restrict hVy   = rhs   + wav->Vy_pos ; 
  float *restrict hVz   = rhs   + wav->Vz_pos ; 
  float *restrict hTxx  = rhs   + wav->Txx_pos; 
  float *restrict hTyy  = rhs   + wav->Tyy_pos; 
  float *restrict hTzz  = rhs   + wav->Tzz_pos; 
  float *restrict hTxz  = rhs   + wav->Txz_pos; 
  float *restrict hTyz  = rhs   + wav->Tyz_pos; 
  float *restrict hTxy  = rhs   + wav->Txy_pos; 


  float **Jxx = (float **) malloc(nmaxwell * sizeof(float *));
  float **Jyy = (float **) malloc(nmaxwell * sizeof(float *));
  float **Jzz = (float **) malloc(nmaxwell * sizeof(float *));
  float **Jxy = (float **) malloc(nmaxwell * sizeof(float *));
  float **Jxz = (float **) malloc(nmaxwell * sizeof(float *));
  float **Jyz = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJxx = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJyy = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJzz = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJxy = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJxz = (float **) malloc(nmaxwell * sizeof(float *));
  float **hJyz = (float **) malloc(nmaxwell * sizeof(float *));

  for(int n=0; n < nmaxwell; n++)
  {
    Jxx[n]  = w_cur + wav->Jxx_pos[n];
    Jyy[n]  = w_cur + wav->Jyy_pos[n];
    Jzz[n]  = w_cur + wav->Jzz_pos[n];
    Jxy[n]  = w_cur + wav->Jxy_pos[n];
    Jyz[n]  = w_cur + wav->Jyz_pos[n];
    Jxz[n]  = w_cur + wav->Jxz_pos[n];
    hJxx[n] = rhs   + wav->Jxx_pos[n];
    hJyy[n] = rhs   + wav->Jyy_pos[n];
    hJzz[n] = rhs   + wav->Jzz_pos[n];
    hJxy[n] = rhs   + wav->Jxy_pos[n];
    hJyz[n] = rhs   + wav->Jyz_pos[n];
    hJxz[n] = rhs   + wav->Jxz_pos[n];
  }

  float *restrict xi_x  = metric->xi_x;
  float *restrict xi_y  = metric->xi_y;
  float *restrict xi_z  = metric->xi_z;
  float *restrict et_x  = metric->eta_x;
  float *restrict et_y  = metric->eta_y;
  float *restrict et_z  = metric->eta_z;
  float *restrict zt_x  = metric->zeta_x;
  float *restrict zt_y  = metric->zeta_y;
  float *restrict zt_z  = metric->zeta_z;
  float *restrict jac3d = metric->jac;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;
  float *restrict slw3d = md->rho;
  float *restrict wl    = md->wl;
  float **restrict Ylam = md->Ylam;
  float **restrict Ymu  = md->Ymu;

  float *restrict TxSrc = src->TxSrc;
  float *restrict TySrc = src->TySrc;
  float *restrict TzSrc = src->TzSrc;
  float *restrict VxSrc = src->VxSrc;
  float *restrict VySrc = src->VySrc;
  float *restrict VzSrc = src->VzSrc;

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

  float *matVx2Vz = bdry->matVx2Vz2;
  float *matVy2Vz = bdry->matVy2Vz2;
  float *matF2Vz  = bdry->matF2Vz2;

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
  sv_curv_col_el_iso_rhs_inner(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
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
  if (bdry->is_sides_free[2][1] == 1)
  {
    // tractiong
    sv_curv_col_el_rhs_timg_z2(Txx,Tyy,Tzz,Txz,Tyz,Txy,hVx,hVy,hVz,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        jac3d, slw3d,
                                        TxSrc,TySrc,TzSrc,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        myid, verbose);

    sv_curv_col_el_iso_rhs_vlow_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                      xi_x, xi_y, xi_z, et_x, et_y,et_z, zt_x, zt_y, zt_z,
                                      lam3d, mu3d, slw3d,
                                      matVx2Vz,matVy2Vz,matF2Vz,
                                      VxSrc,VySrc,VzSrc,
                                      ni1,ni2,nj1,nj2,nk1,nk2,
                                      siz_line,siz_slice,fdx_inn_len, fdx_inn_indx,
                                      fdx_inn_coef,fdy_inn_len, fdy_inn_indx,
                                      fdy_inn_coef,num_of_fdz_op,fdz_op,
                                      fdz_max_len,myid, verbose);
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_curv_col_el_iso_rhs_cfspml(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                       xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                       lam3d, mu3d, slw3d,
                                       nk2, siz_line,siz_slice,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       bdry, 
                                       myid, verbose);
  }

  sv_curv_col_vis_iso_atten(hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                            Jxx,Jyy,Jzz,Jxy,Jxz,Jyz,
                            hJxx,hJyy,hJzz,hJxy,hJxz,hJyz,
                            lam3d, mu3d, slw3d,
                            wl, Ylam, Ymu,
                            ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                            nmaxwell,
                            myid, verbose);

  // add source term
  if (src->total_number > 0)
  {
    sv_curv_col_el_rhs_src(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    jac3d, slw3d, 
                                    gdinfo,src,
                                    myid, verbose);
  }

  // not do if pass time range
  if (src->dd_is_valid ==1)
  {
    sv_curv_col_el_rhs_srcdd(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    jac3d, slw3d, 
                                    src,
                                    myid, verbose);
  }
  // end func
}


/*******************************************************************************
 * add the attenuation term
******************************************************************************/
void
sv_curv_col_vis_iso_atten(
    float *restrict hTxx,  float *restrict hTyy,  float *restrict hTzz,
    float *restrict hTxz,  float *restrict hTyz,  float *restrict hTxy,
    float **restrict Jxx,  float **restrict Jyy,  float **restrict Jzz,
    float **restrict Jxy,  float **restrict Jxz,  float **restrict Jyz,
    float **restrict hJxx, float **restrict hJyy, float **restrict hJzz,
    float **restrict hJxy, float **restrict hJxz, float **restrict hJyz,
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    float *restrict wl,    float **restrict Ylam, float **restrict Ymu,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int nmaxwell,
    const int myid, const int verbose)
{
  float lam,mu;
  float mem_Txx,mem_Tyy,mem_Tzz,mem_Txy,mem_Txz,mem_Tyz;
  float sum_Jxyz,sum_Jxx,sum_Jyy,sum_Jzz,sum_Jxy,sum_Jxz,sum_Jyz;
  float EVxx,EVyy,EVzz,EVxy,EVxz,EVyz;
  float sum_hxyz;

  // loop all points
  for (size_t k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];

        sum_hxyz = (hTxx[iptr]+hTyy[iptr]+hTzz[iptr])/(3*lam+2*mu);

        EVxx = ((2.0*hTxx[iptr]-hTyy[iptr]-hTzz[iptr])/(2*mu) + sum_hxyz)/3;

        EVyy = ((2.0*hTyy[iptr]-hTxx[iptr]-hTzz[iptr])/(2*mu) + sum_hxyz)/3;

        EVzz = ((2.0*hTzz[iptr]-hTxx[iptr]-hTyy[iptr])/(2*mu) + sum_hxyz)/3;

        EVxy = hTxy[iptr]/mu*0.5;
        EVxz = hTxz[iptr]/mu*0.5;
        EVyz = hTyz[iptr]/mu*0.5;
        
        for(int n=0; n<nmaxwell; n++)
        {
          hJxx[n][iptr] = wl[n] * (EVxx - Jxx[n][iptr]);
          hJyy[n][iptr] = wl[n] * (EVyy - Jyy[n][iptr]);
          hJzz[n][iptr] = wl[n] * (EVzz - Jzz[n][iptr]);
          hJxy[n][iptr] = wl[n] * (EVxy - Jxy[n][iptr]);
          hJxz[n][iptr] = wl[n] * (EVxz - Jxz[n][iptr]);
          hJyz[n][iptr] = wl[n] * (EVyz - Jyz[n][iptr]);
        }

        // sum of memory variable for attenuation
        sum_Jxyz = 0.0;
        sum_Jxx = 0.0;
        sum_Jyy = 0.0;
        sum_Jzz = 0.0;
        sum_Jxy = 0.0;
        sum_Jxz = 0.0;
        sum_Jyz = 0.0;
    
        for(int n=0; n<nmaxwell; n++)
        {
          sum_Jxyz += Ylam[n][iptr] * (Jxx[n][iptr]+Jyy[n][iptr]+Jzz[n][iptr]);
          sum_Jxx  += Ymu[n][iptr]  * Jxx[n][iptr];
          sum_Jyy  += Ymu[n][iptr]  * Jyy[n][iptr];
          sum_Jzz  += Ymu[n][iptr]  * Jzz[n][iptr];
          sum_Jxy  += Ymu[n][iptr]  * Jxy[n][iptr];
          sum_Jxz  += Ymu[n][iptr]  * Jxz[n][iptr];
          sum_Jyz  += Ymu[n][iptr]  * Jyz[n][iptr];
        }
    
        mem_Txx = lam*sum_Jxyz + 2.0*mu*sum_Jxx;
        mem_Tyy = lam*sum_Jxyz + 2.0*mu*sum_Jyy;
        mem_Tzz = lam*sum_Jxyz + 2.0*mu*sum_Jzz;
        mem_Txy = 2.0*mu*sum_Jxy;
        mem_Txz = 2.0*mu*sum_Jxz;
        mem_Tyz = 2.0*mu*sum_Jyz;

        hTxx[iptr] -= mem_Txx;
        hTyy[iptr] -= mem_Tyy;
        hTzz[iptr] -= mem_Tzz;
        hTxy[iptr] -= mem_Txy;
        hTxz[iptr] -= mem_Txz;
        hTyz[iptr] -= mem_Tyz;

        iptr += 1;
      }
    }
  }
}

/*******************************************************************************
 * free surface coef
 ******************************************************************************/

int
sv_curv_col_vis_iso_dvh2dvz(gd_t            *gdinfo,
                            gdcurv_metric_t *metric,
                            md_t            *md,
                            bdry_t          *bdryfree,
                            int fd_len,
                            int *restrict fd_indx,
                            float *restrict fd_coef,
                            const int verbose)
{
  int ierr = 0;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nx  = gdinfo->nx;
  int ny  = gdinfo->ny;
  int nz  = gdinfo->nz;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;
  size_t siz_volume = gdinfo->siz_icmp;

  // point to each var
  float *restrict x3d  = gdinfo->x3d;
  float *restrict y3d  = gdinfo->y3d;
  float *restrict z3d  = gdinfo->z3d;
  float *restrict xi_x = metric->xi_x;
  float *restrict xi_y = metric->xi_y;
  float *restrict xi_z = metric->xi_z;
  float *restrict et_x = metric->eta_x;
  float *restrict et_y = metric->eta_y;
  float *restrict et_z = metric->eta_z;
  float *restrict zt_x = metric->zeta_x;
  float *restrict zt_y = metric->zeta_y;
  float *restrict zt_z = metric->zeta_z;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;
  float *matD = bdryfree->matD;
  
  float A[3][3], B[3][3], C[3][3], D[3][3];
  float AB[3][3], AC[3][3];

  float xix, xiy, xiz, etx, ety, etz, ztx, zty, ztz;
  float lam2mu, lam, mu;

  float xet,yet,zet;
  float e_n,e_m,e_nm;
  int n_fd;

  // use local stack array for speedup
  float  lfd_coef [fd_len];
  int    lfdy_shift[fd_len];
  // put fd op into local array
  for (int k=0; k < fd_len; k++) {
    lfd_coef [k] = fd_coef[k];
    lfdy_shift[k] = fd_indx[k] * siz_line ;
  }

  int k = nk2;

  for (size_t j = nj1; j <= nj2; j++)
  {
    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + j * siz_line + k * siz_slice;

      xix = xi_x[iptr];
      xiy = xi_y[iptr];
      xiz = xi_z[iptr];
      etx = et_x[iptr];
      ety = et_y[iptr];
      etz = et_z[iptr];
      ztx = zt_x[iptr];
      zty = zt_y[iptr];
      ztz = zt_z[iptr];

      xet = 0.0; yet = 0.0; zet = 0.0;

      M_FD_SHIFT(xet, x3d, iptr, fd_len, lfdy_shift, lfd_coef, n_fd);
      M_FD_SHIFT(yet, y3d, iptr, fd_len, lfdy_shift, lfd_coef, n_fd);
      M_FD_SHIFT(zet, z3d, iptr, fd_len, lfdy_shift, lfd_coef, n_fd);

      e_n = 1.0/sqrt(xet*xet+yet*yet+zet*zet);
      e_m = 1.0/sqrt(ztx*ztx+zty*zty+ztz*ztz);
      e_nm = e_n*e_m;

      lam    = lam3d[iptr];
      mu     =  mu3d[iptr];
      lam2mu = lam + 2.0f * mu;

      // first dim: irow; sec dim: jcol, as Fortran code
      A[0][0] = lam2mu*ztx*ztx + mu*(zty*zty+ztz*ztz);
      A[0][1] = lam*ztx*zty + mu*zty*ztx;
      A[0][2] = lam*ztx*ztz + mu*ztz*ztx;
      A[1][0] = lam*zty*ztx + mu*ztx*zty;
      A[1][1] = lam2mu*zty*zty + mu*(ztx*ztx+ztz*ztz);
      A[1][2] = lam*zty*ztz + mu*ztz*zty;
      A[2][0] = lam*ztz*ztx + mu*ztx*ztz;
      A[2][1] = lam*ztz*zty + mu*zty*ztz;
      A[2][2] = lam2mu*ztz*ztz + mu*(ztx*ztx+zty*zty);
      fdlib_math_invert3x3(A);

      B[0][0] = -lam2mu*ztx*xix - mu*(zty*xiy+ztz*xiz);
      B[0][1] = -lam*ztx*xiy - mu*zty*xix;
      B[0][2] = -lam*ztx*xiz - mu*ztz*xix;
      B[1][0] = -lam*zty*xix - mu*ztx*xiy;
      B[1][1] = -lam2mu*zty*xiy - mu*(ztx*xix+ztz*xiz);
      B[1][2] = -lam*zty*xiz - mu*ztz*xiy;
      B[2][0] = -lam*ztz*xix - mu*ztx*xiz;
      B[2][1] = -lam*ztz*xiy - mu*zty*xiz;
      B[2][2] = -lam2mu*ztz*xiz - mu*(ztx*xix+zty*xiy);

      C[0][0] = -lam2mu*ztx*etx - mu*(zty*ety+ztz*etz);
      C[0][1] = -lam*ztx*ety - mu*zty*etx;
      C[0][2] = -lam*ztx*etz - mu*ztz*etx;
      C[1][0] = -lam*zty*etx - mu*ztx*ety;
      C[1][1] = -lam2mu*zty*ety - mu*(ztx*etx+ztz*etz);
      C[1][2] = -lam*zty*etz - mu*ztz*ety;
      C[2][0] = -lam*ztz*etx - mu*ztx*etz;
      C[2][1] = -lam*ztz*ety - mu*zty*etz;
      C[2][2] = -lam2mu*ztz*etz - mu*(ztx*etx+zty*ety);

      fdlib_math_matmul3x3(A, B, AB);
      fdlib_math_matmul3x3(A, C, AC);

      D[0][0] = (yet*ztz-zet*zty)*e_nm;
      D[0][1] = (zet*ztx-xet*ztz)*e_nm;
      D[0][2] = (xet*zty-yet*ztx)*e_nm;
      D[1][0] = xet*e_n;
      D[1][1] = yet*e_n;
      D[1][2] = zet*e_n;
      D[2][0] = ztx*e_m;
      D[2][1] = zty*e_m;
      D[2][2] = ztz*e_m;
      size_t ij = (j * siz_line + i) * 9;

      // save into mat
      for(int irow = 0; irow < 3; irow++)
        for(int jcol = 0; jcol < 3; jcol++){
          matVx2Vz[ij + irow*3 + jcol] = AB[irow][jcol];
          matVy2Vz[ij + irow*3 + jcol] = AC[irow][jcol];
          matD[ij + irow*3 + jcol] = D[irow][jcol]; 
        }
    }
  }

  return ierr;
}

int
sv_curv_col_vis_iso_free(float *restrict w_end,
                         wav_t  *wav,
                         gd_t   *gdinfo,
                         gdcurv_metric_t  *metric,
                         md_t *md,
                         bdry_t      *bdryfree,const int myid, 
                         const int verbose)
{
  int ierr = 0;

  float *restrict Txx   = w_end + wav->Txx_pos;
  float *restrict Tyy   = w_end + wav->Tyy_pos;
  float *restrict Txy   = w_end + wav->Txy_pos;
  float *restrict Tzz   = w_end + wav->Tzz_pos;
  float *restrict Txz   = w_end + wav->Txz_pos;
  float *restrict Tyz   = w_end + wav->Tyz_pos;

  float *restrict zt_x  = metric->zeta_x;
  float *restrict zt_y  = metric->zeta_y;
  float *restrict zt_z  = metric->zeta_z;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *matD = bdryfree->matD;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk2 = gdinfo->nk2;

  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;

  float D[3][3], DT[3][3], Tl[3][3], DTl[3][3], Tg[3][3];
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
  float lam,mu,lam2mu;
  float tzz;

  size_t iptr_k = nk2 * siz_slice;
  for (size_t j=nj1; j<=nj2; j++)
  {
      size_t iptr_j = iptr_k + j * siz_line;
      size_t iptr = iptr_j + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
          size_t ij = (i + j * siz_line)*9;

        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        lam2mu = lam + 2.0 * mu;

          D[0][0] = matD[ij+3*0+0];
          D[0][1] = matD[ij+3*0+1];
          D[0][2] = matD[ij+3*0+2];
          D[1][0] = matD[ij+3*1+0];
          D[1][1] = matD[ij+3*1+1];
          D[1][2] = matD[ij+3*1+2];
          D[2][0] = matD[ij+3*2+0];
          D[2][1] = matD[ij+3*2+1];
          D[2][2] = matD[ij+3*2+2];

          DT[0][0] = D[0][0];
          DT[0][1] = D[1][0];
          DT[0][2] = D[2][0];
          DT[1][0] = D[0][1];
          DT[1][1] = D[1][1];
          DT[1][2] = D[2][1];
          DT[2][0] = D[0][2];
          DT[2][1] = D[1][2];
          DT[2][2] = D[2][2];

          d11 = D[0][0];
          d12 = D[0][1];
          d13 = D[0][2];
          d21 = D[1][0];
          d22 = D[1][1];
          d23 = D[1][2];
          d31 = D[2][0];
          d32 = D[2][1];
          d33 = D[2][2];

          Tl[0][0] =    d11*d11*Txx[iptr] + d12*d12*Tyy[iptr] + d13*d13*Tzz[iptr]
                   +2*(d11*d12*Txy[iptr] + d11*d13*Txz[iptr] + d12*d13*Tyz[iptr]);

          Tl[0][1] =    d11*d21*Txx[iptr] + d12*d22*Tyy[iptr] + d13*d23*Tzz[iptr]
                   +   (d11*d22+d12*d21)*Txy[iptr] + (d11*d23+d21*d13)*Txz[iptr] + (d12*d23+d22*d13)*Tyz[iptr];
         
          Tl[1][1] =    d21*d21*Txx[iptr] + d22*d22*Tyy[iptr] + d23*d23*Tzz[iptr]
                   +2*(d21*d22*Txy[iptr] + d21*d23*Txz[iptr] + d22*d23*Tyz[iptr]);

          Tl[1][0] = Tl[0][1];
          Tl[0][2] = 0.0;
          Tl[1][2] = 0.0;
          Tl[2][0] = 0.0;
          Tl[2][1] = 0.0;
          Tl[2][2] = 0.0;

          fdlib_math_matmul3x3(DT, Tl, DTl);
          fdlib_math_matmul3x3(DTl, D, Tg);

          Txx[iptr] = Tg[0][0];
          Tyy[iptr] = Tg[1][1];
          Tzz[iptr] = Tg[2][2];
          Txy[iptr] = Tg[0][1];
          Txz[iptr] = Tg[0][2];
          Tyz[iptr] = Tg[1][2];
          
          iptr+=1;
      }
  }
  return ierr;
}


               


