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
  float *matA     = bdry->matA;

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
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                        myid, verbose);

    // velocity: vlow
    sv_curv_col_vis_iso_rhs_vlow_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                         Jxx,Jyy,Jzz,Jxy,Jxz,Jyz,
                                         hJxx,hJyy,hJzz,hJxy,hJxz,hJyz,
                                         xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                         lam3d, mu3d, slw3d,
                                         wl, Ylam, Ymu,
                                         matVx2Vz,matVy2Vz,matA,
                                         ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                         fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                         fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                         num_of_fdz_op,fdz_op,fdz_max_len,
                                         nmaxwell,
                                         myid, verbose);
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
                                    src,
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
 * free surface boundary
 ******************************************************************************/

/*
 * implement vlow boundary
 */

void
sv_curv_col_vis_iso_rhs_vlow_z2(
    float *restrict  Vx ,  float *restrict  Vy ,  float *restrict  Vz ,
    float *restrict hTxx,  float *restrict hTyy,  float *restrict hTzz,
    float *restrict hTxz,  float *restrict hTyz,  float *restrict hTxy,
    float **restrict Jxx,  float **restrict Jyy,  float **restrict Jzz,
    float **restrict Jxy,  float **restrict Jxz,  float **restrict Jyz,
    float **restrict hJxx, float **restrict hJyy, float **restrict hJzz,
    float **restrict hJxy, float **restrict hJxz, float **restrict hJyz,
    float *restrict xi_x,  float *restrict xi_y,  float *restrict xi_z,
    float *restrict et_x,  float *restrict et_y,  float *restrict et_z,
    float *restrict zt_x,  float *restrict zt_y,  float *restrict zt_z,
    float *restrict lam3d, float *restrict mu3d,  float *restrict slw3d,
    float *restrict wl,    float **restrict Ylam, float **restrict Ymu,
    float *restrict matVx2Vz, float *restrict matVy2Vz, float *restrict matA,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int num_of_fdz_op, fd_op_t *fdz_op, int fdz_max_len,
    int nmaxwell,
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
  float mem_Txx,mem_Tyy,mem_Tzz,mem_Txy,mem_Txz,mem_Tyz;
  float sum_Jxyz,sum_Jxx,sum_Jyy,sum_Jzz,sum_Jxy,sum_Jxz,sum_Jyz;

  float A[3][3], M[3][3], AM[3][3];

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

          size_t ij = (i + j * siz_line)*9;

          // coef cal. first dim: irow; sec dim: jcol, as Fortran code
          A[0][0] = matA[ij+3+0+0]; 
          A[0][1] = matA[ij+3+0+1]; 
          A[0][2] = matA[ij+3+0+2]; 
          A[1][0] = matA[ij+3+1+0]; 
          A[1][1] = matA[ij+3+1+1]; 
          A[1][2] = matA[ij+3+1+2]; 
          A[2][0] = matA[ij+3+2+0]; 
          A[2][1] = matA[ij+3+2+1]; 
          A[2][2] = matA[ij+3+2+2]; 

          M[0][0] = mem_Txx;
          M[0][1] = mem_Txy;
          M[0][2] = mem_Txz;
          M[1][0] = mem_Txy;
          M[1][1] = mem_Tyy;
          M[1][2] = mem_Tyz;
          M[2][0] = mem_Txz;
          M[2][1] = mem_Tyz;
          M[2][2] = mem_Tzz;
          
          fdlib_math_matmul3x3(A, M, AM);

          DzVx = matVx2Vz[ij+3*0+0] * DxVx
               + matVx2Vz[ij+3*0+1] * DxVy
               + matVx2Vz[ij+3*0+2] * DxVz
               + matVy2Vz[ij+3*0+0] * DyVx
               + matVy2Vz[ij+3*0+1] * DyVy
               + matVy2Vz[ij+3*0+2] * DyVz
               + AM[0][0] * ztx
               + AM[0][1] * zty
               + AM[0][2] * ztz
               ;

          DzVy = matVx2Vz[ij+3*1+0] * DxVx
               + matVx2Vz[ij+3*1+1] * DxVy
               + matVx2Vz[ij+3*1+2] * DxVz
               + matVy2Vz[ij+3*1+0] * DyVx
               + matVy2Vz[ij+3*1+1] * DyVy
               + matVy2Vz[ij+3*1+2] * DyVz
               + AM[1][0] * ztx
               + AM[1][1] * zty
               + AM[1][2] * ztz
               ;

          DzVz = matVx2Vz[ij+3*2+0] * DxVx
               + matVx2Vz[ij+3*2+1] * DxVy
               + matVx2Vz[ij+3*2+2] * DxVz
               + matVy2Vz[ij+3*2+0] * DyVx
               + matVy2Vz[ij+3*2+1] * DyVy
               + matVy2Vz[ij+3*2+2] * DyVz
               + AM[2][0] * ztx
               + AM[2][1] * zty
               + AM[2][2] * ztz
               ;
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
  float lam,mu,lam2mu,slw;
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
        slw = slw3d[iptr];
        lam2mu = lam + 2.0 * mu;

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
sv_curv_col_vis_iso_dvh2dvz(gd_t        *gdinfo,
                                 gdcurv_metric_t *metric,
                                 md_t       *md,
                                 bdry_t      *bdryfree,
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
  float *matA = bdryfree->matA;
  
  float A[3][3], B[3][3], C[3][3];
  float AB[3][3], AC[3][3];

  float e11, e12, e13, e21, e22, e23, e31, e32, e33;
  float lam2mu, lam, mu;
 
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
          matA[ij + irow*3 + jcol] = A[irow][jcol]; //A is DZ invert
        }
    }
  }

  return ierr;
}
