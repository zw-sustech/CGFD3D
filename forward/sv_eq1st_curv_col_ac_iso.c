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
#include "sv_eq1st_curv_col_ac_iso.h"

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eq1st_curv_col_ac_iso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  gdcurv_metric_t  *metric,
  md_t *md,
  bdry_t    *bdry,
  src_t *src,
  // include different order/stentil
  int num_of_fdx_op, fd_op_t *fdx_op,
  int num_of_fdy_op, fd_op_t *fdy_op,
  int num_of_fdz_op, fd_op_t *fdz_op,
  int fdz_max_len, 
  const int myid, const int verbose)
{
  // local pointer get each vars
  float *restrict Vx    = w_cur + wav->Vx_pos ;
  float *restrict Vy    = w_cur + wav->Vy_pos ;
  float *restrict Vz    = w_cur + wav->Vz_pos ;
  float *restrict P     = w_cur + wav->Txx_pos;
  float *restrict hVx   = rhs   + wav->Vx_pos ; 
  float *restrict hVy   = rhs   + wav->Vy_pos ; 
  float *restrict hVz   = rhs   + wav->Vz_pos ; 
  float *restrict hP    = rhs   + wav->Txx_pos; 

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

  float *restrict kappa3d = md->kappa;
  float *restrict slw3d = md->rho;

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

  // free surface at z2 for pressure
  if (bdry->is_sides_free[2][1] == 1)
  {
    // imaging
    sv_eq1st_curv_col_ac_iso_rhs_timg_z2(P,
                                     ni1,ni2,nj1,nj2,nk1,nk2,nz,
                                     siz_line,siz_slice,
                                     myid, verbose);
  }

  // inner points
  sv_eq1st_curv_col_ac_iso_rhs_inner(Vx,Vy,Vz,P,
                                    hVx,hVy,hVz,hP,
                                    xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                    kappa3d, slw3d,
                                    ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                    fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                    fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                    fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                    myid, verbose);

  // free, abs, source in turn

  // free surface at z2 for velocity
  if (bdry->is_sides_free[2][1] == 1)
  {
    // velocity: vlow
    sv_eq1st_curv_col_ac_iso_rhs_vlow_z2(Vx,Vy,Vz,hP,
                                        xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                        kappa3d, slw3d,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        num_of_fdz_op,fdz_op,fdz_max_len,
                                        myid, verbose);
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_eq1st_curv_col_ac_iso_rhs_cfspml(Vx,Vy,Vz,P,
                                    hVx,hVy,hVz,hP,
                                       xi_x, xi_y, xi_z, et_x, et_y, et_z, zt_x, zt_y, zt_z,
                                       kappa3d, slw3d,
                                       nk2, siz_line,siz_slice,
                                       fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                       fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                       fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                       bdry,
                                       myid, verbose);
    
  }

  // add source term
  if (src->total_number > 0)
  {
    sv_eq1st_curv_col_ac_iso_rhs_src(hVx,hVy,hVz,hP,
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
sv_eq1st_curv_col_ac_iso_rhs_inner(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  P, 
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hP, 
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict kappa3d, float *restrict slw3d,
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
  float DxP,DxVx,DxVy,DxVz;
  float DyP,DyVx,DyVy,DyVz;
  float DzP,DzVx,DzVy,DzVz;
  float kappa,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;

  float *restrict Vx_ptr;
  float *restrict Vy_ptr;
  float *restrict Vz_ptr;
  float *restrict P_ptr;

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
        P_ptr  = P  + iptr;

        // Vx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVx, Vx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyVx, Vx_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVx, Vx_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVy, Vy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyVy, Vy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVy, Vy_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxVz, Vz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyVz, Vz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzVz, Vz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // P derivatives
        M_FD_SHIFT_PTR_MACDRP(DxP, P_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyP, P_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzP, P_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

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
        kappa = kappa3d[iptr];
        slw = slw3d[iptr];

        // moment equation
        hVx[iptr] = - slw*( xix*DxP  
                           +etx*DyP 
                           +ztx*DzP );
        hVy[iptr] = - slw*( xiy*DxP
                           +ety*DyP
                           +zty*DzP );
        hVz[iptr] = - slw*( xiz*DxP 
                           +etz*DyP
                           +ztz*DzP );

        // Hooke's equatoin
        hP[iptr] = -kappa *  ( xix*DxVx  +etx*DyVx + ztx*DzVx
                              +xiy*DxVy + ety*DyVy + zty*DzVy
                              +xiz*DxVz + etz*DyVz + ztz*DzVz);

        iptr += 1;
      }
    }
  }

  return;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_eq1st_curv_col_ac_iso_rhs_timg_z2(
    float *restrict  P,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2, int nz,
    size_t siz_line, size_t siz_slice,
    const int myid, const int verbose)
{
  // nk2
  size_t iptr_k = nk2 * siz_slice;
  for (size_t j=nj1; j<=nj2; j++)
  {
    size_t iptr_j = iptr_k + j * siz_line;
    size_t iptr = iptr_j + ni1;
    for (size_t i=ni1; i<=ni2; i++)
    {
      P[iptr] = 0.0;

      // next
      iptr += 1;
    }
  }

  // mirror point
  for (size_t k=nk2+1; k<nz; k++)
  {
    int k_phy = nk2 - (k-nk2);
    for (size_t j=nj1; j<=nj2; j++)
    {
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr_gho = i + j * siz_line + k     * siz_slice;
        size_t iptr_phy = i + j * siz_line + k_phy * siz_slice;

        P[iptr_gho] = -P[iptr_phy];
      }
    }
  }

  return;
}

/*
 * implement vlow boundary
 */

void
sv_eq1st_curv_col_ac_iso_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hP, 
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict kappa3d, float *restrict slw3d,
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
  float kappa;
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
        kappa = kappa3d[iptr];

        // Vx derivatives
        M_FD_SHIFT(DxVx, Vx, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVx, Vx, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Vy derivatives
        M_FD_SHIFT(DxVy, Vy, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVy, Vy, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Vz derivatives
        M_FD_SHIFT(DxVz, Vz, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT(DyVz, Vz, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        if (k==nk2) // at surface, zero
        {
          hP[iptr] =  0.0;
        }
        else // lower than surface, lower order
        {
          M_FD_SHIFT(DzVx, Vx, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVy, Vy, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);
          M_FD_SHIFT(DzVz, Vz, iptr, lfdz_len, lfdz_shift, lfdz_coef, n_fd);

          hP[iptr] = -kappa *  ( xix*DxVx  +etx*DyVx + ztx*DzVx
                                +xiy*DxVy + ety*DyVy + zty*DzVy
                                +xiz*DxVz + etz*DyVz + ztz*DzVz);
        }

        iptr += 1;
      }
    }
  }

  return;
}

/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * cfspml, reference to each pml var inside function
 */

void
sv_eq1st_curv_col_ac_iso_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  P, 
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hP,
    float *restrict xi_x, float *restrict xi_y, float *restrict xi_z,
    float *restrict et_x, float *restrict et_y, float *restrict et_z,
    float *restrict zt_x, float *restrict zt_y, float *restrict zt_z,
    float *restrict kappa3d, float *restrict slw3d,
    int nk2, size_t siz_line, size_t siz_slice,
    int fdx_len, int *restrict fdx_indx, float *restrict fdx_coef,
    int fdy_len, int *restrict fdy_indx, float *restrict fdy_coef,
    int fdz_len, int *restrict fdz_indx, float *restrict fdz_coef,
    bdry_t    *bdry,
    const int myid, const int verbose)
{
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
  float DxP,DxVx,DxVy,DxVz;
  float DyP,DyVx,DyVy,DyVz;
  float DzP,DzVx,DzVy,DzVz;
  float kappa,slw;
  float xix,xiy,xiz,etx,ety,etz,ztx,zty,ztz;
  float hVx_rhs,hVy_rhs,hVz_rhs;
  float hP_rhs;

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
      if (bdry->is_sides_pml[idim][iside] == 0) continue;

      // get index into local var
      int abs_ni1 = bdry->ni1[idim][iside];
      int abs_ni2 = bdry->ni2[idim][iside];
      int abs_nj1 = bdry->nj1[idim][iside];
      int abs_nj2 = bdry->nj2[idim][iside];
      int abs_nk1 = bdry->nk1[idim][iside];
      int abs_nk2 = bdry->nk2[idim][iside];

      // get coef for this face
      float *restrict ptr_coef_A = bdry->A[idim][iside];
      float *restrict ptr_coef_B = bdry->B[idim][iside];
      float *restrict ptr_coef_D = bdry->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);

      // get pml vars
      float *restrict abs_vars_cur = auxvar->cur;
      float *restrict abs_vars_rhs = auxvar->rhs;

      float *restrict pml_Vx   = abs_vars_cur + auxvar->Vx_pos;
      float *restrict pml_Vy   = abs_vars_cur + auxvar->Vy_pos;
      float *restrict pml_Vz   = abs_vars_cur + auxvar->Vz_pos;
      float *restrict pml_P    = abs_vars_cur + auxvar->Txx_pos;

      float *restrict pml_hVx  = abs_vars_rhs + auxvar->Vx_pos;
      float *restrict pml_hVy  = abs_vars_rhs + auxvar->Vy_pos;
      float *restrict pml_hVz  = abs_vars_rhs + auxvar->Vz_pos;
      float *restrict pml_hP   = abs_vars_rhs + auxvar->Txx_pos;

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
              kappa = kappa3d[iptr];
              slw = slw3d[iptr];

              // xi derivatives
              M_FD_SHIFT(DxVx , Vx , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVy , Vy , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxVz , Vz , iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
              M_FD_SHIFT(DxP, P, iptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = -slw * ( xix*DxP                         );
               hVy_rhs = -slw * (             xiy*DxP             );
               hVz_rhs = -slw * (                         xiz*DxP );
                hP_rhs = -kappa*( xix*DxVx + xiy*DxVy + xiz*DxVz  );

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hP[iptr] += coef_B_minus_1 * hP_rhs - coef_B * pml_P[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hP[iptr_a] = coef_D * hP_rhs - coef_A * pml_P[iptr_a];

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
              kappa = kappa3d[iptr];
              slw = slw3d[iptr];

              // et derivatives
              M_FD_SHIFT(DyVx , Vx , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyVy , Vy , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyVz , Vz , iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
              M_FD_SHIFT(DyP, P, iptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = -slw * ( etx*DyP                         );
               hVy_rhs = -slw * (             ety*DyP             );
               hVz_rhs = -slw * (                         etz*DyP );
                hP_rhs = -kappa*(etx*DyVx + ety*DyVy + etz*DyVz);

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hP[iptr] += coef_B_minus_1 * hP_rhs - coef_B * pml_P[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hP[iptr_a] = coef_D * hP_rhs - coef_A * pml_P[iptr_a];

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
              slw = slw3d[iptr];
              kappa = kappa3d[iptr];

              // zt derivatives
              M_FD_SHIFT(DzVx , Vx , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVy , Vy , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzVz , Vz , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);
              M_FD_SHIFT(DzP , P , iptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

              // combine for corr and aux vars
               hVx_rhs = -slw * ( ztx*DzP                         );
               hVy_rhs = -slw * (             zty*DzP             );
               hVz_rhs = -slw * (                         ztz*DzP );
                hP_rhs = -kappa*(ztx*DzVx + zty*DzVy + ztz*DzVz);

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hP[iptr] += coef_B_minus_1 * hP_rhs - coef_B * pml_P[iptr_a];
              
              // 2: aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_hVx[iptr_a]  = coef_D * hVx_rhs  - coef_A * pml_Vx[iptr_a];
              pml_hVy[iptr_a]  = coef_D * hVy_rhs  - coef_A * pml_Vy[iptr_a];
              pml_hVz[iptr_a]  = coef_D * hVz_rhs  - coef_A * pml_Vz[iptr_a];
              pml_hP[iptr_a] = coef_D * hP_rhs - coef_A * pml_P[iptr_a];

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
 * add source terms
 ******************************************************************************/

int
sv_eq1st_curv_col_ac_iso_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hP, 
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
  float Mii;

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
        Mii = src->Mxx[iptr_cur_stage];
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
          hP[iptr] -= Mii * rjac;
        }
      } // i_ext

    } // it
  } // is

  return ierr;
}

