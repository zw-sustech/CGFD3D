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
#include "sv_eq1st_cart_col_el_iso.h"

//#define SV_EQ1ST_CURV_COLGRD_ISO_DEBUG

/*******************************************************************************
 * perform one stage calculation of rhs
 ******************************************************************************/

void
sv_eq1st_cart_col_el_iso_onestage(
  float *restrict w_cur,
  float *restrict rhs, 
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  gd_t    *gdcart,
  md_t *md,
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

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;
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

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;

  float dx = gdcart->dx;
  float dy = gdcart->dy;
  float dz = gdcart->dz;

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

  // free surface at z2
  if (bdryfree->is_at_sides[2][1] == 1)
  {
    // tractiong
    sv_eq1st_cart_col_el_iso_rhs_timg_z2(Tzz,Txz,Tyz,
                                        ni1,ni2,nj1,nj2,nk1,nk2,nz,
                                        siz_line,siz_slice,
                                        myid, verbose);
  }

  // inner points
  sv_eq1st_cart_col_el_iso_rhs_inner(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                    hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    lam3d, mu3d, slw3d,
                                    dx, dy, dz,
                                    ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                    fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                    fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                    fdz_inn_len, fdz_inn_indx, fdz_inn_coef,
                                    myid, verbose);

  // free, abs, source in turn

  // free surface at z2
  if (bdryfree->is_at_sides[2][1] == 1)
  {
    // velocity: vlow
    sv_eq1st_cart_col_el_iso_rhs_vlow_z2(Vx,Vy,Vz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                        lam3d, mu3d, slw3d,
                                        matVx2Vz,matVy2Vz,
                                        dx, dy, dz,
                                        ni1,ni2,nj1,nj2,nk1,nk2,siz_line,siz_slice,
                                        fdx_inn_len, fdx_inn_indx, fdx_inn_coef,
                                        fdy_inn_len, fdy_inn_indx, fdy_inn_coef,
                                        num_of_fdz_op,fdz_op,fdz_max_len,
                                        myid, verbose);
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_cart_col_el_iso_rhs_cfspml(Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                       lam3d, mu3d, slw3d,
                                       dx, dy, dz,
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
    sv_eq1st_cart_col_el_iso_rhs_src(hVx,hVy,hVz,hTxx,hTyy,hTzz,hTxz,hTyz,hTxy,
                                    slw3d, dx, dy, dz,
                                    src,
                                    myid, verbose);
  }
  // end func

  return;
}

/*******************************************************************************
 * calculate all points without boundaries treatment
 ******************************************************************************/

void
sv_eq1st_cart_col_el_iso_rhs_inner(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float dx, float dy, float dz,
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
  float DxTxx,            DxTxy,DxTxz,      DxVx,DxVy,DxVz;
  float       DyTyy,      DyTxy,      DyTyz,DyVx,DyVy,DyVz;
  float             DzTzz,      DzTxz,DzTyz,DzVx,DzVy,DzVz;
  float lam,mu,mu2,slw,Eii;

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
    lfdx_coef [i] = fdx_coef[i] / dx;
    lfdx_shift[i] = fdx_indx[i];
  }
  for (int j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j] / dy;
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (int k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k] / dz;
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
        Txx_ptr = Txx + iptr;
        Tyy_ptr = Tyy + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;
        Tyz_ptr = Tyz + iptr;
        Txy_ptr = Txy + iptr;

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

        // Txx derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxx, Txx_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);

        // Tyy derivatives
        M_FD_SHIFT_PTR_MACDRP(DyTyy, Tyy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // Tzz derivatives
        M_FD_SHIFT_PTR_MACDRP(DzTzz, Tzz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txz derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxz, Txz_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTxz, Txz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Tyz derivatives
        M_FD_SHIFT_PTR_MACDRP(DyTyz, Tyz_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DzTyz, Tyz_ptr, fdz_len, lfdz_shift, lfdz_coef, n_fd);

        // Txy derivatives
        M_FD_SHIFT_PTR_MACDRP(DxTxy, Txy_ptr, fdx_len, lfdx_shift, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR_MACDRP(DyTxy, Txy_ptr, fdy_len, lfdy_shift, lfdy_coef, n_fd);

        // medium
        slw = slw3d[iptr];
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        // moment equation
        hVx[iptr] = slw*( DxTxx + DyTxy + DzTxz );
        hVy[iptr] = slw*( DxTxy + DyTyy + DzTyz );
        hVz[iptr] = slw*( DxTxz + DyTyz + DzTzz );

        // Hooke's equatoin
        Eii = lam * (DxVx + DyVy + DzVz);

        hTxx[iptr] = Eii + mu2 * DxVx;
        hTyy[iptr] = Eii + mu2 * DyVy;
        hTzz[iptr] = Eii + mu2 * DzVz;
        hTxz[iptr] = mu *( DxVz + DzVx );
        hTyz[iptr] = mu *( DyVz + DzVy );
        hTxy[iptr] = mu *( DxVy + DyVx );

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
sv_eq1st_cart_col_el_iso_rhs_timg_z2(
    float *restrict  Tzz, float *restrict  Txz, float *restrict  Tyz,
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
      Tzz[iptr] = 0.0;
      Tyz[iptr] = 0.0;
      Txz[iptr] = 0.0;

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

        Tzz[iptr_gho] = -Tzz[iptr_phy];
        Tyz[iptr_gho] = -Tyz[iptr_phy];
        Txz[iptr_gho] = -Txz[iptr_phy];
      }
    }
  }

  return;
}

/*
 * implement vlow boundary
 */

void
sv_eq1st_cart_col_el_iso_rhs_vlow_z2(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict lam3d, float *restrict mu3d, float *restrict slw3d,
    float *restrict matVx2Vz, float *restrict matVy2Vz,
    float dx, float dy, float dz,
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
  float lam,mu,mu2,slw,Eii;

  // put fd op into local array
  for (i=0; i < fdx_len; i++) {
    lfdx_coef [i] = fdx_coef[i] / dx;
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j] / dy;
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
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    }

    // for index
    size_t iptr_k = k * siz_slice;

    for (j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;

      size_t iptr = iptr_j + ni1;

      for (i=ni1; i<=ni2; i++)
      {
        // medium
        slw = slw3d[iptr];
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

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
        Eii = lam * (DxVx + DyVy + DzVz);

        hTxx[iptr] = Eii + mu2 * DxVx;
        hTyy[iptr] = Eii + mu2 * DyVy;
        hTzz[iptr] = Eii + mu2 * DzVz;
        hTxz[iptr] = mu *( DxVz + DzVx );
        hTyz[iptr] = mu *( DyVz + DzVy );
        hTxy[iptr] = mu *( DxVy + DyVx );

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
sv_eq1st_cart_col_el_iso_rhs_cfspml(
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict lam3d, float *restrict  mu3d, float *restrict slw3d,
    float dx, float dy, float dz,
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
    lfdx_coef [i] = fdx_coef[i] / dx;
    lfdx_shift[i] = fdx_indx[i];
  }
  for (j=0; j < fdy_len; j++) {
    lfdy_coef [j] = fdy_coef[j] / dy;
    lfdy_shift[j] = fdy_indx[j] * siz_line;
  }
  for (k=0; k < fdz_len; k++) {
    lfdz_coef [k] = fdz_coef[k] / dz;
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

              // medium
              slw = slw3d[iptr];
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
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
               hVx_rhs = slw * ( DxTxx  );
               hVy_rhs = slw * ( DxTxy  );
               hVz_rhs = slw * ( DxTxz  );
              hTxx_rhs = lam2mu * DxVx;
              hTyy_rhs = lam * DxVx;
              hTzz_rhs = lam * DxVx;
              hTxz_rhs = mu*( DxVz );
              hTxy_rhs = mu*( DxVy );

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
              hTxz[iptr] += coef_B_minus_1 * hTxz_rhs - coef_B * pml_Txz[iptr_a];
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
              pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

              // add contributions from free surface condition
              //  not consider timg because conflict with main cfspml,
              //     need to revise in the future if required
              if (bdryfree->is_at_sides[CONST_NDIM-1][1]==1 && k==nk2)
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
                hTxx_rhs = lam    * Dx_DzVz;

                hTyy_rhs = lam    * Dx_DzVz;

                hTzz_rhs = lam2mu * Dx_DzVz;

                hTxz_rhs = mu     * Dx_DzVx ;

                // make corr to Hooke's equatoin
                hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
                hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs;
                hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
                hTxz[iptr] += (coef_B - 1.0) * hTxz_rhs;

                // aux var
                //   a1 = alpha + d / beta, dealt in abs_set_cfspml
                pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
                pml_hTyy[iptr_a] += coef_D * hTyy_rhs;
                pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
                pml_hTxz[iptr_a] += coef_D * hTxz_rhs;
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
               hVx_rhs = slw * DyTxy;
               hVy_rhs = slw * DyTyy;
               hVz_rhs = slw * DyTyz;
              hTxx_rhs = lam * DyVy;
              hTyy_rhs = lam2mu * DyVy;
              hTzz_rhs = lam * DyVy;
              hTxy_rhs = mu  * DyVx;
              hTyz_rhs = mu  * DyVz;

              // 1: make corr to moment equation
              hVx[iptr] += coef_B_minus_1 * hVx_rhs - coef_B * pml_Vx[iptr_a];
              hVy[iptr] += coef_B_minus_1 * hVy_rhs - coef_B * pml_Vy[iptr_a];
              hVz[iptr] += coef_B_minus_1 * hVz_rhs - coef_B * pml_Vz[iptr_a];

              // make corr to Hooke's equatoin
              hTxx[iptr] += coef_B_minus_1 * hTxx_rhs - coef_B * pml_Txx[iptr_a];
              hTyy[iptr] += coef_B_minus_1 * hTyy_rhs - coef_B * pml_Tyy[iptr_a];
              hTzz[iptr] += coef_B_minus_1 * hTzz_rhs - coef_B * pml_Tzz[iptr_a];
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
              pml_hTyz[iptr_a] = coef_D * hTyz_rhs - coef_A * pml_Tyz[iptr_a];
              pml_hTxy[iptr_a] = coef_D * hTxy_rhs - coef_A * pml_Txy[iptr_a];

              // add contributions from free surface condition
              if (bdryfree->is_at_sides[CONST_NDIM-1][1]==1 && k==nk2)
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

                hTxx_rhs =   lam * Dy_DzVz;

                hTyy_rhs =   lam * Dy_DzVz;

                hTzz_rhs =   lam2mu * Dy_DzVz;

                hTyz_rhs = mu * Dy_DzVy;

                // make corr to Hooke's equatoin
                hTxx[iptr] += (coef_B - 1.0) * hTxx_rhs;
                hTyy[iptr] += (coef_B - 1.0) * hTyy_rhs;
                hTzz[iptr] += (coef_B - 1.0) * hTzz_rhs;
                hTyz[iptr] += (coef_B - 1.0) * hTyz_rhs;

                // aux var
                //   a1 = alpha + d / beta, dealt in abs_set_cfspml
                pml_hTxx[iptr_a] += coef_D * hTxx_rhs;
                pml_hTyy[iptr_a] += coef_D * hTyy_rhs;
                pml_hTzz[iptr_a] += coef_D * hTzz_rhs;
                pml_hTyz[iptr_a] += coef_D * hTyz_rhs;
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
               hVx_rhs = slw * DzTxz;
               hVy_rhs = slw * DzTyz;
               hVz_rhs = slw * DzTzz;
              hTxx_rhs = lam    * DzVz;
              hTyy_rhs = lam    * DzVz;
              hTzz_rhs = lam2mu * DzVz;
              hTxz_rhs = mu     * DzVx; 
              hTyz_rhs = mu     * DzVy;

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
 * free surface coef
 ******************************************************************************/

int
sv_eq1st_cart_col_el_iso_dvh2dvz(gdinfo_t   *gdinfo,
                                 md_t       *md,
                                 bdryfree_t      *bdryfree,
                                 const int verbose)
{
  int ierr = 0;

  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *matVx2Vz = bdryfree->matVx2Vz2;
  float *matVy2Vz = bdryfree->matVy2Vz2;
  
  float lam2mu, lam, mu;
 
  int k = nk2;

  for (size_t j = nj1; j <= nj2; j++)
  {
    for (size_t i = ni1; i <= ni2; i++)
    {
      size_t iptr = i + j * siz_line + k * siz_slice;

      lam    = lam3d[iptr];
      mu     =  mu3d[iptr];
      lam2mu = lam + 2.0f * mu;

      size_t ij = (j * siz_line + i) * 9;

      // save into mat
      for(int irow = 0; irow < 3; irow++)
        for(int jcol = 0; jcol < 3; jcol++){
          matVx2Vz[ij + irow*3 + jcol] = 0.0;
          matVy2Vz[ij + irow*3 + jcol] = 0.0;
        }

      // DzVx = -DxVz
      int DzVx_j = 0;
      int DxVz_i = 2;
      float coef = -1.0;
      matVx2Vz[ij + DzVx_j * CONST_NDIM + DxVz_i] = coef;

      // DzVy = -DyVz
      int DzVy_j = 1;
      int DyVz_i = 2;
      coef = -1.0;
      matVy2Vz[ij + DzVy_j * CONST_NDIM + DyVz_i] = coef;

      // DzVz = - (lam/lam2mu) * DxVx - (lam/lam2mu) * DyVy
      int DzVz_j = 2;
      int DxVx_i = 0;
      int DyVy_i = 1;
      coef = - lam / lam2mu;
      matVx2Vz[ij + DzVz_j * CONST_NDIM + DxVx_i] = coef;
      matVy2Vz[ij + DzVz_j * CONST_NDIM + DyVy_i] = coef;
    }
  }

  return ierr;
}

/*******************************************************************************
 * add source terms
 ******************************************************************************/

int
sv_eq1st_cart_col_el_iso_rhs_src(
    float *restrict hVx , float *restrict hVy , float *restrict hVz ,
    float *restrict hTxx, float *restrict hTyy, float *restrict hTzz,
    float *restrict hTxz, float *restrict hTyz, float *restrict hTxy,
    float *restrict slw3d,
    float dx, float dy, float dz,
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

  float vol = dx * dy * dz;

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
          float V = coef * slw3d[iptr] / vol;
          hVx[iptr] += fx * V;
          hVy[iptr] += fy * V;
          hVz[iptr] += fz * V;
        }

        if (src->moment_actived == 1) {
          float rjac = coef / vol;
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

