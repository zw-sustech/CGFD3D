/*******************************************************************************
 * solver of isotropic elastic 1st-order eqn using cart grid and staggerd scheme
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"
#include "sv_eq1st_cart_stg_ac_iso.h"

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 *  simple MPI exchange without computing-communication overlapping
 ******************************************************************************/

void
sv_eq1st_cart_stg_ac_iso_allstep(
  fdstg_t    *fd,
  gdinfo_t   *gdinfo,
  gd_t       *gdcart,
  md_t       *md,
  src_t      *src,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  wav_t      *wav,
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
  float dx = gdcart->dx;
  float dy = gdcart->dy;
  float dz = gdcart->dz;

  // mpi
  int myid = mympi->myid;
  int *topoid = mympi->topoid;
  MPI_Comm comm = mympi->comm;
  float *restrict sbuff = mympi->sbuff;
  float *restrict rbuff = mympi->rbuff;

  float t_cur;
  float t_end; // time after this loop for nc output

  // create slice nc output files
  if (myid==0 && verbose>0) fprintf(stdout,"prepare slice nc output ...\n"); 
  ioslice_nc_t ioslice_nc;
  io_slice_nc_create(ioslice, wav->ncmp, wav->cmp_name,
                     gdinfo->ni, gdinfo->nj, gdinfo->nk, topoid,
                     &ioslice_nc);

  // create snapshot nc output files
  if (myid==0 && verbose>0) fprintf(stdout,"prepare snap nc output ...\n"); 
  iosnap_nc_t  iosnap_nc;
  io_snap_nc_create_ac(iosnap, &iosnap_nc, topoid);

  // alloc working space for slice and snap ouotput
  size_t siz_max_wrk = iosnap->siz_max_wrk;
  if (ioslice->siz_max_wrk > siz_max_wrk) {
    siz_max_wrk = ioslice->siz_max_wrk;
  }
  float *restrict wrk = (float *) fdlib_mem_calloc_1d_float(siz_max_wrk,
                        0.0, "wrk in stg el iso");

  // only x/y mpi
  int num_of_r_reqs = 4;
  int num_of_s_reqs = 4;

  //--------------------------------------------------------
  // time loop
  //--------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"start time loop ...\n"); 

  for (int it=0; it<nt_total; it++)
  {
    //
    // Update Tij first, considering force Green function 
    //

    t_cur  = it * dt + t0;
    t_end = t_cur + 0.5 * dt;

    if (myid==0 && verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // set src_t time
    src_set_time(src, it, 0);

    sv_eq1st_cart_stg_ac_iso_hook(fd, wav,
        gdinfo, md, bdryfree, bdrypml, src,
        dt, dx, dy, dz,
        myid, verbose);

    //// tij recv mesg
    for (int n=0; n < mympi->siz_rbuff; n++) {
      rbuff[n] = 0.0;
    }
    MPI_Startall(num_of_r_reqs, mympi->r_reqs_stress);

    blk_stg_ac1st_pack_mesg_pressure(fd, gdinfo, wav, sbuff);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs_stress);

    // add IO here
    // write slice, use w_rhs as buff
    io_slice_nc_put(ioslice,&ioslice_nc,gdinfo,wav->v5d,wrk,it,t_end,
                     CONST_NDIM, wav->ncmp-1);

    // snapshot
    io_snap_nc_put(iosnap, &iosnap_nc, gdinfo, md, wav, 
                   wav->v5d, wrk, nt_total, it, t_end, 0, 1, 0);

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs_stress, MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs_stress, MPI_STATUS_IGNORE);

    //for (int n=0; n < mympi->siz_rbuff; n++) {
    //  rbuff[n] = 0.0;
    //}
    blk_stg_ac1st_unpack_mesg_pressure(fd,mympi,gdinfo,wav,rbuff,mympi->siz_rbuff);

    //
    // Update Vi
    //

    t_cur  = it * dt + t0 + 0.5 * dt;
    t_end = t_cur + 0.5 * dt;
    if (myid==0 && verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // set src_t time
    src_set_time(src, it, 0);

    // from Tij to Vi
    sv_eq1st_cart_stg_ac_iso_moment(fd,wav,
        gdinfo, md, bdryfree, bdrypml, src,
        dt, dx, dy, dz,
        myid, verbose);

    // recv mesg
    for (int n=0; n < mympi->siz_rbuff; n++) {
      rbuff[n] = 0.0;
    }
    MPI_Startall(num_of_r_reqs, mympi->r_reqs_vel);

    // pack and isend
    blk_stg_ac1st_pack_mesg_vel(fd,gdinfo,wav,sbuff);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs_vel);

    // add IO here

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs_vel, MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs_vel, MPI_STATUS_IGNORE);

    blk_stg_ac1st_unpack_mesg_vel(fd,mympi,gdinfo,wav,rbuff,mympi->siz_rbuff);

    //--------------------------------------------
    // QC
    //--------------------------------------------

    if (qc_check_nan_num_of_step >0  && (it % qc_check_nan_num_of_step) == 0) {
      if (myid==0 && verbose>10) fprintf(stdout,"-> check value nan\n");
        //wav_check_value(wav->v5d);
    }

    //--------------------------------------------
    // save results
    //--------------------------------------------

    //-- recv by interp
    io_recv_keep(iorecv, wav->v5d, it, wav->ncmp, wav->siz_icmp);

    //-- line values
    io_line_keep(ioline, wav->v5d, it, wav->ncmp, wav->siz_icmp);

    // write slice, 
    io_slice_nc_put(ioslice,&ioslice_nc,gdinfo,wav->v5d,wrk,it,t_end,
                     0, CONST_NDIM-1);

    // snapshot
    io_snap_nc_put_ac(iosnap, &iosnap_nc, gdinfo, wav, 
                   wav->v5d, wrk, nt_total, it, t_end, 1, 0, 1);

    // debug output
    if (output_all==1)
    {
      char ou_file[CONST_MAX_STRLEN];

        io_build_fname_time(output_dir,"w3d",".nc",topoid,it,ou_file);
        io_var3d_export_nc(ou_file,
                           wav->v5d,
                           wav->cmp_pos,
                           wav->cmp_name,
                           wav->ncmp,
                           gdinfo->index_name,
                           gdinfo->nx,
                           gdinfo->ny,
                           gdinfo->nz);
    }

  } // time loop

  // postproc

  // close nc
  io_slice_nc_close(&ioslice_nc);
  io_snap_nc_close(&iosnap_nc);

  if (wrk) free(wrk);

  return;
}

/*******************************************************************************
 * perform one stage calculation
 ******************************************************************************/

void
sv_eq1st_cart_stg_ac_iso_hook(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose)
{
  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  float *restrict P     = wav->v5d + wav->Txx_pos;

  float *restrict kappa3d = md->kappa;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nz  = gdinfo->nz ;

  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  // local var
  float kappa;
  float Eii;
  float DxVx,DyVx,DzVx;
  float DxVy,DyVy,DzVy;
  float DxVz,DyVz,DzVz;

  // allocate max_len because fdz may have different lens
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdy_coef   [fd->fdy_max_len];
  int    lfdy_shift_F[fd->fdy_max_len];
  int    lfdy_shift_B[fd->fdy_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  int n_fd;
  float *restrict Vx_ptr;
  float *restrict Vy_ptr;
  float *restrict Vz_ptr;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op

  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdy_op = fd->num_of_fdy_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  int   lfdy_len    = fd->lay_fdy_op[num_of_fdy_op-1].total_len;
  int   *p_fdy_indx = fd->lay_fdy_op[num_of_fdy_op-1].indx;
  float *p_fdy_coef = fd->lay_fdy_op[num_of_fdy_op-1].coef;
  for (n_fd = 0; n_fd < lfdy_len ; n_fd++)
  {
    lfdy_coef[n_fd]  = p_fdy_coef[n_fd] / dy;
    lfdy_shift_F[n_fd] = (p_fdy_indx[n_fd] + 1) * siz_line;
    lfdy_shift_B[n_fd] = (p_fdy_indx[n_fd] + 0) * siz_line;
  }

  // free surface at z2
  if (bdryfree->is_at_sides[2][1] == 1)
  {
    sv_eq1st_cart_stg_ac_iso_free_vimg(Vx,Vy,Vz,
                                  ni1,ni2,nj1,nj2,nk1,nk2,nz,
                                  siz_line,siz_slice,
                                  myid, verbose);
  }

  // at surface layers
  for (size_t n=0; n < num_of_fdz_op-1; n++)
  {
    // conver to k index, from surface to inner
    int k = nk2 - n;

    // use normal op if not free surface
    int lfdz_op_n = num_of_fdz_op - 1; 
    // use lower order free surface at z2
    if (bdryfree->is_at_sides[2][1] == 1) {
      lfdz_op_n = n;
    }

    // get pos and len for this point
    int  lfdz_len      = fd->lay_fdz_op[lfdz_op_n].total_len;
    int   *p_fdz_indx  = fd->lay_fdz_op[lfdz_op_n].indx;
    float *p_fdz_coef  = fd->lay_fdz_op[lfdz_op_n].coef;
    for (n_fd = 0; n_fd < lfdz_len ; n_fd++) {
      lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
      lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_slice;
      lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_slice;
    }

    size_t iptr_k = k * siz_slice;
    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      size_t iptr = iptr_j + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        kappa = kappa3d[iptr];

        Vx_ptr = Vx + iptr;
        Vy_ptr = Vy + iptr;
        Vz_ptr = Vz + iptr;

        // Derivatives for Tii
        M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyVy, Vy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // Hooke's equatoin
        P[iptr] -= dt * kappa * (DxVx + DyVy + DzVz);

        iptr += 1;
      }
    }
  }

  // for inner layers
  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]    = p_fdz_coef[n_fd] / dz;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_slice;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_slice;
  }

  //for (size_t k=nk1; k <= nk2; k++)
  for (size_t k=nk1; k <= nk2 - num_of_fdz_op + 1; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      size_t iptr = iptr_j + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        kappa = kappa3d[iptr];

        Vx_ptr = Vx + iptr;
        Vy_ptr = Vy + iptr;
        Vz_ptr = Vz + iptr;

        // Derivatives for Tii
        M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyVy, Vy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // Hooke's equatoin
        P[iptr] -= dt * kappa * (DxVx + DyVy + DzVz);

        iptr += 1;
      }
    }
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_cart_stg_ac_iso_hook_cfspml(fd,Vx,Vy,Vz,P,
                                       kappa3d, 
                                       dt, dx, dy, dz,
                                       nk2, siz_line,siz_slice,
                                       bdrypml, bdryfree,
                                       myid, verbose);
  }

  // add source term
  if (src->moment_actived == 1)
  {
    sv_eq1st_cart_stg_ac_iso_hook_src(P,
                                     dt, dx, dy, dz, src,
                                     myid, verbose);
  }
  // end func

  return;
}

void
sv_eq1st_cart_stg_ac_iso_moment(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose)
{
  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  float *restrict P    = wav->v5d + wav->Txx_pos;

  float *restrict slw3d = md->rho;

  // grid size
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  int nz  = gdinfo->nz ;

  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  // local var
  float slwdt;
  float DxP, DyP, DzP;

  // allocate max_len because fdz may have different lens
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdy_coef   [fd->fdy_max_len];
  int    lfdy_shift_F[fd->fdy_max_len];
  int    lfdy_shift_B[fd->fdy_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  int n_fd;
  float *restrict P_ptr;

  // for get a op from 1d array, currently use num_of_fdz_op as index
  // length, index, coef of a op

  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdy_op = fd->num_of_fdy_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  int   lfdy_len    = fd->lay_fdy_op[num_of_fdy_op-1].total_len;
  int   *p_fdy_indx = fd->lay_fdy_op[num_of_fdy_op-1].indx;
  float *p_fdy_coef = fd->lay_fdy_op[num_of_fdy_op-1].coef;
  for (n_fd = 0; n_fd < lfdy_len ; n_fd++)
  {
    lfdy_coef[n_fd]  = p_fdy_coef[n_fd] / dy;
    lfdy_shift_F[n_fd] = (p_fdy_indx[n_fd] + 1) * siz_line;
    lfdy_shift_B[n_fd] = (p_fdy_indx[n_fd] + 0) * siz_line;
  }

  // for inner layers
  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_slice;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_slice;
  }

  // free surface at z2
  if (bdryfree->is_at_sides[2][1] == 1)
  {
    sv_eq1st_cart_stg_ac_iso_free_simg(P,
                                  ni1,ni2,nj1,nj2,nk1,nk2,nz,
                                  siz_line,siz_slice,
                                  myid, verbose);
  }

  for (size_t k=nk1; k <= nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (size_t j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      size_t iptr = iptr_j + ni1;
      for (size_t i=ni1; i<=ni2; i++)
      {
        // medium
        slwdt = slw3d[iptr] * dt;

        P_ptr = P + iptr;

        // Tii deriv
        M_FD_SHIFT_PTR(DxP, P_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyP, P_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzP, P_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // Moment equatoin
        Vx[iptr] -= slwdt * DxP;
        Vy[iptr] -= slwdt * DyP;
        Vz[iptr] -= slwdt * DzP;

        iptr += 1;
      }
    }
  }

  // cfs-pml, loop face inside
  if (bdrypml->is_enable == 1)
  {
    sv_eq1st_cart_stg_ac_iso_moment_cfspml(fd,Vx,Vy,Vz,P,
                                       slw3d, dt, dx, dy, dz,
                                       nk2, siz_line,siz_slice,
                                       bdrypml,
                                       myid, verbose);
  }

  // add source term
  if (src->force_actived == 1)
  {
    sv_eq1st_cart_stg_ac_iso_moment_src(Vx, Vy, Vz,
                                     slw3d, dt, dx, dy, dz, src,
                                    myid, verbose);
  }
  // end func

  return;
}


/*******************************************************************************
 * CFS-PML boundary
 ******************************************************************************/

/*
 * updating Tij
 */

void
sv_eq1st_cart_stg_ac_iso_hook_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  P,
    float *restrict kappa3d, 
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdrypml_t *bdrypml, bdryfree_t *bdryfree,
    const int myid, const int verbose)
{
  // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdy_coef   [fd->fdy_max_len];
  int    lfdy_shift_F[fd->fdy_max_len];
  int    lfdy_shift_B[fd->fdy_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  int n_fd;
  float *restrict Vx_ptr;
  float *restrict Vy_ptr;
  float *restrict Vz_ptr;

  // val on point
  float kappa;

  // local
  int i,j,k;
  int iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_rB_minus_1;

  // put fd op into local array
  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdy_op = fd->num_of_fdy_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  int   lfdy_len    = fd->lay_fdy_op[num_of_fdy_op-1].total_len;
  int   *p_fdy_indx = fd->lay_fdy_op[num_of_fdy_op-1].indx;
  float *p_fdy_coef = fd->lay_fdy_op[num_of_fdy_op-1].coef;
  for (n_fd = 0; n_fd < lfdy_len ; n_fd++)
  {
    lfdy_coef[n_fd]  = p_fdy_coef[n_fd] / dy;
    lfdy_shift_F[n_fd] = (p_fdy_indx[n_fd] + 1) * siz_line;
    lfdy_shift_B[n_fd] = (p_fdy_indx[n_fd] + 0) * siz_line;
  }

  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_slice;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_slice;
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

      // get coef for this face
      float *restrict ptr_coef_A = bdrypml->A[idim][iside];
      float *restrict ptr_coef_B = bdrypml->B[idim][iside];
      float *restrict ptr_coef_D = bdrypml->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);

      // for each dim
      if (idim == 0 ) // x direction
      {
        float DxVx;
        // get pml vars, some vars are not used
        float *restrict pml_DxVx   = auxvar->var + auxvar->Vx_pos;

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
              coef_rB_minus_1 = coef_B - 1.0;

              // medium
              kappa = kappa3d[iptr];

              Vx_ptr = Vx + iptr;

              // for x deriv
              M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;

              // correct P
              float Vx1 = (2.0 * pml_DxVx[iptr_a] + dt* coef_D *DxVx ) * r_dt_AD;
              float Vx1_corr = coef_rB_minus_1 * DxVx - Vx1 * coef_B;
              P[iptr] -= dt * kappa * Vx1_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DxVx[iptr_a]  = 2.0 * Vx1  - pml_DxVx[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else if (idim == 1) // y direction
      {
        float DyVy;
        // get pml vars, some vars are not used
        float *restrict pml_DyVy   = auxvar->var + auxvar->Vy_pos;

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
            coef_rB_minus_1 = coef_B - 1.0;

            float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;

            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              kappa = kappa3d[iptr];

              Vy_ptr = Vy + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DyVy, Vy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);

              // correct Txx,Tyy,Tzz
              float Vy2 =  (2.0 * pml_DyVy[iptr_a] + dt * coef_D * DyVy ) * r_dt_AD;
              float Vy2_corr = coef_rB_minus_1 * DyVy - Vy2 * coef_B;
              P[iptr] -= dt * kappa * Vy2_corr;
              pml_DyVy[iptr_a] = 2.0 * Vy2 - pml_DyVy[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else // z direction
      {
        float DzVz;

        // get pml vars, some vars are not used
        float *restrict pml_DzVz   = auxvar->var + auxvar->Vz_pos;

        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_slice;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_rB_minus_1 = coef_B - 1.0;

          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              kappa = kappa3d[iptr];

              Vz_ptr = Vz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;

              // correct Txx,Tyy,Tzz
              float Vz3 =  (2.0 * pml_DzVz[iptr_a] + dt * coef_D * DzVz ) * r_dt_AD;
              float Vz3_corr = coef_rB_minus_1 * DzVz - Vz3 * coef_B;

              P[iptr] -= dt * kappa * Vz3_corr;

              pml_DzVz[iptr_a] = 2.0 * Vz3 - pml_DzVz[iptr_a];

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

void
sv_eq1st_cart_stg_ac_iso_moment_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  P,
    float *restrict slw3d,
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdrypml_t *bdrypml,
    const int myid, const int verbose)
{
  // loop var for fd
  int n_fd; // loop var for fd
  // use local stack array for better speed
  float  lfdx_coef   [fd->fdx_max_len];
  int    lfdx_shift_F[fd->fdx_max_len];
  int    lfdx_shift_B[fd->fdx_max_len];

  float  lfdy_coef   [fd->fdy_max_len];
  int    lfdy_shift_F[fd->fdy_max_len];
  int    lfdy_shift_B[fd->fdy_max_len];

  float  lfdz_coef   [fd->fdz_max_len];
  int    lfdz_shift_F[fd->fdz_max_len];
  int    lfdz_shift_B[fd->fdz_max_len];

  // local
  int i,j,k;
  int iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_rB_minus_1;
  float slw;

  float *restrict P_ptr;

  // put fd op into local array

  int num_of_fdx_op = fd->num_of_fdx_op;
  int num_of_fdy_op = fd->num_of_fdy_op;
  int num_of_fdz_op = fd->num_of_fdz_op;

  int   lfdx_len    = fd->lay_fdx_op[num_of_fdx_op-1].total_len;
  int   *p_fdx_indx = fd->lay_fdx_op[num_of_fdx_op-1].indx;
  float *p_fdx_coef = fd->lay_fdx_op[num_of_fdx_op-1].coef;
  for (n_fd = 0; n_fd < lfdx_len ; n_fd++)
  {
    lfdx_coef[n_fd]  = p_fdx_coef[n_fd] / dx;
    lfdx_shift_F[n_fd] = (p_fdx_indx[n_fd] + 1);
    lfdx_shift_B[n_fd] = (p_fdx_indx[n_fd]    );
  }

  int   lfdy_len    = fd->lay_fdy_op[num_of_fdy_op-1].total_len;
  int   *p_fdy_indx = fd->lay_fdy_op[num_of_fdy_op-1].indx;
  float *p_fdy_coef = fd->lay_fdy_op[num_of_fdy_op-1].coef;
  for (n_fd = 0; n_fd < lfdy_len ; n_fd++)
  {
    lfdy_coef[n_fd]  = p_fdy_coef[n_fd] / dy;
    lfdy_shift_F[n_fd] = (p_fdy_indx[n_fd] + 1) * siz_line;
    lfdy_shift_B[n_fd] = (p_fdy_indx[n_fd] + 0) * siz_line;
  }

  int   lfdz_len    = fd->lay_fdz_op[num_of_fdz_op-1].total_len;
  int   *p_fdz_indx = fd->lay_fdz_op[num_of_fdz_op-1].indx;
  float *p_fdz_coef = fd->lay_fdz_op[num_of_fdz_op-1].coef;
  for (n_fd = 0; n_fd < lfdz_len ; n_fd++)
  {
    lfdz_coef[n_fd]  = p_fdz_coef[n_fd] / dz;
    lfdz_shift_F[n_fd] = (p_fdz_indx[n_fd] + 1) * siz_slice;
    lfdz_shift_B[n_fd] = (p_fdz_indx[n_fd] + 0) * siz_slice;
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

      // get coef for this face
      float *restrict ptr_coef_A = bdrypml->A[idim][iside];
      float *restrict ptr_coef_B = bdrypml->B[idim][iside];
      float *restrict ptr_coef_D = bdrypml->D[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdrypml->auxvar[idim][iside]);

      // for each dim
      if (idim == 0 ) // x direction
      {
        float DxP;

        // get pml vars, some vars are not used
        float *restrict pml_DxP   = auxvar->var + auxvar->Txx_pos;

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
              coef_rB_minus_1 = coef_B - 1.0;

              // medium
              slw = slw3d[iptr];

              P_ptr = P + iptr;

              //  x deriv
              M_FD_SHIFT_PTR(DxP, P_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;

              // correct Vx
              float Txx1 = (2.0 * pml_DxP[iptr_a] + dt* coef_D *DxP ) * r_dt_AD;
              float Txx1_corr = coef_rB_minus_1 * DxP - Txx1 * coef_B;
              Vx[iptr] -= dt * slw * Txx1_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DxP[iptr_a]  = 2.0 * Txx1  - pml_DxP[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else if (idim == 1) // y direction
      {
        float DyP;

        // get pml vars, some vars are not used
        float *restrict pml_DyP   = auxvar->var + auxvar->Txx_pos;

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
            coef_rB_minus_1 = coef_B - 1.0;

            float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;

            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              slw = slw3d[iptr];

              P_ptr = P + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DyP, P_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);

              // correct Vy
              float Tyy2 = (2.0 * pml_DyP[iptr_a] + dt* coef_D *DyP ) * r_dt_AD;
              float Tyy2_corr = coef_rB_minus_1 * DyP - Tyy2 * coef_B;
              Vy[iptr] -= dt * slw * Tyy2_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DyP[iptr_a]  = 2.0 * Tyy2  - pml_DyP[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else // z direction
      {
        float DzP;

        // get pml vars, some vars are not used
        float *restrict pml_DzP   = auxvar->var + auxvar->Txx_pos;

        iptr_a = 0;
        for (k=abs_nk1; k<=abs_nk2; k++)
        {
          iptr_k = k * siz_slice;

          // pml coefs
          int abs_k = k - abs_nk1;
          coef_D = ptr_coef_D[abs_k];
          coef_A = ptr_coef_A[abs_k];
          coef_B = ptr_coef_B[abs_k];
          coef_rB_minus_1 = coef_B - 1.0;

          float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;

          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              slw = slw3d[iptr];

              P_ptr = P + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DzP, P_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

              // correct Vz
              float Tzz3 = (2.0 * pml_DzP[iptr_a] + dt* coef_D *DzP ) * r_dt_AD;
              float Tzz3_corr = coef_rB_minus_1 * DzP - Tzz3 * coef_B;
              Vz[iptr] -= dt * slw * Tzz3_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DzP[iptr_a]  = 2.0 * Tzz3  - pml_DzP[iptr_a];

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
sv_eq1st_cart_stg_ac_iso_hook_src(
    float *restrict P,
    float dt, float dx, float dy, float dz,
    src_t *src, // short nation for reference member
    const int myid, const int verbose)
{
  int ierr = 0;

  // local var
  int si,sj,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi / mij
  float Mii;

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

      Mii = src->Mxx[iptr_cur_stage];
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        float rjac = dt * coef / vol;

        P[iptr] -= Mii * rjac;
      } // i_ext

    } // it
  } // is

  return ierr;
}

int
sv_eq1st_cart_stg_ac_iso_moment_src(
    float *restrict Vx, float *restrict Vy, float *restrict Vz,
    float *restrict slw3d, float dt, float dx, float dy, float dz,
    src_t *src, // short nation for reference member
    const int myid, const int verbose)
{
  int ierr = 0;

  // local var
  int si,sj,sk, iptr;

  // for easy coding and efficiency
  int max_ext = src->max_ext;

  // get fi 
  float fx, fy, fz;

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

      fx  = src->Fx [iptr_cur_stage];
      fy  = src->Fy [iptr_cur_stage];
      fz  = src->Fz [iptr_cur_stage];
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        float V = dt * coef * slw3d[iptr] / vol;

        Vx[iptr] += fx * V;
        Vy[iptr] += fy * V;
        Vz[iptr] += fz * V;

      } // i_ext

    } // it
  } // is

  return ierr;
}

/*******************************************************************************
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_eq1st_cart_stg_ac_iso_free_simg(
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
 * image Vx Vy above free surface as Vx=0, Vy=0
 *    note Vz != 0, set Vz to zero above surface
 */

void
sv_eq1st_cart_stg_ac_iso_free_vimg(
    float *restrict  Vx, float *restrict  Vy, float *restrict  Vz,
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
      Vx[iptr] = 0.0;
      Vy[iptr] = 0.0;

      Vz[iptr] = 0.0;

      // next
      iptr += 1;
    }
  }

  // above surface
  for (size_t k=nk2+1; k<nz; k++)
  {
    int k_phy = nk2 - (k-nk2);
    for (size_t j=nj1; j<=nj2; j++)
    {
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr_gho = i + j * siz_line + k     * siz_slice;
        size_t iptr_phy = i + j * siz_line + k_phy * siz_slice;

        Vz[iptr_gho] = 0.0;

        Vx[iptr_gho] = -Vx[iptr_phy];
        Vy[iptr_gho] = -Vy[iptr_phy];
      }
    }
  }

  return;
}
