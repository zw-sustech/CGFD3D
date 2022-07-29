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
#include "sv_eq1st_cart_stg_el_iso.h"

/*******************************************************************************
 * one simulation over all time steps, could be used in imaging or inversion
 *  simple MPI exchange without computing-communication overlapping
 ******************************************************************************/

void
sv_eq1st_cart_stg_el_iso_allstep(
  fdstg_t    *fd,
  gdinfo_t   *gdinfo,
  gd_t       *gdcart,
  md_t       *md,
  src_t      *src,
  bdry_t    *bdry,
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
  io_snap_nc_create(iosnap, &iosnap_nc, topoid);

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

  // calculate conversion matrix for free surface
  if (bdry->is_sides_free[CONST_NDIM-1][1] == 1)
  {
     sv_eq1st_cart_stg_el_iso_dvh2dvz(gdinfo,md,bdry,verbose);
  }

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

    sv_eq1st_cart_stg_el_iso_hook(fd, wav,
        gdinfo, md, bdry, src,
        dt, dx, dy, dz,
        myid, verbose);

    //// tij recv mesg
    for (int n=0; n < mympi->siz_rbuff; n++) {
      rbuff[n] = 0.0;
    }
    MPI_Startall(num_of_r_reqs, mympi->r_reqs_stress);

    // pack and isend
    //for (int k=gdinfo->nk1; k <= gdinfo->nk2; k++) {
    //for (int j=gdinfo->nj1; j <= gdinfo->nj2; j++) {
    //for (int i=gdinfo->ni1; i <= gdinfo->ni2; i++) {
    //  int iptr = i + j * gdinfo->siz_iy + k * gdinfo->siz_iz;
    //  //wav->v5d[ wav->Txx_pos + iptr ] = (myid + 1) * 10000 + 1000 + i;
    //  //wav->v5d[ wav->Tyy_pos + iptr ] = (myid + 1) * 10000 + 2000 + i;
    //  //wav->v5d[ wav->Tzz_pos + iptr ] = (myid + 1) * 10000 + 3000 + i;
    //  //wav->v5d[ wav->Tyz_pos + iptr ] = (myid + 1) * 10000 + 4000 + i;
    //  //wav->v5d[ wav->Txz_pos + iptr ] = (myid + 1) * 10000 + 5000 + i;
    //  //wav->v5d[ wav->Txy_pos + iptr ] = (myid + 1) * 10000 + 6000 + i;
    //  wav->v5d[ wav->Txx_pos + iptr ] = iptr + (myid + 1) * 100000;
    //  wav->v5d[ wav->Tyy_pos + iptr ] = iptr + (myid + 1) * 200000;
    //  wav->v5d[ wav->Tzz_pos + iptr ] = iptr + (myid + 1) * 300000;
    //  wav->v5d[ wav->Tyz_pos + iptr ] = iptr + (myid + 1) * 400000;
    //  wav->v5d[ wav->Txz_pos + iptr ] = iptr + (myid + 1) * 500000;
    //  wav->v5d[ wav->Txy_pos + iptr ] = iptr + (myid + 1) * 600000;
    //}
    //}
    //}
    blk_stg_el1st_pack_mesg_stress(fd, gdinfo, wav, sbuff);

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
    blk_stg_el1st_unpack_mesg_stress(fd,mympi,gdinfo,wav,rbuff,mympi->siz_rbuff);

    //
    // Update Vi
    //

    t_cur  = it * dt + t0 + 0.5 * dt;
    t_end = t_cur + 0.5 * dt;
    if (myid==0 && verbose>10) fprintf(stdout,"-> it=%d, t=%f\n", it, t_cur);

    // set src_t time
    src_set_time(src, it, 0);

    // from Tij to Vi
    sv_eq1st_cart_stg_el_iso_moment(fd,wav,
        gdinfo, md, bdry, src,
        dt, dx, dy, dz,
        myid, verbose);

    // recv mesg
    for (int n=0; n < mympi->siz_rbuff; n++) {
      rbuff[n] = 0.0;
    }
    MPI_Startall(num_of_r_reqs, mympi->r_reqs_vel);

    // pack and isend
    //for (int k=gdinfo->nk1; k <= gdinfo->nk2; k++) {
    //for (int j=gdinfo->nj1; j <= gdinfo->nj2; j++) {
    //for (int i=gdinfo->ni1; i <= gdinfo->ni2; i++) {
    //  int iptr = i + j * gdinfo->siz_iy + k * gdinfo->siz_iz;
    //  //wav->v5d[ wav->Vx_pos + iptr ] = (myid + 1) * 20000 + 1000 + i;
    //  //wav->v5d[ wav->Vy_pos + iptr ] = (myid + 1) * 20000 + 2000 + i;
    //  //wav->v5d[ wav->Vz_pos + iptr ] = (myid + 1) * 20000 + 3000 + i;
    //  wav->v5d[ wav->Vx_pos + iptr ] = iptr + (myid + 1) * 10000; 
    //  wav->v5d[ wav->Vy_pos + iptr ] = iptr + (myid + 1) * 20000; 
    //  wav->v5d[ wav->Vz_pos + iptr ] = iptr + (myid + 1) * 30000; 
    //}
    //}
    //}
    blk_stg_el1st_pack_mesg_vel(fd,gdinfo,wav,sbuff);

    MPI_Startall(num_of_s_reqs, mympi->s_reqs_vel);

    // add IO here

    MPI_Waitall(num_of_s_reqs, mympi->s_reqs_vel, MPI_STATUS_IGNORE);
    MPI_Waitall(num_of_r_reqs, mympi->r_reqs_vel, MPI_STATUS_IGNORE);
    //int siz_x_per_lay = gdinfo->nj * gdinfo->nk;
    //int siz_y_per_lay = gdinfo->ni * gdinfo->nk;
    //int siz_to_x1 = siz_x_per_lay * ( fd->fdx_nghosts  - 1                // Vx
    //                                 +fd->fdx_nghosts * (CONST_NDIM-1) ); // Vy,Vz

    //int siz_to_x2 = siz_x_per_lay * ( fd->fdx_nghosts                         // Vx
    //                                +(fd->fdx_nghosts-1) * (CONST_NDIM-1) ); // Vy,Vz
    //float *sbuff_x1 = mympi->sbuff;
    //float *sbuff_x2 = sbuff_x1 + siz_to_x1;
    //for (int n=0; n < siz_to_x2; n++) {
    //  rbuff[n] = sbuff[n];
    //}

    //for (int n=0; n < mympi->siz_rbuff; n++) {
    //  rbuff[n] = -999.0;
    //}
    blk_stg_el1st_unpack_mesg_vel(fd,mympi,gdinfo,wav,rbuff,mympi->siz_rbuff);

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
    io_snap_nc_put(iosnap, &iosnap_nc, gdinfo, md, wav,
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
sv_eq1st_cart_stg_el_iso_hook(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdry_t    *bdry,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose)
{
  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  float *restrict Txx   = wav->v5d + wav->Txx_pos;
  float *restrict Tyy   = wav->v5d + wav->Tyy_pos;
  float *restrict Tzz   = wav->v5d + wav->Tzz_pos;
  float *restrict Txz   = wav->v5d + wav->Txz_pos;
  float *restrict Tyz   = wav->v5d + wav->Tyz_pos;
  float *restrict Txy   = wav->v5d + wav->Txy_pos;

  float *restrict lam3d = md->lambda;
  float *restrict  mu3d = md->mu;

  float *matVx2Vz = bdry->matVx2Vz2;
  float *matVy2Vz = bdry->matVy2Vz2;

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
  float lam, mu, mu2;
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
  if (bdry->is_sides_free[2][1] == 1)
  {
    sv_eq1st_cart_stg_el_iso_free_vzero(Vx,Vy,Vz,
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
    if (bdry->is_sides_free[2][1] == 1) {
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
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        Vx_ptr = Vx + iptr;
        Vy_ptr = Vy + iptr;
        Vz_ptr = Vz + iptr;

        // Derivatives for Tii
        M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyVy, Vy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // derivatives for Txz
        M_FD_SHIFT_PTR(DxVz, Vz_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzVx, Vx_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // derivatives for Tyz
        M_FD_SHIFT_PTR(DyVz, Vz_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzVy, Vy_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // derivatives for Txy
        M_FD_SHIFT_PTR(DxVy, Vy_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyVx, Vx_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);

        // if surface
        if (k==nk2 && bdry->is_sides_free[CONST_NDIM-1][1]) // at surface, convert
        {
          size_t ij = (i + j * siz_line)*9;
          DzVz = matVx2Vz[ij+3*2+0] * DxVx
               + matVy2Vz[ij+3*2+1] * DyVy;
        }

        // Hooke's equatoin
        Eii = lam * (DxVx + DyVy + DzVz);

        Txx[iptr] += dt * (Eii + mu2 * DxVx);
        Tyy[iptr] += dt * (Eii + mu2 * DyVy);
        Tzz[iptr] += dt * (Eii + mu2 * DzVz);

        Txy[iptr] += dt * mu *( DxVy + DyVx );
        Txz[iptr] += dt * mu *( DxVz + DzVx );
        Tyz[iptr] += dt * mu *( DyVz + DzVy );

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
        lam = lam3d[iptr];
        mu  =  mu3d[iptr];
        mu2 =  mu3d[iptr] * 2.0;

        Vx_ptr = Vx + iptr;
        Vy_ptr = Vy + iptr;
        Vz_ptr = Vz + iptr;

        // Derivatives for Tii
        M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyVy, Vy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // derivatives for Txz
        M_FD_SHIFT_PTR(DxVz, Vz_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzVx, Vx_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // derivatives for Tyz
        M_FD_SHIFT_PTR(DyVz, Vz_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzVy, Vy_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // derivatives for Txy
        M_FD_SHIFT_PTR(DxVy, Vy_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyVx, Vx_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);

        // Hooke's equatoin
        Eii = lam * (DxVx + DyVy + DzVz);

        Txx[iptr] += dt * (Eii + mu2 * DxVx);
        Tyy[iptr] += dt * (Eii + mu2 * DyVy);
        Tzz[iptr] += dt * (Eii + mu2 * DzVz);

        Txy[iptr] += dt * mu *( DxVy + DyVx );
        Txz[iptr] += dt * mu *( DxVz + DzVx );
        Tyz[iptr] += dt * mu *( DyVz + DzVy );

        iptr += 1;
      }
    }
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_eq1st_cart_stg_el_iso_hook_cfspml(fd,Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       lam3d, mu3d, 
                                       dt, dx, dy, dz,
                                       nk2, siz_line,siz_slice,
                                       bdry,
                                       myid, verbose);
  }

  // add source term
  if (src->moment_actived == 1)
  {
    sv_eq1st_cart_stg_el_iso_hook_src(Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                     dt, dx, dy, dz, src,
                                     myid, verbose);
  }
  // end func

  return;
}

void
sv_eq1st_cart_stg_el_iso_moment(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdry_t    *bdry,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose)
{
  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  float *restrict Txx   = wav->v5d + wav->Txx_pos;
  float *restrict Tyy   = wav->v5d + wav->Tyy_pos;
  float *restrict Tzz   = wav->v5d + wav->Tzz_pos;
  float *restrict Txz   = wav->v5d + wav->Txz_pos;
  float *restrict Tyz   = wav->v5d + wav->Tyz_pos;
  float *restrict Txy   = wav->v5d + wav->Txy_pos;

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
  float DxTxx,            DxTxy,DxTxz      ;
  float       DyTyy,      DyTxy,      DyTyz;
  float             DzTzz,      DzTxz,DzTyz;

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
  float *restrict Txx_ptr;
  float *restrict Tyy_ptr;
  float *restrict Tzz_ptr;
  float *restrict Tyz_ptr;
  float *restrict Txz_ptr;
  float *restrict Txy_ptr;

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
  if (bdry->is_sides_free[2][1] == 1)
  {
    sv_eq1st_cart_stg_el_iso_free_simg(Tzz,Txz,Tyz,
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

        Txx_ptr = Txx + iptr;
        Tyy_ptr = Tyy + iptr;
        Tzz_ptr = Tzz + iptr;
        Txz_ptr = Txz + iptr;
        Tyz_ptr = Tyz + iptr;
        Txy_ptr = Txy + iptr;

        // Tii deriv
        M_FD_SHIFT_PTR(DxTxx, Txx_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyTyy, Tyy_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzTzz, Tzz_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

        // Txz deriv
        M_FD_SHIFT_PTR(DxTxz, Txz_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DzTxz, Txz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // Tyz deriv
        M_FD_SHIFT_PTR(DyTyz, Tyz_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
        M_FD_SHIFT_PTR(DzTyz, Tyz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

        // Txy deriv
        M_FD_SHIFT_PTR(DxTxy, Txy_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
        M_FD_SHIFT_PTR(DyTxy, Txy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);

        // Moment equatoin
        Vx[iptr] += slwdt * (DxTxx + DyTxy + DzTxz);
        Vy[iptr] += slwdt * (DxTxy + DyTyy + DzTyz);
        Vz[iptr] += slwdt * (DxTxz + DyTyz + DzTzz);

        iptr += 1;
      }
    }
  }

  // cfs-pml, loop face inside
  if (bdry->is_enable_pml == 1)
  {
    sv_eq1st_cart_stg_el_iso_moment_cfspml(fd,Vx,Vy,Vz,Txx,Tyy,Tzz,Txz,Tyz,Txy,
                                       slw3d, dt, dx, dy, dz,
                                       nk2, siz_line,siz_slice,
                                       bdry,
                                       myid, verbose);
  }

  // add source term
  if (src->force_actived == 1)
  {
    sv_eq1st_cart_stg_el_iso_moment_src(Vx, Vy, Vz,
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
sv_eq1st_cart_stg_el_iso_hook_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict lam3d, float *restrict  mu3d,
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdry_t    *bdry,
    const int myid, const int verbose)
{
  float *matVx2Vz = bdry->matVx2Vz2;
  float *matVy2Vz = bdry->matVy2Vz2;

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
  float lam,mu,lam2mu;

  // local
  int i,j,k;
  int iptr, iptr_j, iptr_k, iptr_a;
  float coef_A, coef_B, coef_D, coef_rB_minus_1;
  float coef_Am, coef_Bm, coef_Dm, coef_rBm_minus_1;

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
      float *restrict ptr_coef_Am = bdry->Am[idim][iside];
      float *restrict ptr_coef_Bm = bdry->Bm[idim][iside];
      float *restrict ptr_coef_Dm = bdry->Dm[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);

      // for each dim
      if (idim == 0 ) // x direction
      {
        float DxVx,DxVy,DxVz;
        // get pml vars, some vars are not used
        float *restrict pml_DxVx   = auxvar->var + auxvar->Vx_pos;
        float *restrict pml_DxVy   = auxvar->var + auxvar->Vy_pos;
        float *restrict pml_DxVz   = auxvar->var + auxvar->Vz_pos;

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
              coef_Dm = ptr_coef_Dm[abs_i];
              coef_Am = ptr_coef_Am[abs_i];
              coef_Bm = ptr_coef_Bm[abs_i];

              coef_rB_minus_1 = coef_B - 1.0;
              coef_rBm_minus_1 = coef_Bm - 1.0;

              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              lam2mu = lam + 2.0 * mu;

              Vx_ptr = Vx + iptr;
              Vy_ptr = Vy + iptr;
              Vz_ptr = Vz + iptr;

              // for x deriv
              M_FD_SHIFT_PTR(DxVx, Vx_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
              M_FD_SHIFT_PTR(DxVy, Vy_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
              M_FD_SHIFT_PTR(DxVz, Vz_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
              float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

              // correct Txx,Tyy,Tzz
              float Vx1 = (2.0 * pml_DxVx[iptr_a] + dt* coef_D *DxVx ) * r_dt_AD;
              float Vx1_corr = coef_rB_minus_1 * DxVx - Vx1 * coef_B;
              Txx[iptr] += dt * lam2mu * Vx1_corr;
              Tyy[iptr] += dt * lam    * Vx1_corr;
              Tzz[iptr] += dt * lam    * Vx1_corr;

              // correct Tij
              float Vz1 =  (2.0 * pml_DxVz[iptr_a] + dt * coef_Dm * DxVz ) * r_dt_ADm; 
              float Vz1_corr = coef_rBm_minus_1 * DxVz - Vz1 * coef_Bm;
              Txz[iptr] += dt * mu * Vz1_corr;

              float Vy1 =  (2.0 * pml_DxVy[iptr_a] + dt * coef_Dm * DxVy ) * r_dt_ADm; 
              float Vy1_corr = coef_rBm_minus_1 * DxVy - Vy1 * coef_Bm;
              Txy[iptr] += dt * mu * Vy1_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DxVx[iptr_a]  = 2.0 * Vx1  - pml_DxVx[iptr_a];
              pml_DxVy[iptr_a]  = 2.0 * Vy1  - pml_DxVy[iptr_a];
              pml_DxVz[iptr_a]  = 2.0 * Vz1  - pml_DxVz[iptr_a];

              // add contributions from free surface condition
              //     only consider Tii on surface
              if (bdry->is_sides_free[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i + j * siz_line)*9;
                float Dx_DzVz = matVx2Vz[ij+3*2+0] * DxVx
                              + matVx2Vz[ij+3*2+1] * DxVy
                              + matVx2Vz[ij+3*2+2] * DxVz;

                // make corr to Hooke's equatoin
                float lamsq_lam2u = lam*lam/lam2mu;
                Txx[iptr] +=   dt * lam * Dx_DzVz * coef_rB_minus_1
                             + dt * lamsq_lam2u * Vx1 * coef_B;
                Tyy[iptr] +=   dt * lam * Dx_DzVz * coef_rB_minus_1
                             + dt * lamsq_lam2u * Vx1 * coef_B;
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
        float DyVx,DyVy,DyVz;
        // get pml vars, some vars are not used
        float *restrict pml_DyVx   = auxvar->var + auxvar->Vx_pos;
        float *restrict pml_DyVy   = auxvar->var + auxvar->Vy_pos;
        float *restrict pml_DyVz   = auxvar->var + auxvar->Vz_pos;

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

            coef_Dm = ptr_coef_Dm[abs_j];
            coef_Am = ptr_coef_Am[abs_j];
            coef_Bm = ptr_coef_Bm[abs_j];
            coef_rBm_minus_1 = coef_Bm - 1.0;

            float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
            float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              lam2mu = lam + 2.0 * mu;

              Vx_ptr = Vx + iptr;
              Vy_ptr = Vy + iptr;
              Vz_ptr = Vz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DyVx, Vx_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);
              M_FD_SHIFT_PTR(DyVy, Vy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
              M_FD_SHIFT_PTR(DyVz, Vz_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);

              // correct Txx,Tyy,Tzz
              float Vy2 =  (2.0 * pml_DyVy[iptr_a] + dt * coef_D * DyVy ) * r_dt_AD;
              float Vy2_corr = coef_rB_minus_1 * DyVy - Vy2 * coef_B;
              Txx[iptr] += dt * lam    * Vy2_corr;
              Tyy[iptr] += dt * lam2mu * Vy2_corr;
              Tzz[iptr] += dt * lam    * Vy2_corr;
              pml_DyVy[iptr_a] = 2.0 * Vy2 - pml_DyVy[iptr_a];

              // correct Tyz
              float Vz2 =  (2.0 * pml_DyVz[iptr_a] + dt * coef_Dm * DyVz ) * r_dt_ADm;
              float Vz2_corr = coef_rBm_minus_1 * DyVz - Vz2 * coef_Bm;
              Tyz[iptr] += dt * mu * Vz2_corr;
              pml_DyVz[iptr_a] = 2.0 * Vz2 - pml_DyVz[iptr_a];

              // correct Txy
              float Vx2 =  (2.0 * pml_DyVx[iptr_a] + dt * coef_Dm * DyVx ) * r_dt_ADm;
              float Vx2_corr = coef_rBm_minus_1 * DyVx - Vx2 * coef_Bm;
              Txy[iptr] += dt * mu * Vx2_corr;
              pml_DyVx[iptr_a] = 2.0 * Vx2 - pml_DyVx[iptr_a];

              // add contributions from free surface condition
              if (bdry->is_sides_free[CONST_NDIM-1][1]==1 && k==nk2)
              {
                // zeta derivatives
                int ij = (i + j * siz_line)*9;
                float Dy_DzVz = matVy2Vz[ij+3*2+0] * DyVx
                              + matVy2Vz[ij+3*2+1] * DyVy
                              + matVy2Vz[ij+3*2+2] * DyVz;

                // make corr to Hooke's equatoin
                float lamsq_lam2u = lam*lam/lam2mu;
                Txx[iptr] +=   dt * lam * Dy_DzVz * coef_rB_minus_1
                             + dt * lamsq_lam2u * Vy2 * coef_B;
                Tyy[iptr] +=   dt * lam * Dy_DzVz * coef_rB_minus_1
                             + dt * lamsq_lam2u * Vy2 * coef_B;
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
        float DzVx,DzVy,DzVz;

        // get pml vars, some vars are not used
        float *restrict pml_DzVx   = auxvar->var + auxvar->Vx_pos;
        float *restrict pml_DzVy   = auxvar->var + auxvar->Vy_pos;
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

          coef_Dm = ptr_coef_Dm[abs_k];
          coef_Am = ptr_coef_Am[abs_k];
          coef_Bm = ptr_coef_Bm[abs_k];
          coef_rBm_minus_1 = coef_Bm - 1.0;

          float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
          float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              lam = lam3d[iptr];
              mu  =  mu3d[iptr];
              lam2mu = lam + 2.0 * mu;

              Vx_ptr = Vx + iptr;
              Vy_ptr = Vy + iptr;
              Vz_ptr = Vz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DzVx, Vx_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);
              M_FD_SHIFT_PTR(DzVy, Vy_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);
              M_FD_SHIFT_PTR(DzVz, Vz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);

              // correct Txx,Tyy,Tzz
              float Vz3 =  (2.0 * pml_DzVz[iptr_a] + dt * coef_D * DzVz ) * r_dt_AD;
              float Vz3_corr = coef_rB_minus_1 * DzVz - Vz3 * coef_B;

              Txx[iptr] += dt *lam   * Vz3_corr;
              Tyy[iptr] += dt *lam   * Vz3_corr;
              Tzz[iptr] += dt *lam2mu* Vz3_corr;

              pml_DzVz[iptr_a] = 2.0 * Vz3 - pml_DzVz[iptr_a];

              // correct Tyz
              float Vy3 =  (2.0 * pml_DzVy[iptr_a] + dt * coef_Dm * DzVy ) * r_dt_ADm;
              float Vy3_corr = coef_rBm_minus_1 * DzVy - Vy3 * coef_Bm;

              Tyz[iptr] += dt * mu * Vy3_corr;

              pml_DzVy[iptr_a] = 2.0 * Vy3 - pml_DzVy[iptr_a];

              // correct Txz
              float Vx3 =  (2.0 * pml_DzVx[iptr_a] + dt * coef_Dm * DzVx ) * r_dt_ADm;
              float Vx3_corr = coef_rBm_minus_1 * DzVx - Vx3 * coef_Bm;

              Txz[iptr] += dt * mu * Vx3_corr;

              pml_DzVx[iptr_a] = 2.0 * Vx3 - pml_DzVx[iptr_a];

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
sv_eq1st_cart_stg_el_iso_moment_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict slw3d,
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdry_t *bdry,
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
  float coef_Am, coef_Bm, coef_Dm, coef_rBm_minus_1;
  float slw;

  float *restrict Txx_ptr;
  float *restrict Tyy_ptr;
  float *restrict Tzz_ptr;
  float *restrict Tyz_ptr;
  float *restrict Txz_ptr;
  float *restrict Txy_ptr;

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
      float *restrict ptr_coef_Am = bdry->Am[idim][iside];
      float *restrict ptr_coef_Bm = bdry->Bm[idim][iside];
      float *restrict ptr_coef_Dm = bdry->Dm[idim][iside];

      bdrypml_auxvar_t *auxvar = &(bdry->auxvar[idim][iside]);

      // for each dim
      if (idim == 0 ) // x direction
      {
        float DxTxx, DxTxy, DxTxz;

        // get pml vars, some vars are not used
        float *restrict pml_DxTxx   = auxvar->var + auxvar->Txx_pos;
        float *restrict pml_DxTxy   = auxvar->var + auxvar->Txy_pos;
        float *restrict pml_DxTxz   = auxvar->var + auxvar->Txz_pos;

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

              coef_Dm = ptr_coef_Dm[abs_i];
              coef_Am = ptr_coef_Am[abs_i];
              coef_Bm = ptr_coef_Bm[abs_i];
              coef_rBm_minus_1 = coef_Bm - 1.0;

              float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
              float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

              // medium
              slw = slw3d[iptr];

              Txx_ptr = Txx + iptr;
              Txz_ptr = Txz + iptr;
              Txy_ptr = Txy + iptr;

              //  x deriv
              M_FD_SHIFT_PTR(DxTxx, Txx_ptr, lfdx_len, lfdx_shift_F, lfdx_coef, n_fd);
              M_FD_SHIFT_PTR(DxTxy, Txy_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);
              M_FD_SHIFT_PTR(DxTxz, Txz_ptr, lfdx_len, lfdx_shift_B, lfdx_coef, n_fd);

              // correct Vx
              float Txx1 = (2.0 * pml_DxTxx[iptr_a] + dt* coef_Dm *DxTxx ) * r_dt_ADm;
              float Txx1_corr = coef_rBm_minus_1 * DxTxx - Txx1 * coef_Bm;
              Vx[iptr] += dt * slw * Txx1_corr;

              // correct Vy
              float Txy1 = (2.0 * pml_DxTxy[iptr_a] + dt* coef_D *DxTxy ) * r_dt_AD;
              float Txy1_corr = coef_rB_minus_1 * DxTxy - Txy1 * coef_B;
              Vy[iptr] += dt * slw * Txy1_corr;

              // correct Vz
              float Txz1 = (2.0 * pml_DxTxz[iptr_a] + dt* coef_D *DxTxz ) * r_dt_AD;
              float Txz1_corr = coef_rB_minus_1 * DxTxz - Txz1 * coef_B;
              Vz[iptr] += dt * slw * Txz1_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DxTxx[iptr_a]  = 2.0 * Txx1  - pml_DxTxx[iptr_a];
              pml_DxTxy[iptr_a]  = 2.0 * Txy1  - pml_DxTxy[iptr_a];
              pml_DxTxz[iptr_a]  = 2.0 * Txz1  - pml_DxTxz[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else if (idim == 1) // y direction
      {
        float DyTxy, DyTyy, DyTyz;

        // get pml vars, some vars are not used
        float *restrict pml_DyTxy   = auxvar->var + auxvar->Txy_pos;
        float *restrict pml_DyTyy   = auxvar->var + auxvar->Tyy_pos;
        float *restrict pml_DyTyz   = auxvar->var + auxvar->Tyz_pos;

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

            coef_Dm = ptr_coef_Dm[abs_j];
            coef_Am = ptr_coef_Am[abs_j];
            coef_Bm = ptr_coef_Bm[abs_j];
            coef_rBm_minus_1 = coef_Bm - 1.0;

            float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
            float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              slw = slw3d[iptr];

              Txy_ptr = Txy + iptr;
              Tyy_ptr = Tyy + iptr;
              Tyz_ptr = Tyz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DyTxy, Txy_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);
              M_FD_SHIFT_PTR(DyTyy, Tyy_ptr, lfdy_len, lfdy_shift_F, lfdy_coef, n_fd);
              M_FD_SHIFT_PTR(DyTyz, Tyz_ptr, lfdy_len, lfdy_shift_B, lfdy_coef, n_fd);

              // correct Vx
              float Txy2 = (2.0 * pml_DyTxy[iptr_a] + dt* coef_D *DyTxy ) * r_dt_AD;
              float Txy2_corr = coef_rB_minus_1 * DyTxy - Txy2 * coef_B;
              Vx[iptr] += dt * slw * Txy2_corr;

              // correct Vy
              float Tyy2 = (2.0 * pml_DyTyy[iptr_a] + dt* coef_Dm *DyTyy ) * r_dt_ADm;
              float Tyy2_corr = coef_rBm_minus_1 * DyTyy - Tyy2 * coef_Bm;
              Vy[iptr] += dt * slw * Tyy2_corr;

              // correct Vz
              float Tyz2 = (2.0 * pml_DyTyz[iptr_a] + dt* coef_D *DyTyz ) * r_dt_AD;
              float Tyz2_corr = coef_rB_minus_1 * DyTyz - Tyz2 * coef_B;
              Vz[iptr] += dt * slw * Tyz2_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DyTxy[iptr_a]  = 2.0 * Txy2  - pml_DyTxy[iptr_a];
              pml_DyTyy[iptr_a]  = 2.0 * Tyy2  - pml_DyTyy[iptr_a];
              pml_DyTyz[iptr_a]  = 2.0 * Tyz2  - pml_DyTyz[iptr_a];

              // incr index
              iptr   += 1;
              iptr_a += 1;
            } // i
          } // j
        } // k
      }
      else // z direction
      {
        float DzTxz, DzTyz, DzTzz;

        // get pml vars, some vars are not used
        float *restrict pml_DzTxz   = auxvar->var + auxvar->Txz_pos;
        float *restrict pml_DzTyz   = auxvar->var + auxvar->Tyz_pos;
        float *restrict pml_DzTzz   = auxvar->var + auxvar->Tzz_pos;

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

          coef_Dm = ptr_coef_Dm[abs_k];
          coef_Am = ptr_coef_Am[abs_k];
          coef_Bm = ptr_coef_Bm[abs_k];
          coef_rBm_minus_1 = coef_Bm - 1.0;

          float r_dt_AD = 1.0/(2.0 + dt * ( coef_A + coef_D)) ;
          float r_dt_ADm = 1.0/(2.0 + dt * ( coef_Am + coef_Dm)) ;

          for (j=abs_nj1; j<=abs_nj2; j++)
          {
            iptr_j = iptr_k + j * siz_line;
            iptr = iptr_j + abs_ni1;
            for (i=abs_ni1; i<=abs_ni2; i++)
            {
              // medium
              slw = slw3d[iptr];

              Txz_ptr = Txz + iptr;
              Tyz_ptr = Tyz + iptr;
              Tzz_ptr = Tzz + iptr;

              // derivatives
              M_FD_SHIFT_PTR(DzTxz, Txz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);
              M_FD_SHIFT_PTR(DzTyz, Tyz_ptr, lfdz_len, lfdz_shift_B, lfdz_coef, n_fd);
              M_FD_SHIFT_PTR(DzTzz, Tzz_ptr, lfdz_len, lfdz_shift_F, lfdz_coef, n_fd);

              // correct Vx
              float Txz3 = (2.0 * pml_DzTxz[iptr_a] + dt* coef_D *DzTxz ) * r_dt_AD;
              float Txz3_corr = coef_rB_minus_1 * DzTxz - Txz3 * coef_B;
              Vx[iptr] += dt * slw * Txz3_corr;

              // correct Vy
              float Tyz3 = (2.0 * pml_DzTyz[iptr_a] + dt* coef_D *DzTyz ) * r_dt_AD;
              float Tyz3_corr = coef_rB_minus_1 * DzTyz - Tyz3 * coef_B;
              Vy[iptr] += dt * slw * Tyz3_corr;

              // correct Vz
              float Tzz3 = (2.0 * pml_DzTzz[iptr_a] + dt* coef_Dm *DzTzz ) * r_dt_ADm;
              float Tzz3_corr = coef_rBm_minus_1 * DzTzz - Tzz3 * coef_Bm;
              Vz[iptr] += dt * slw * Tzz3_corr;

              // updat aux var
              //   a1 = alpha + d / beta, dealt in abs_set_cfspml
              pml_DzTxz[iptr_a]  = 2.0 * Txz3  - pml_DzTxz[iptr_a];
              pml_DzTyz[iptr_a]  = 2.0 * Tyz3  - pml_DzTyz[iptr_a];
              pml_DzTzz[iptr_a]  = 2.0 * Tzz3  - pml_DzTzz[iptr_a];

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
sv_eq1st_cart_stg_el_iso_hook_src(
    float *restrict Txx, float *restrict Tyy, float *restrict Tzz,
    float *restrict Txz, float *restrict Tyz, float *restrict Txy,
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

      Mxx = src->Mxx[iptr_cur_stage];
      Myy = src->Myy[iptr_cur_stage];
      Mzz = src->Mzz[iptr_cur_stage];
      Mxz = src->Mxz[iptr_cur_stage];
      Myz = src->Myz[iptr_cur_stage];
      Mxy = src->Mxy[iptr_cur_stage];
      
      // for extend points
      for (int i_ext=0; i_ext < src->ext_num[is]; i_ext++)
      {
        int   iptr = ptr_ext_indx[i_ext];
        float coef = ptr_ext_coef[i_ext];
        float rjac = dt * coef / vol;

        Txx[iptr] -= Mxx * rjac;
        Tyy[iptr] -= Myy * rjac;
        Tzz[iptr] -= Mzz * rjac;
        Txz[iptr] -= Mxz * rjac;
        Tyz[iptr] -= Myz * rjac;
        Txy[iptr] -= Mxy * rjac;
      } // i_ext

    } // it
  } // is

  return ierr;
}

int
sv_eq1st_cart_stg_el_iso_moment_src(
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
 * free surface coef
 ******************************************************************************/

int
sv_eq1st_cart_stg_el_iso_dvh2dvz(gdinfo_t   *gdinfo,
                                 md_t       *md,
                                 bdry_t      *bdryfree,
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
 * free surface boundary
 ******************************************************************************/

/*
 * implement traction image boundary 
 */

void
sv_eq1st_cart_stg_el_iso_free_simg(
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

      Tyz[iptr] = -Tyz[iptr - siz_slice];
      Txz[iptr] = -Txz[iptr - siz_slice];

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
        Tyz[iptr_gho] = -Tyz[iptr_phy - siz_slice];
        Txz[iptr_gho] = -Txz[iptr_phy - siz_slice];
      }
    }
  }

  return;
}

/*
 * set zero velocity above free surface
 */

void
sv_eq1st_cart_stg_el_iso_free_vzero(
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
      Vz[iptr] = 0.0;

      // next
      iptr += 1;
    }
  }

  // above surface
  for (size_t k=nk2+1; k<nz; k++)
  {
    for (size_t j=nj1; j<=nj2; j++)
    {
      for (size_t i=ni1; i<=ni2; i++)
      {
        size_t iptr = i + j * siz_line + k     * siz_slice;

        Vz[iptr] = 0.0;
        Vx[iptr] = 0.0;
        Vy[iptr] = 0.0;
      }
    }
  }

  return;
}
