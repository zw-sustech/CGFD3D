/*
********************************************************************************
* Curve grid metric calculation using MacCormack scheme                        *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "fd_t.h"
#include "gd_t.h"
#include "isPointInHexahedron.h"

// used in read grid file
#define M_gd_INDEX( i, j, k, ni, nj ) ( ( i ) + ( j ) * ( ni ) + ( k ) * ( ni ) * ( nj ) )

void 
gd_curv_init(gdinfo_t *gdinfo, gd_t *gdcurv)
{
  /*
   * 0-2: x3d, y3d, z3d
   */

  gdcurv->type = GD_TYPE_CURV;

  gdcurv->nx   = gdinfo->nx;
  gdcurv->ny   = gdinfo->ny;
  gdcurv->nz   = gdinfo->nz;
  gdcurv->ncmp = CONST_NDIM;

  gdcurv->siz_iy   = gdcurv->nx;
  gdcurv->siz_iz   = gdcurv->nx * gdcurv->ny;
  gdcurv->siz_icmp = gdcurv->nx * gdcurv->ny * gdcurv->nz;
  
  // vars
  gdcurv->v4d = (float *) fdlib_mem_calloc_1d_float(
                  gdcurv->siz_icmp * gdcurv->ncmp, 0.0, "gd_curv_init");
  if (gdcurv->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // position of each v4d
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(gdcurv->ncmp,
                                                         0,
                                                         "gd_curv_init");
  
  // name of each v4d
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(gdcurv->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "gd_curv_init");
  
  // set value
  int icmp = 0;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","x");
  gdcurv->x3d = gdcurv->v4d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","y");
  gdcurv->y3d = gdcurv->v4d + cmp_pos[icmp];

  icmp += 1;
  cmp_pos[icmp] = icmp * gdcurv->siz_icmp;
  sprintf(cmp_name[icmp],"%s","z");
  gdcurv->z3d = gdcurv->v4d + cmp_pos[icmp];
  
  // set pointer
  gdcurv->cmp_pos  = cmp_pos;
  gdcurv->cmp_name = cmp_name;

  return;
}

void 
gd_curv_metric_init(gdinfo_t        *gdinfo,
                    gdcurv_metric_t *metric)
{
  const int num_grid_vars = 10;
  /*
   * 0: jac
   * 1-3: xi_x, xi_y, xi_z
   * 4-6: eta_x, eta_y, eta_z
   * 7-9: zeta_x, zeta_y, zeta_z
   */

  metric->nx   = gdinfo->nx;
  metric->ny   = gdinfo->ny;
  metric->nz   = gdinfo->nz;
  metric->ncmp = num_grid_vars;

  metric->siz_iy   = metric->nx;
  metric->siz_iz   = metric->nx * metric->ny;
  metric->siz_icmp = metric->nx * metric->ny * metric->nz;
  
  // vars
  metric->v4d = (float *) fdlib_mem_calloc_1d_float(
                  metric->siz_icmp * metric->ncmp, 0.0, "gd_curv_init_g4d");
  if (metric->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc metric vars\n");
      fflush(stderr);
  }
  
  // position of each v4d
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(metric->ncmp,
                                                         0, 
                                                         "gd_curv_metric_init");
  
  // name of each v4d
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(metric->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "gd_curv_metric_init");
  
  // set value
  for (int icmp=0; icmp < metric->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * metric->siz_icmp;
  }

  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","jac");
  metric->jac = metric->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","xi_x");
  metric->xi_x = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","xi_y");
  metric->xi_y = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","xi_z");
  metric->xi_z = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","eta_x");
  metric->eta_x = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","eta_y");
  metric->eta_y = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","eta_z");
  metric->eta_z = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","zeta_x");
  metric->zeta_x = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","zeta_y");
  metric->zeta_y = metric->v4d + cmp_pos[icmp];
  
  icmp += 1;
  sprintf(cmp_name[icmp],"%s","zeta_z");
  metric->zeta_z = metric->v4d + cmp_pos[icmp];
  
  // set pointer
  metric->cmp_pos  = cmp_pos;
  metric->cmp_name = cmp_name;
}

//
// need to change to use fdlib_math.c
//
void
gd_curv_metric_cal(gdinfo_t        *gdinfo,
                   gd_t        *gdcurv,
                   gdcurv_metric_t *metric,
                   int fd_len, int *restrict fd_indx, float *restrict fd_coef)
{
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
  float *restrict x3d  = gdcurv->x3d;
  float *restrict y3d  = gdcurv->y3d;
  float *restrict z3d  = gdcurv->z3d;
  float *restrict jac3d= metric->jac;
  float *restrict xi_x = metric->xi_x;
  float *restrict xi_y = metric->xi_y;
  float *restrict xi_z = metric->xi_z;
  float *restrict et_x = metric->eta_x;
  float *restrict et_y = metric->eta_y;
  float *restrict et_z = metric->eta_z;
  float *restrict zt_x = metric->zeta_x;
  float *restrict zt_y = metric->zeta_y;
  float *restrict zt_z = metric->zeta_z;

  float x_xi, x_et, x_zt;
  float y_xi, y_et, y_zt;
  float z_xi, z_et, z_zt;
  float jac;
  float vec1[3], vec2[3], vec3[3], vecg[3];
  int n_fd;

  // use local stack array for speedup
  float  lfd_coef [fd_len];
  int    lfdx_shift[fd_len];
  int    lfdy_shift[fd_len];
  int    lfdz_shift[fd_len];
  // put fd op into local array
  for (int k=0; k < fd_len; k++) {
    lfd_coef [k] = fd_coef[k];
    lfdx_shift[k] = fd_indx[k]            ;
    lfdy_shift[k] = fd_indx[k] * siz_line ;
    lfdz_shift[k] = fd_indx[k] * siz_slice;
  }

  for (size_t k = nk1; k <= nk2; k++){
    for (size_t j = nj1; j <= nj2; j++) {
      for (size_t i = ni1; i <= ni2; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;

        x_xi = 0.0; x_et = 0.0; x_zt = 0.0;
        y_xi = 0.0; y_et = 0.0; y_zt = 0.0;
        z_xi = 0.0; z_et = 0.0; z_zt = 0.0;

        M_FD_SHIFT(x_xi, x3d, iptr, fd_len, lfdx_shift, lfd_coef, n_fd);
        M_FD_SHIFT(y_xi, y3d, iptr, fd_len, lfdx_shift, lfd_coef, n_fd);
        M_FD_SHIFT(z_xi, z3d, iptr, fd_len, lfdx_shift, lfd_coef, n_fd);

        M_FD_SHIFT(x_et, x3d, iptr, fd_len, lfdy_shift, lfd_coef, n_fd);
        M_FD_SHIFT(y_et, y3d, iptr, fd_len, lfdy_shift, lfd_coef, n_fd);
        M_FD_SHIFT(z_et, z3d, iptr, fd_len, lfdy_shift, lfd_coef, n_fd);

        M_FD_SHIFT(x_zt, x3d, iptr, fd_len, lfdz_shift, lfd_coef, n_fd);
        M_FD_SHIFT(y_zt, y3d, iptr, fd_len, lfdz_shift, lfd_coef, n_fd);
        M_FD_SHIFT(z_zt, z3d, iptr, fd_len, lfdz_shift, lfd_coef, n_fd);

        vec1[0] = x_xi; vec1[1] = y_xi; vec1[2] = z_xi;
        vec2[0] = x_et; vec2[1] = y_et; vec2[2] = z_et;
        vec3[0] = x_zt; vec3[1] = y_zt; vec3[2] = z_zt;

        fdlib_math_cross_product(vec1, vec2, vecg);
        jac = fdlib_math_dot_product(vecg, vec3);
        jac3d[iptr]  = jac;

        fdlib_math_cross_product(vec2, vec3, vecg);
        xi_x[iptr] = vecg[0] / jac;
        xi_y[iptr] = vecg[1] / jac;
        xi_z[iptr] = vecg[2] / jac;

        fdlib_math_cross_product(vec3, vec1, vecg);
        et_x[iptr] = vecg[0] / jac;
        et_y[iptr] = vecg[1] / jac;
        et_z[iptr] = vecg[2] / jac;

        fdlib_math_cross_product(vec1, vec2, vecg);
        zt_x[iptr] = vecg[0] / jac;
        zt_y[iptr] = vecg[1] / jac;
        zt_z[iptr] = vecg[2] / jac;
      }
    }
  }
    
  // extend to ghosts. may replaced by mpi exchange
  // x1, mirror
  for (size_t k = 0; k < nz; k++){
    for (size_t j = 0; j < ny; j++) {
      for (size_t i = 0; i < ni1; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = jac3d[iptr + (ni1-i)*2 -1 ];
         xi_x[iptr] =  xi_x[iptr + (ni1-i)*2 -1 ];
         xi_y[iptr] =  xi_y[iptr + (ni1-i)*2 -1 ];
         xi_z[iptr] =  xi_z[iptr + (ni1-i)*2 -1 ];
         et_x[iptr] =  et_x[iptr + (ni1-i)*2 -1 ];
         et_y[iptr] =  et_y[iptr + (ni1-i)*2 -1 ];
         et_z[iptr] =  et_z[iptr + (ni1-i)*2 -1 ];
         zt_x[iptr] =  zt_x[iptr + (ni1-i)*2 -1 ];
         zt_y[iptr] =  zt_y[iptr + (ni1-i)*2 -1 ];
         zt_z[iptr] =  zt_z[iptr + (ni1-i)*2 -1 ];
      }
    }
  }
  // x2, mirror
  for (size_t k = 0; k < nz; k++){
    for (size_t j = 0; j < ny; j++) {
      for (size_t i = ni2+1; i < nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = jac3d[iptr - (i-ni2)*2 +1 ];
         xi_x[iptr] =  xi_x[iptr - (i-ni2)*2 +1 ];
         xi_y[iptr] =  xi_y[iptr - (i-ni2)*2 +1 ];
         xi_z[iptr] =  xi_z[iptr - (i-ni2)*2 +1 ];
         et_x[iptr] =  et_x[iptr - (i-ni2)*2 +1 ];
         et_y[iptr] =  et_y[iptr - (i-ni2)*2 +1 ];
         et_z[iptr] =  et_z[iptr - (i-ni2)*2 +1 ];
         zt_x[iptr] =  zt_x[iptr - (i-ni2)*2 +1 ];
         zt_y[iptr] =  zt_y[iptr - (i-ni2)*2 +1 ];
         zt_z[iptr] =  zt_z[iptr - (i-ni2)*2 +1 ];
      }
    }
  }
  // y1, mirror
  for (size_t k = 0; k < nz; k++){
    for (size_t j = 0; j < nj1; j++) {
      for (size_t i = 0; i < nx; i++) {
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = jac3d[iptr + ((nj1-j)*2 -1) * siz_line ];
         xi_x[iptr] =  xi_x[iptr + ((nj1-j)*2 -1) * siz_line ];
         xi_y[iptr] =  xi_y[iptr + ((nj1-j)*2 -1) * siz_line ];
         xi_z[iptr] =  xi_z[iptr + ((nj1-j)*2 -1) * siz_line ];
         et_x[iptr] =  et_x[iptr + ((nj1-j)*2 -1) * siz_line ];
         et_y[iptr] =  et_y[iptr + ((nj1-j)*2 -1) * siz_line ];
         et_z[iptr] =  et_z[iptr + ((nj1-j)*2 -1) * siz_line ];
         zt_x[iptr] =  zt_x[iptr + ((nj1-j)*2 -1) * siz_line ];
         zt_y[iptr] =  zt_y[iptr + ((nj1-j)*2 -1) * siz_line ];
         zt_z[iptr] =  zt_z[iptr + ((nj1-j)*2 -1) * siz_line ];
      }
    }
  }
  // y2, mirror
  for (size_t k = 0; k < nz; k++){
    for (size_t j = nj2+1; j < ny; j++) {
      for (size_t i = 0; i < nx; i++) {
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = jac3d[iptr - ((j-nj2)*2 -1) * siz_line ];
         xi_x[iptr] =  xi_x[iptr - ((j-nj2)*2 -1) * siz_line ];
         xi_y[iptr] =  xi_y[iptr - ((j-nj2)*2 -1) * siz_line ];
         xi_z[iptr] =  xi_z[iptr - ((j-nj2)*2 -1) * siz_line ];
         et_x[iptr] =  et_x[iptr - ((j-nj2)*2 -1) * siz_line ];
         et_y[iptr] =  et_y[iptr - ((j-nj2)*2 -1) * siz_line ];
         et_z[iptr] =  et_z[iptr - ((j-nj2)*2 -1) * siz_line ];
         zt_x[iptr] =  zt_x[iptr - ((j-nj2)*2 -1) * siz_line ];
         zt_y[iptr] =  zt_y[iptr - ((j-nj2)*2 -1) * siz_line ];
         zt_z[iptr] =  zt_z[iptr - ((j-nj2)*2 -1) * siz_line ];
      }
    }
  }
  // z1, mirror
  for (size_t k = 0; k < nk1; k++) {
    for (size_t j = 0; j < ny; j++){
      for (size_t i = 0; i < nx; i++) {
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = jac3d[iptr + ((nk1-k)*2 -1) * siz_slice ];
         xi_x[iptr] =  xi_x[iptr + ((nk1-k)*2 -1) * siz_slice ];
         xi_y[iptr] =  xi_y[iptr + ((nk1-k)*2 -1) * siz_slice ];
         xi_z[iptr] =  xi_z[iptr + ((nk1-k)*2 -1) * siz_slice ];
         et_x[iptr] =  et_x[iptr + ((nk1-k)*2 -1) * siz_slice ];
         et_y[iptr] =  et_y[iptr + ((nk1-k)*2 -1) * siz_slice ];
         et_z[iptr] =  et_z[iptr + ((nk1-k)*2 -1) * siz_slice ];
         zt_x[iptr] =  zt_x[iptr + ((nk1-k)*2 -1) * siz_slice ];
         zt_y[iptr] =  zt_y[iptr + ((nk1-k)*2 -1) * siz_slice ];
         zt_z[iptr] =  zt_z[iptr + ((nk1-k)*2 -1) * siz_slice ];
      }
    }
  }
  // z2, mirror
  for (size_t k = nk2+1; k < nz; k++) {
    for (size_t j = 0; j < ny; j++){
      for (size_t i = 0; i < nx; i++) {
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = jac3d[iptr - ((k-nk2)*2 -1) * siz_slice ];
         xi_x[iptr] =  xi_x[iptr - ((k-nk2)*2 -1) * siz_slice ];
         xi_y[iptr] =  xi_y[iptr - ((k-nk2)*2 -1) * siz_slice ];
         xi_z[iptr] =  xi_z[iptr - ((k-nk2)*2 -1) * siz_slice ];
         et_x[iptr] =  et_x[iptr - ((k-nk2)*2 -1) * siz_slice ];
         et_y[iptr] =  et_y[iptr - ((k-nk2)*2 -1) * siz_slice ];
         et_z[iptr] =  et_z[iptr - ((k-nk2)*2 -1) * siz_slice ];
         zt_x[iptr] =  zt_x[iptr - ((k-nk2)*2 -1) * siz_slice ];
         zt_y[iptr] =  zt_y[iptr - ((k-nk2)*2 -1) * siz_slice ];
         zt_z[iptr] =  zt_z[iptr - ((k-nk2)*2 -1) * siz_slice ];
      }
    }
  }
}

//
// exchange metics
//
void
gd_curv_metric_exchange(gdinfo_t        *gdinfo,
                        gdcurv_metric_t *metric,
                        int             *neighid,
                        MPI_Comm        topocomm)
{
  int nx  = gdinfo->nx;
  int ny  = gdinfo->ny;
  int nz  = gdinfo->nz;
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;

  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;
  size_t siz_volume = gdinfo->siz_icmp;
  float *restrict g3d = metric->v4d;

  // extend to ghosts, using mpi exchange
  // NOTE in different myid, nx(or ny) may not equal
  // so send type DTypeXL not equal recv type DTypeXL
  int s_iptr;
  int r_iptr;

  MPI_Status status;
  MPI_Datatype DTypeXL, DTypeYL;

  MPI_Type_vector(ny*nz*metric->ncmp,
                  3,
                  nx,
                  MPI_FLOAT,
                  &DTypeXL);
  MPI_Type_vector(nz*metric->ncmp,
                  3*nx,
                  nx*ny,
                  MPI_FLOAT,
                  &DTypeYL);
  MPI_Type_commit(&DTypeXL);
  MPI_Type_commit(&DTypeYL);
  // to X1
  s_iptr = ni1 + 0 * siz_line + 0 * siz_slice;        //sendbuff point (ni1,ny1,nz1)
  r_iptr = (ni2+1) + 0 * siz_line + 0 * siz_slice;    //recvbuff point (ni2+1,ny1,nz1)
  MPI_Sendrecv(&g3d[s_iptr],1,DTypeXL,neighid[0],110,
               &g3d[r_iptr],1,DTypeXL,neighid[1],110,
               topocomm,&status);
  // to X2
  s_iptr = (ni2-3+1) + 0 * siz_line + 0 * siz_slice;    //sendbuff point (ni2-3+1,ny1,nz1)
  r_iptr = (ni1-3) + 0 * siz_line + 0 * siz_slice;      //recvbuff point (ni1-3,ny1,nz1)
  MPI_Sendrecv(&g3d[s_iptr],1,DTypeXL,neighid[1],120,
               &g3d[r_iptr],1,DTypeXL,neighid[0],120,
               topocomm,&status);
  // to Y1
  s_iptr = 0 + nj1 * siz_line + 0 * siz_slice;        //sendbuff point (nx1,nj1,nz1)
  r_iptr = 0 + (nj2+1) * siz_line + 0 * siz_slice;    //recvbuff point (nx1,nj2+1,nz1)
  MPI_Sendrecv(&g3d[s_iptr],1,DTypeYL,neighid[2],210,
               &g3d[r_iptr],1,DTypeYL,neighid[3],210,
               topocomm,&status);
  // to Y2
  s_iptr = 0 + (nj2-3+1) * siz_line + 0 * siz_slice;   //sendbuff point (nx1,nj2-3+1,nz1)
  r_iptr = 0 + (nj1-3) * siz_line + 0 * siz_slice;     //recvbuff point (nx1,nj1-3,nz1)
  MPI_Sendrecv(&g3d[s_iptr],1,DTypeYL,neighid[3],220,
               &g3d[r_iptr],1,DTypeYL,neighid[2],220,
               topocomm,&status);
/*
  // for test.  check send_recv whether success
  for (size_t i = ni1; i < ni1+3; i++)
  {
    int j = nj1;
    int k = nk1;
    size_t iptr = i + j * siz_line + k * siz_slice;
    if(myid2[0]==1 && myid2[1]==0 ) fprintf(stdout,"aaa  %f %f %f aaa\n",jac3d[iptr],xi_x[iptr], zt_y[iptr]);
  }
  for (size_t i = ni2+1; i <= ni2+3; i++)
  {
    int j = nj1;
    int k = nk1;
    size_t iptr = i + j * siz_line + k * siz_slice;
    if(myid2[0]==0 && myid2[1]==0) fprintf(stdout,"**a %f %f %f **a\n",jac3d[iptr],xi_x[iptr], zt_y[iptr]);
  }

  for (size_t i = ni2-2; i <= ni2; i++)
  {
    int j = nj1;
    int k = nk1;
    size_t iptr = i + j * siz_line + k * siz_slice;
    if(myid2[0]==0 && myid2[1]==1) fprintf(stdout,"bbb %f %f %f bbb\n",jac3d[iptr],xi_x[iptr], zt_y[iptr]);
  }
  for (size_t i = ni1-3; i < ni1; i++)
  {
    int j = nj1;
    int k = nk1;
    size_t iptr = i + j * siz_line + k * siz_slice;
    if(myid2[0]==1 && myid2[1]==1) fprintf(stdout,"**b %f %f %f **b\n",jac3d[iptr],xi_x[iptr], zt_y[iptr]);
  }

  for (size_t j = nj2-2; j <= nj2; j++)
  {
    int i = ni1+1;
    int k = nk1+1;
    size_t iptr = i + j * siz_line + k * siz_slice;
    if(myid2[0]==0 && myid2[1]==0) fprintf(stdout,"ccc %f %f %f ccc\n",jac3d[iptr],xi_x[iptr], zt_y[iptr]);
  }
  for (size_t j = nj1-3; j < nj1; j++)
  {
    int i = ni1+1;
    int k = nk1+1;
    size_t iptr = i + j * siz_line + k * siz_slice;
    if(myid2[0]==0 && myid2[1]==1) fprintf(stdout,"**c %f %f %f **c\n",jac3d[iptr],xi_x[iptr], zt_y[iptr]);
  }
*/

  // for test
  /*
  for (size_t k = 0; k < nz; k++){
    for (size_t j = 0; j < ny; j++) {
      for (size_t i = 0; i < nx; i++)
      {
        float dh=100.0;
        size_t iptr = i + j * siz_line + k * siz_slice;
        jac3d[iptr] = dh*dh*dh;
         xi_x[iptr] = 1.0/dh;
         xi_y[iptr] =  0.0;
         xi_z[iptr] =  0.0;
         et_x[iptr] =  0.0;
         et_y[iptr] = 1.0/dh;
         et_z[iptr] =  0.0;
         zt_x[iptr] =  0.0;
         zt_y[iptr] =  0.0;
         zt_z[iptr] = 1.0/dh;
      }
    }
  }
  */
}

/*
 * generate cartesian grid for curv struct
 */
void
gd_curv_gen_cart(
  gdinfo_t *gdinfo,
  gd_t *gdcurv,
  float dx, float x0_glob,
  float dy, float y0_glob,
  float dz, float z0_glob)
{
  float *x3d = gdcurv->x3d;
  float *y3d = gdcurv->y3d;
  float *z3d = gdcurv->z3d;

  float x0 = x0_glob + (gdinfo->ni1_to_glob_phys0 - gdinfo->fdx_nghosts) * dx;
  float y0 = y0_glob + (gdinfo->nj1_to_glob_phys0 - gdinfo->fdy_nghosts) * dy;
  float z0 = z0_glob + (gdinfo->nk1_to_glob_phys0 - gdinfo->fdz_nghosts) * dz;

  size_t iptr = 0;
  for (size_t k=0; k<gdcurv->nz; k++)
  {
    for (size_t j=0; j<gdcurv->ny; j++)
    {
      for (size_t i=0; i<gdcurv->nx; i++)
      {
        x3d[iptr] = x0 + i * dx;
        y3d[iptr] = y0 + j * dy;
        z3d[iptr] = z0 + k * dz;

        iptr++;
      }
    }
  }

  return;
}

/*
 * generate cartesian grid
 */

void 
gd_cart_init_set(gdinfo_t *gdinfo, gd_t *gdcart,
  float dx, float x0_glob,
  float dy, float y0_glob,
  float dz, float z0_glob)
{
  /*
   * 0-2: x3d, y3d, z3d
   */

  gdcart->type = GD_TYPE_CART;

  gdcart->nx   = gdinfo->nx;
  gdcart->ny   = gdinfo->ny;
  gdcart->nz   = gdinfo->nz;
  gdcart->ncmp = CONST_NDIM;

  gdcart->siz_iy   = gdcart->nx;
  gdcart->siz_iz   = gdcart->nx * gdcart->ny;
  gdcart->siz_icmp = gdcart->nx * gdcart->ny * gdcart->nz;
  
  // vars
  float *x1d = (float *) fdlib_mem_calloc_1d_float(
                  gdcart->nx, 0.0, "gd_cart_init");
  float *y1d = (float *) fdlib_mem_calloc_1d_float(
                  gdcart->ny, 0.0, "gd_cart_init");
  float *z1d = (float *) fdlib_mem_calloc_1d_float(
                  gdcart->nz, 0.0, "gd_cart_init");
  if (z1d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }

  float x0 = x0_glob + (gdinfo->ni1_to_glob_phys0 - gdinfo->fdx_nghosts) * dx;
  float y0 = y0_glob + (gdinfo->nj1_to_glob_phys0 - gdinfo->fdy_nghosts) * dy;
  float z0 = z0_glob + (gdinfo->nk1_to_glob_phys0 - gdinfo->fdz_nghosts) * dz;

  for (size_t k=0; k< gdcart->nz; k++)
  {
        z1d[k] = z0 + k * dz;
  }
  for (size_t j=0; j< gdcart->ny; j++)
  {
        y1d[j] = y0 + j * dy;
  }
  for (size_t i=0; i< gdcart->nx; i++)
  {
        x1d[i] = x0 + i * dx;
  }

  gdcart->dx = dx;
  gdcart->dy = dy;
  gdcart->dz = dz;

  gdcart->xmin = x0;
  gdcart->ymin = y0;
  gdcart->zmin = z0;
  gdcart->xmax = x0 + (gdcart->nx-1) * dx;
  gdcart->ymax = y0 + (gdcart->ny-1) * dy;
  gdcart->zmax = z0 + (gdcart->nz-1) * dz;

  gdcart->x0_glob = x0_glob;
  gdcart->y0_glob = y0_glob;
  gdcart->z0_glob = z0_glob;

  gdcart->x1d = x1d;
  gdcart->y1d = y1d;
  gdcart->z1d = z1d;
  
  return;
}

//
// input/output
//
void
gd_curv_metric_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/metric_%s.nc", in_dir, fname_coords);
  
  // read in nc
  int ncid, varid;
  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR) {
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
  
  for (int ivar=0; ivar<number_of_vars; ivar++) {
      ierr = nc_inq_varid(ncid, g3d_name[ivar], &varid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  
      ierr = nc_get_var_float(ncid,varid,g3d+g3d_pos[ivar]);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
}

void
gd_curv_coord_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords)
{
  // construct file name
  char in_file[CONST_MAX_STRLEN];
  sprintf(in_file, "%s/coord_%s.nc", in_dir, fname_coords);
  
  // read in nc
  int ncid, varid;
  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
  
  for (int ivar=0; ivar<number_of_vars; ivar++) {
      ierr = nc_inq_varid(ncid, g3d_name[ivar], &varid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  
      ierr = nc_get_var_float(ncid,varid,g3d+g3d_pos[ivar]);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
}

void
gd_curv_coord_export(
  gdinfo_t *gdinfo,
  gd_t *gdcurv,
  char *fname_coords,
  char *output_dir)
{
  size_t *restrict c3d_pos   = gdcurv->cmp_pos;
  char  **restrict c3d_name  = gdcurv->cmp_name;
  int number_of_vars = gdcurv->ncmp;
  int  nx = gdcurv->nx;
  int  ny = gdcurv->ny;
  int  nz = gdcurv->nz;
  int  ni1 = gdinfo->ni1;
  int  nj1 = gdinfo->nj1;
  int  nk1 = gdinfo->nk1;
  int  ni  = gdinfo->ni;
  int  nj  = gdinfo->nj;
  int  nk  = gdinfo->nk;
  int  gni1 = gdinfo->ni1_to_glob_phys0;
  int  gnj1 = gdinfo->nj1_to_glob_phys0;
  int  gnk1 = gdinfo->nk1_to_glob_phys0;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord_%s.nc", output_dir, fname_coords);
  
  // read in nc
  int ncid;
  int varid[gdcurv->ncmp];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    ierr = nc_def_var(ncid, gdcurv->cmp_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nj1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int g_start[] = { gni1, gnj1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,g_start);

  int l_count[] = { ni, nj, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<gdcurv->ncmp; ivar++) {
    float *ptr = gdcurv->v4d + gdcurv->cmp_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return;
}

void
gd_cart_coord_export(
  gdinfo_t *gdinfo,
  gd_t *gdcart,
  char *fname_coords,
  char *output_dir)
{
  int  nx = gdcart->nx;
  int  ny = gdcart->ny;
  int  nz = gdcart->nz;
  int  ni1 = gdinfo->ni1;
  int  nj1 = gdinfo->nj1;
  int  nk1 = gdinfo->nk1;
  int  ni  = gdinfo->ni;
  int  nj  = gdinfo->nj;
  int  nk  = gdinfo->nk;
  int  gni1 = gdinfo->ni1_to_glob_phys0;
  int  gnj1 = gdinfo->nj1_to_glob_phys0;
  int  gnk1 = gdinfo->nk1_to_glob_phys0;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/coord_%s.nc", output_dir, fname_coords);
  
  // read in nc
  int ncid;
  int varid[CONST_NDIM];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  ierr = nc_def_var(ncid, "x", NC_FLOAT, 1, dimid+2, &varid[0]);
  ierr = nc_def_var(ncid, "y", NC_FLOAT, 1, dimid+1, &varid[1]);
  ierr = nc_def_var(ncid, "z", NC_FLOAT, 1, dimid+0, &varid[2]);

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nj1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int g_start[] = { gni1, gnj1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,g_start);

  int l_count[] = { ni, nj, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  ierr = nc_put_var_float(ncid, varid[0], gdcart->x1d);
  ierr = nc_put_var_float(ncid, varid[1], gdcart->y1d);
  ierr = nc_put_var_float(ncid, varid[2], gdcart->z1d);
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return;
}

void
gd_curv_metric_export(gdinfo_t        *gdinfo,
                      gdcurv_metric_t *metric,
                      char *fname_coords,
                      char *output_dir)
{
  size_t *restrict g3d_pos   = metric->cmp_pos;
  char  **restrict g3d_name  = metric->cmp_name;
  int  number_of_vars = metric->ncmp;
  int  nx = metric->nx;
  int  ny = metric->ny;
  int  nz = metric->nz;
  int  ni1 = gdinfo->ni1;
  int  nj1 = gdinfo->nj1;
  int  nk1 = gdinfo->nk1;
  int  ni  = gdinfo->ni;
  int  nj  = gdinfo->nj;
  int  nk  = gdinfo->nk;
  int  gni1 = gdinfo->ni1_to_glob_phys0;
  int  gnj1 = gdinfo->nj1_to_glob_phys0;
  int  gnk1 = gdinfo->nk1_to_glob_phys0;

  // construct file name
  char ou_file[CONST_MAX_STRLEN];
  sprintf(ou_file, "%s/metric_%s.nc", output_dir, fname_coords);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, g3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nj1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,l_start);

  int g_start[] = { gni1, gnj1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,CONST_NDIM,g_start);

  int l_count[] = { ni, nj, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,CONST_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = metric->v4d + g3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
}

/*
 * generate grid from a input grid layer file
 * code by SunWenliang  2021/06/14
 */
int
gd_curv_gen_layer(char *in_grid_layer_file,
                      int *grid_layer_resample_factor,
                      int *grid_layer_start,
                      int n_total_grid_x,
                      int n_total_grid_y,
                      int n_total_grid_z,
                      float *restrict x3d,
                      float *restrict y3d,
                      float *restrict z3d,
                      int nx, int ni, int gni1, int fdx_nghosts, 
                      int ny, int nj, int gnj1, int fdy_nghosts, 
                      int nz, int nk, int gnk1, int fdz_nghosts)
{
  int ierr = 0;

  int nLayers;
  int n_Interfaces;
  int nx_layers;
  int ny_layers;
  int nx_first;
  int ny_first;
  int nz_first;
  int nx_interp;
  int ny_interp;
  int nz_interp;
  size_t iptr1 ;
  size_t iptr2 ;
  size_t iptr3 ;

//
// read data from layer file
//

  FILE * fp;
  fp = fopen( in_grid_layer_file, "r");
  if (!fp){
    fprintf(stderr,"Failed to open input file of interface!!\n");
    exit(-1);
  }

  // read number of interface
  fscanf(fp,"%d", &n_Interfaces);

  nLayers = n_Interfaces - 1;
  int NCellPerlay_IN[nLayers];
  int NCellPerlay[nLayers];
  int VmapSpacingIsequal[nLayers];

  // read cells and is_equal of each layer
  for ( int i=0; i<nLayers; i++)
  {
    fscanf(fp,"%d",&NCellPerlay_IN[i]);
  }
  for ( int i=0; i<nLayers; i++)
  {
    fscanf(fp,"%d",&VmapSpacingIsequal[i]);
  }

  // read number of horizontal sampling
  fscanf(fp,"%d",&nx_layers);
  fscanf(fp,"%d",&ny_layers);

  size_t siz_volume_layerIn = nx_layers  * ny_layers  * (nLayers+1) ;  
  float * layer3d_In        = NULL;
  layer3d_In = (float *)fdlib_mem_malloc_1d(siz_volume_layerIn * 3 * sizeof(float), 
                                                              "gd_curv_gen_layer");
  for ( int i=0; i<siz_volume_layerIn; i++)
  {
    // read x,y,z of each sample
    fscanf(fp,"%f",&layer3d_In[i                       ]);
    fscanf(fp,"%f",&layer3d_In[i + siz_volume_layerIn  ]);
    fscanf(fp,"%f",&layer3d_In[i + siz_volume_layerIn*2]);
  }
  fclose( fp );

//
// resample of input interface grid nodes
//

  int x_interp_factor = abs(grid_layer_resample_factor[0]);
  int y_interp_factor = abs(grid_layer_resample_factor[1]);
  int z_interp_factor = abs(grid_layer_resample_factor[2]);

  // effective layer nx after downsampling
  if ( grid_layer_resample_factor[0] < 0 ) {
    nx_interp = floor((nx_layers              + 1) / x_interp_factor);
    nx_first  = floor((grid_layer_start[0] + 1) / x_interp_factor);
  }
  // effective layer nx after upsampling
  else
  { 
    nx_interp = (nx_layers               -1) * x_interp_factor + 1;
    nx_first  = (grid_layer_start[0]-1) * x_interp_factor + 1;
  }

  if ( grid_layer_resample_factor[1] < 0 ) 
  {
    ny_interp = floor( (ny_layers               +1) / y_interp_factor );
    ny_first  = floor( (grid_layer_start[1]+1) / y_interp_factor );
  }
  else
  { 
    ny_interp = (ny_layers               -1) * y_interp_factor + 1;
    ny_first  = (grid_layer_start[1]-1) * y_interp_factor + 1;
  }

  nz_interp = 1;
  if ( grid_layer_resample_factor[2] < 0 ) 
  {
    for ( int kk=0; kk<nLayers; kk++)
    {
      NCellPerlay[kk] = NCellPerlay_IN[kk] / z_interp_factor ;
      nz_interp = nz_interp + NCellPerlay[kk];
    }
  }
  else
  { 
    for ( int kk=0; kk<nLayers; kk++)
    {
      NCellPerlay[kk] = NCellPerlay_IN[kk] * z_interp_factor ;
      nz_interp = nz_interp + NCellPerlay[kk];
    }
  }

  nz_first = nz_interp - grid_layer_start[2] - n_total_grid_z;
  int x_gd_first = nx_first + gni1 - fdx_nghosts; //Starting x index at input interface
  int y_gd_first = ny_first + gnj1 - fdy_nghosts;
  int z_gd_first = nz_first + gnk1 - fdz_nghosts;

  if (gni1 + gnj1 == 0)
  {
    fprintf(stdout, "-------------------------------------------------------\n");
    fprintf(stdout, "--> Input layer information.\n");
    fprintf(stdout, "-------------------------------------------------------\n");
    fprintf(stdout, "n_Interfaces = %d\n", n_Interfaces);
    fprintf(stdout, "nx_layers = %d\nny_layers = %d\n", nx_layers, ny_layers);

    fprintf(stdout, "     resample_factor = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", grid_layer_resample_factor[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "         NCellPerlay = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", NCellPerlay_IN[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "NCellPerlay_resample = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", NCellPerlay[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "  VmapSpacingIsequal = [");
    for (int ii = 0; ii < nLayers; ii++) {
      fprintf(stdout, " %d", VmapSpacingIsequal[ii]);
    }
    fprintf(stdout, " ]\n");

    fprintf(stdout, "n_input_resample = [ %d %d %d ]\n", nx_interp, ny_interp, nz_interp);
    fprintf(stdout, "          n_c3d  = [ %d %d %d ]\n", n_total_grid_x + fdx_nghosts * 2, 
                      n_total_grid_y + fdy_nghosts * 2, n_total_grid_z + fdz_nghosts * 2);
    fprintf(stdout, "horizontal_index = [ %d %d ] --> [ %d %d ]\n", x_gd_first, y_gd_first, 
    x_gd_first + n_total_grid_x + fdx_nghosts * 2, y_gd_first + n_total_grid_y + fdy_nghosts * 2);
    fprintf(stdout, "  vertical_index = %d --> %d \n", z_gd_first, z_gd_first + n_total_grid_z+ 
                                                                           fdz_nghosts);
    fprintf(stdout, "-------------------------------------------------------\n");
    fflush(stdout);
  }

  //chech interface parameters
  if (nx_first - fdx_nghosts < 0 || ny_first - fdy_nghosts < 0)
  {
    fprintf(stdout, "Input Parameter Error: horizontal_start_index.\n");
    fflush(stdout);
    fprintf(stdout, "horizontal_start_index is related to fdx_nghosts, fdy_nghosts and refine_factor\n");
    fflush(stdout);
    exit(-1);
  }
  else if (nx_first + n_total_grid_x + fdx_nghosts > nx_interp || ny_first + n_total_grid_y + fdy_nghosts > ny_interp)
  {
    fprintf(stdout, "Input Parameter Error:\n"); 
    fflush(stdout);
    fprintf(stdout, "please check the number_of_total_grid_points， refine_factor and in_grid_layer_file file !!\n"); 
    fflush(stdout);
    exit(-1);
  }
  else if ( z_gd_first < 0)
  {
    fprintf(stdout, "Input Parameter Error:\n"); 
    fflush(stdout);
    fprintf(stdout, "please check the vertical_FreeSurf_AfterRefine_start_index\n number_of_total_grid_point_z\n refine_factor \nin_grid_layer_file file\n"); 
    fflush(stdout);
    exit(-1);
  }

  // resample based on interp_factor on input grid layer
  size_t siz_volume_layer_interp = nx_interp * ny_interp * (nLayers + 1);
  float *layer3d_interp = NULL;
  layer3d_interp = (float *)fdlib_mem_malloc_1d(siz_volume_layer_interp * 3 * sizeof(float), 
                                                                        "gd_curv_gen_layer");
  
  // 
  float   *xi_lenX = NULL;
  float   *xn_lenX = NULL;
  float *yi_x_lenX = NULL;
  float *yn_x_lenX = NULL;
  float *yi_y_lenX = NULL;
  float *yn_y_lenX = NULL;
  float *yi_z_lenX = NULL;
  float *yn_z_lenX = NULL;
  float   *xi_lenY = NULL;
  float   *xn_lenY = NULL;
  float *yi_x_lenY = NULL;
  float *yn_x_lenY = NULL;
  float *yi_y_lenY = NULL;
  float *yn_y_lenY = NULL;
  float *yi_z_lenY = NULL;
  float *yn_z_lenY = NULL;

  xi_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  xn_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");
  yi_x_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  yn_x_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");
  yi_y_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  yn_y_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");
  yi_z_lenX = (float *)fdlib_mem_malloc_1d(nx_interp * sizeof(float), "gd_curv_gen_layer");
  yn_z_lenX = (float *)fdlib_mem_malloc_1d(nx_layers * sizeof(float), "gd_curv_gen_layer");

  xi_lenY = (float *)fdlib_mem_malloc_1d(ny_interp * sizeof(float), "gd_curv_gen_layer");
  xn_lenY = (float *)fdlib_mem_malloc_1d(ny_layers * sizeof(float), "gd_curv_gen_layer");
  yi_x_lenY = (float *)fdlib_mem_malloc_1d(ny_interp * sizeof(float), "gd_curv_gen_layer");
  yn_x_lenY = (float *)fdlib_mem_malloc_1d(ny_layers * sizeof(float), "gd_curv_gen_layer");
  yi_y_lenY = (float *)fdlib_mem_malloc_1d(ny_interp * sizeof(float), "gd_curv_gen_layer");
  yn_y_lenY = (float *)fdlib_mem_malloc_1d(ny_layers * sizeof(float), "gd_curv_gen_layer");
  yi_z_lenY = (float *)fdlib_mem_malloc_1d(ny_interp * sizeof(float), "gd_curv_gen_layer");
  yn_z_lenY = (float *)fdlib_mem_malloc_1d(ny_layers * sizeof(float), "gd_curv_gen_layer");

  // Downsampling in X and Y
  if (grid_layer_resample_factor[0] < 2 && grid_layer_resample_factor[1] < 2 ) 
  {
    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int jj = 0; jj < ny_interp; jj++) {
        for (int ii = 0; ii < nx_interp; ii++) {
          iptr1 = M_gd_INDEX(ii, jj, kk, nx_interp, ny_interp);
          iptr2 = M_gd_INDEX(ii * x_interp_factor, jj * y_interp_factor, kk, nx_layers,
                        ny_layers);
          layer3d_interp[ iptr1                           ] = layer3d_In[ iptr2 ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp  ] = layer3d_In[ iptr2 
                          +siz_volume_layerIn   ]; 
          layer3d_interp[ iptr1 +siz_volume_layer_interp*2] = layer3d_In[ iptr2 
                          +siz_volume_layerIn*2 ]; 
        }
      }
    }
  }
  //Downsampling in in X, Interpolating in Y
  else if( grid_layer_resample_factor[0] < 2 && grid_layer_resample_factor[1] >= 2  )  
  {
    for ( int jj1=0; jj1<ny_interp; jj1++) xi_lenY[jj1] = jj1;
    for ( int jj1=0; jj1<ny_layers; jj1++) xn_lenY[jj1] = jj1*y_interp_factor;

    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int ii = 0; ii < nx_interp; ii++) {
        for ( int jj1=0; jj1<ny_layers; jj1++)
        {
          yn_x_lenY[jj1] = layer3d_In[M_gd_INDEX(ii*x_interp_factor, jj1, kk, nx_layers,
                                       ny_layers)                        ];
          yn_y_lenY[jj1] = layer3d_In[M_gd_INDEX(ii*x_interp_factor, jj1, kk, nx_layers,
                                             ny_layers) + siz_volume_layerIn   ];
          yn_z_lenY[jj1] = layer3d_In[M_gd_INDEX(ii*x_interp_factor, jj1, kk, nx_layers,
                                             ny_layers) + siz_volume_layerIn*2 ];
        }

        gd_SPL(ny_layers, xn_lenY, yn_x_lenY, ny_interp, xi_lenY, yi_x_lenY);
        gd_SPL(ny_layers, xn_lenY, yn_y_lenY, ny_interp, xi_lenY, yi_y_lenY);
        gd_SPL(ny_layers, xn_lenY, yn_z_lenY, ny_interp, xi_lenY, yi_z_lenY);

        for (int jj = 0; jj < ny_interp; jj++)
        {
          iptr1 = M_gd_INDEX(ii, jj, kk, nx_interp, ny_interp);
          layer3d_interp[iptr1] = yi_x_lenY[jj];
          layer3d_interp[iptr1 + siz_volume_layer_interp] = yi_y_lenY[jj];
          layer3d_interp[iptr1 + siz_volume_layer_interp * 2] = yi_z_lenY[jj];
        }
      }
    }
  }
  //Interpolating in in X, Downsampling in Y
  else if( grid_layer_resample_factor[0] >=2 && grid_layer_resample_factor[1] < 2  )  
  {
    for ( int ii1=0; ii1<nx_interp; ii1++) xi_lenX[ii1] = ii1;
    for ( int ii1=0; ii1<nx_layers; ii1++) xn_lenX[ii1] = ii1*x_interp_factor;

    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int jj = 0; jj < ny_interp; jj++) {
        for ( int ii1=0; ii1<nx_layers; ii1++)
        {
          yn_x_lenX[ii1] = layer3d_In[M_gd_INDEX(ii1, jj*y_interp_factor, kk, nx_layers,
                                            ny_layers)                        ];
          yn_y_lenX[ii1] = layer3d_In[M_gd_INDEX(ii1, jj*y_interp_factor, kk, nx_layers,
                                            ny_layers) + siz_volume_layerIn   ];
          yn_z_lenX[ii1] = layer3d_In[M_gd_INDEX(ii1, jj*y_interp_factor, kk, nx_layers,
                                            ny_layers) + siz_volume_layerIn*2 ];
        }

        gd_SPL(nx_layers, xn_lenX, yn_x_lenX, nx_interp, xi_lenX, yi_x_lenX);
        gd_SPL(nx_layers, xn_lenX, yn_y_lenX, nx_interp, xi_lenX, yi_y_lenX);
        gd_SPL(nx_layers, xn_lenX, yn_z_lenX, nx_interp, xi_lenX, yi_z_lenX);

        for (int ii = 0; ii < nx_interp; ii++)
        {
            iptr1 = M_gd_INDEX(ii, jj, kk, nx_interp, ny_interp);
            layer3d_interp[iptr1] = yi_x_lenX[ii];
            layer3d_interp[iptr1 + siz_volume_layer_interp] = yi_y_lenX[ii];
            layer3d_interp[iptr1 + siz_volume_layer_interp * 2] = yi_z_lenX[ii];
        }
      }
    }
  }
  // Interpolating in in X and Y
  else if( grid_layer_resample_factor[0] >=2 && grid_layer_resample_factor[1] >= 2  ) 
  {
    // interp x
    for ( int ii1=0; ii1<nx_interp; ii1++) xi_lenX[ii1] = ii1;
    for ( int ii1=0; ii1<nx_layers; ii1++) xn_lenX[ii1] = ii1*x_interp_factor;
    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int jj = 0; jj < ny_layers; jj++) {
        for ( int ii1=0; ii1<nx_layers; ii1++)
        {
          yn_x_lenX[ii1] = layer3d_In[M_gd_INDEX(ii1, jj, kk, nx_layers, ny_layers)];
          yn_y_lenX[ii1] = layer3d_In[M_gd_INDEX(ii1, jj, kk, nx_layers, ny_layers) 
                                       + siz_volume_layerIn   ];
          yn_z_lenX[ii1] = layer3d_In[M_gd_INDEX(ii1, jj, kk, nx_layers, ny_layers) 
                                       + siz_volume_layerIn*2 ];
        }

        gd_SPL( nx_layers, xn_lenX, yn_x_lenX, nx_interp, xi_lenX, yi_x_lenX);
        gd_SPL( nx_layers, xn_lenX, yn_y_lenX, nx_interp, xi_lenX, yi_y_lenX);
        gd_SPL( nx_layers, xn_lenX, yn_z_lenX, nx_interp, xi_lenX, yi_z_lenX);

        for (int ii = 0; ii < nx_interp; ii++)
        {
          iptr1 = M_gd_INDEX( ii, jj*y_interp_factor, kk, nx_interp, ny_interp );
          layer3d_interp[ iptr1                           ] = yi_x_lenX[ ii ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp  ] = yi_y_lenX[ ii ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp*2] = yi_z_lenX[ ii ];
        }
      }
    }

    // interp y
    for ( int jj1=0; jj1<ny_interp; jj1++) xi_lenY[jj1] = jj1;
    for ( int jj1=0; jj1<ny_layers; jj1++) xn_lenY[jj1] = jj1*y_interp_factor;
    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int ii = 0; ii < nx_interp; ii++) {
        for ( int jj1=0; jj1<ny_layers; jj1++) {
          yn_x_lenY[jj1] = layer3d_interp[M_gd_INDEX(ii, jj1 * y_interp_factor, kk,
                               nx_interp, ny_interp)];
          yn_y_lenY[jj1] = layer3d_interp[M_gd_INDEX(ii, jj1 * y_interp_factor, kk,
                               nx_interp, ny_interp) + siz_volume_layer_interp];
          yn_z_lenY[jj1] = layer3d_interp[M_gd_INDEX(ii, jj1 * y_interp_factor, kk,
                               nx_interp, ny_interp) + siz_volume_layer_interp * 2];
          }
          gd_SPL( ny_layers, xn_lenY, yn_x_lenY, ny_interp, xi_lenY, yi_x_lenY);
          gd_SPL( ny_layers, xn_lenY, yn_y_lenY, ny_interp, xi_lenY, yi_y_lenY);
          gd_SPL( ny_layers, xn_lenY, yn_z_lenY, ny_interp, xi_lenY, yi_z_lenY);
        for (int jj = 0; jj < ny_interp; jj++) {
          iptr1 = M_gd_INDEX( ii, jj, kk, nx_interp, ny_interp );
          layer3d_interp[ iptr1                           ] = yi_x_lenY[ jj ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp  ] = yi_y_lenY[ jj ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp*2] = yi_z_lenY[ jj ];
        }
      }
    }  
  }

  free(xi_lenX);
  free(xn_lenX);
  free(yi_x_lenX);
  free(yn_x_lenX);
  free(yi_y_lenX);
  free(yn_y_lenX);
  free(yi_z_lenX);
  free(yn_z_lenX);
  free(xi_lenY);
  free(xn_lenY);
  free(yi_x_lenY);
  free(yn_x_lenY);
  free(yi_y_lenY);
  free(yn_y_lenY);
  free(yi_z_lenY);
  free(yn_z_lenY);

//
// generate FD grid
//

  size_t siz_volume = nx * ny * nz;
  float *layer3d = NULL;
  layer3d = ( float * ) malloc( sizeof( float ) * nx * ny * (nLayers+1) * 3 );

  // suffix ch means:
  float zdiff_two_interface; // thickness between interfaces
  float zdiff_two_interface_ch;

  float *xlayer3d = layer3d + 0 * nx*ny*(nLayers+1);
  float *ylayer3d = layer3d + 1 * nx*ny*(nLayers+1);
  float *zlayer3d = layer3d + 2 * nx*ny*(nLayers+1);
  
  float * zlayerpart  = NULL;
  float * ylayerpart  = NULL;
  float * xlayerpart  = NULL;
  zlayerpart = (float *)fdlib_mem_malloc_1d((nLayers+1) * sizeof(float), "gd_curv_gen_layer");
  ylayerpart = (float *)fdlib_mem_malloc_1d((nLayers+1) * sizeof(float), "gd_curv_gen_layer");
  xlayerpart = (float *)fdlib_mem_malloc_1d((nLayers+1) * sizeof(float), "gd_curv_gen_layer");

  float * z3dpart  = NULL;
  float * y3dpart  = NULL;
  float * x3dpart  = NULL;
  z3dpart = (float *)fdlib_mem_malloc_1d((nz_interp+fdz_nghosts) * sizeof(float), "gd_curv_gen_layer");
  y3dpart = (float *)fdlib_mem_malloc_1d((nz_interp+fdz_nghosts) * sizeof(float), "gd_curv_gen_layer");
  x3dpart = (float *)fdlib_mem_malloc_1d((nz_interp+fdz_nghosts) * sizeof(float), "gd_curv_gen_layer");

  // layer regions of input model covered by this thread
  for (int k = 0; k < (nLayers + 1); k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++)
      {
        iptr1 = M_gd_INDEX(i, j, k, nx, ny);
        iptr2 = M_gd_INDEX(i + x_gd_first, j + y_gd_first, k, nx_interp, ny_interp);
        xlayer3d[iptr1] = layer3d_interp[iptr2];
        ylayer3d[iptr1] = layer3d_interp[iptr2 + siz_volume_layer_interp];
        zlayer3d[iptr1] = layer3d_interp[iptr2 + siz_volume_layer_interp * 2];
      }
    }
  }

  // interp both z and x,y
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      for (int k = 0; k < nLayers + 1; k++)
      {
        xlayerpart[k] = xlayer3d[M_gd_INDEX(i, j, k, nx, ny)];
        ylayerpart[k] = ylayer3d[M_gd_INDEX(i, j, k, nx, ny)];
        zlayerpart[k] = zlayer3d[M_gd_INDEX(i, j, k, nx, ny)];
      }
      // Interpolating the Z coordinate of the grid 
      gd_grid_z_interp(z3dpart, zlayerpart, NCellPerlay, VmapSpacingIsequal,
                       nLayers, nx, ny);

      // Interpolating the X and Y coordinate of the grid ...
      // by cubic spline interpolation method.
      gd_SPL(nLayers + 1, zlayerpart, xlayerpart, nz_interp, z3dpart, x3dpart);
      gd_SPL(nLayers + 1, zlayerpart, ylayerpart, nz_interp, z3dpart, y3dpart);

      //Grids outside the Free Surface.
      for (int k = 1; k < fdz_nghosts + 1; k++)
      {
        iptr1 = nz_interp - 1;
        x3dpart[iptr1 + k] = x3dpart[iptr1] * 2 - x3dpart[iptr1 - k];
        y3dpart[iptr1 + k] = y3dpart[iptr1] * 2 - y3dpart[iptr1 - k];
        z3dpart[iptr1 + k] = z3dpart[iptr1] * 2 - z3dpart[iptr1 - k];
      }

      for ( int k = 0; k < nk + fdz_nghosts*2; k ++)
      {
        iptr1 = z_gd_first + k;
        x3d[M_gd_INDEX(i, j, k, nx, ny)] = x3dpart[iptr1];
        y3d[M_gd_INDEX(i, j, k, nx, ny)] = y3dpart[iptr1];
        z3d[M_gd_INDEX(i, j, k, nx, ny)] = z3dpart[iptr1];
      }
    }
  }

  free(layer3d);
  free(layer3d_In);
  free(layer3d_interp);
  free(zlayerpart);
  free(ylayerpart);
  free(xlayerpart);
  free(z3dpart);
  free(y3dpart);
  free(x3dpart);

  return ierr;
}

//  Interpolating the Z coordinate
int gd_grid_z_interp(float *z3dpart, float *zlayerpart, int *NCellPerlay,
                     int *VmapSpacingIsequal, int nLayers, int nx, int ny)
{
  int ierr = 0;

  // If Grid step size of Z coordinate is not equal, divide the layer into two parts.
  // N1 is first part number of this layer and N2 is second.
  int N1;
  int N2;
  float distance;  //Interface spacing
  float dmidpoint;  //The length of the N1th grid
  float N1_Isochromatic; 
  float N2_Isochromatic;
  int NCellPerlay_max = 0;
  int sumNCellPerlay = 0;
  float LayerDz[nLayers];

  for (int i = 0; i < nLayers; i++)
  {
    if (NCellPerlay_max < NCellPerlay[i])
    {
      NCellPerlay_max = NCellPerlay[i];
    }
  }

  float range1[NCellPerlay_max];
  float range2[NCellPerlay_max];

  for (int i = 0; i < nLayers; i++)
  {
    // LayerDz[i] = (zlayer3d[M_gd_INDEX(xi, yi, i + 1, nx, ny)] -
    //               zlayer3d[M_gd_INDEX(xi, yi, i, nx, ny)]) / NCellPerlay[i];
    LayerDz[i] = (zlayerpart[i + 1] - zlayerpart[i]) / NCellPerlay[i];
  }

  for (int i = 0; i < nLayers; i++)
  {
    // The grid spacing is equal
    z3dpart[sumNCellPerlay] = zlayerpart[i];
    for (int ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i] - 1; ii++)
    {
      z3dpart[ii+1] = zlayerpart[i] + (ii - sumNCellPerlay + 1) * LayerDz[i];
    }
    // The grid spacing is not equal
    if (VmapSpacingIsequal[i] < 1 && i > 0 && i < nLayers - 1)
    {
      range1[0] = zlayerpart[i];
      if ((LayerDz[i] - LayerDz[i - 1]) * (LayerDz[i + 1] - LayerDz[i]) > 0)
      {
        N1 = ((LayerDz[i] - LayerDz[i - 1]) / (LayerDz[i + 1] - LayerDz[i - 1])) * (NCellPerlay[i] + 2) + 1;
        if (N1 > NCellPerlay[i] - 1)
          N1 = NCellPerlay[i] - 1;
        N2 = NCellPerlay[i] + 2 - N1;

        distance = zlayerpart[i + 1] - zlayerpart[i] + LayerDz[i - 1] + LayerDz[i + 1];
        dmidpoint = (2 * distance - LayerDz[i - 1] * N1 - LayerDz[i + 1] * (N2 + 1)) / (N1 + N2 - 1);

        N1_Isochromatic = (dmidpoint - LayerDz[i - 1]) / (N1 - 1);
        N2_Isochromatic = (LayerDz[i + 1] - dmidpoint) / (N2);

        for (int j = 1; j < N1; j++)
        {
          range1[j] = range1[j - 1] + LayerDz[i - 1] + j * N1_Isochromatic;
        }
        for (int j = N1; j < N1 + N2; j++)
        {
          range1[j] = range1[j - 1] + LayerDz[i + 1] - (N2 + N1 - j - 1) * N2_Isochromatic;
        }
        for (int ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i]; ii++)
        {
          z3dpart[ii] = range1[ii - sumNCellPerlay];
        }
      }
    }
    sumNCellPerlay += NCellPerlay[i];
  }
  z3dpart[sumNCellPerlay] = zlayerpart[nLayers];

  return ierr;
}

float gd_seval(int ni, float u,
            int n, float x[], float y[],
            float b[], float c[], float d[],
            int *last)
{
  int i, j, k;
  float w;
  i = *last;
  if (i >= n - 1) i = 0;
  if (i < 0) i = 0;
  if ((x[i] > u) || (x[i + 1] < u))//??
  {
    i = 0;
    j = n;
    do
    {
      k = (i + j) / 2;
      if (u < x[k]) j = k;
      if (u >= x[k]) i = k;
    } while (j > i + 1);
  }
  *last = i;
  w = u - x[i];
  w = y[i] + w * (b[i] + w * (c[i] + w * d[i]));
  return (w);
}

//// Cubic spline difference function
int gd_SPLine( int n, int end1, int end2,
           float slope1, float slope2,
           float x[], float y[],
           float b[], float c[], float d[],
           int *iflag)
{
  int nm1, ib, i, ascend;
  float t;
  nm1 = n - 1;
  *iflag = 0;
  if (n < 2)
  {
    *iflag = 1;
    goto Leavegd_SPLine;
  }
  ascend = 1;
  for (i = 1; i < n; ++i) if (x[i] <= x[i - 1]) ascend = 0;
  if (!ascend)
  {
    *iflag = 2;
    goto Leavegd_SPLine;
  }
  if (n >= 3)
  {
    d[0] = x[1] - x[0];
    c[1] = (y[1] - y[0]) / d[0];
    for (i = 1; i < nm1; ++i)
    {
      d[i] = x[i + 1] - x[i];
      b[i] = 2.0 * (d[i - 1] + d[i]);
      c[i + 1] = (y[i + 1] - y[i]) / d[i];
      c[i] = c[i + 1] - c[i];
    }
    b[0] = -d[0];
    b[nm1] = -d[n - 2];
    c[0] = 0.0;
    c[nm1] = 0.0;
    if (n != 3)
    {
      c[0] = c[2] / (x[3] - x[1]) - c[1] / (x[2] - x[0]);
      c[nm1] = c[n - 2] / (x[nm1] - x[n - 3]) - c[n - 3] / (x[n - 2] - x[n - 4]);
      c[0] = c[0] * d[0] * d[0] / (x[3] - x[0]);
      c[nm1] = -c[nm1] * d[n - 2] * d[n - 2] / (x[nm1] - x[n - 4]);
    }
    if (end1 == 1)
    {
      b[0] = 2.0 * (x[1] - x[0]);
      c[0] = (y[1] - y[0]) / (x[1] - x[0]) - slope1;
    }
    if (end2 == 1)
    {
      b[nm1] = 2.0 * (x[nm1] - x[n - 2]);
      c[nm1] = slope2 - (y[nm1] - y[n - 2]) / (x[nm1] - x[n - 2]);
    }
    // Forward elimination 
    for (i = 1; i < n; ++i)
    {
      t = d[i - 1] / b[i - 1];
      b[i] = b[i] - t * d[i - 1];
      c[i] = c[i] - t * c[i - 1];
    }
    // Back substitution 
    c[nm1] = c[nm1] / b[nm1];
    for (ib = 0; ib < nm1; ++ib)
    {
      i = n - ib - 2;
      c[i] = (c[i] - d[i] * c[i + 1]) / b[i];
    }
    b[nm1] = (y[nm1] - y[n - 2]) / d[n - 2] + d[n - 2] * (c[n - 2] + 2.0 * c[nm1]);
    for (i = 0; i < nm1; ++i)
    {
      b[i] = (y[i + 1] - y[i]) / d[i] - d[i] * (c[i + 1] + 2.0 * c[i]);
      d[i] = (c[i + 1] - c[i]) / d[i];
      c[i] = 3.0 * c[i];
    }
    c[nm1] = 3.0 * c[nm1];
    d[nm1] = d[n - 2];
  }
  else
  {
    b[0] = (y[1] - y[0]) / (x[1] - x[0]);
    c[0] = 0.0;
    d[0] = 0.0;
    b[1] = b[0];
    c[1] = 0.0;
    d[1] = 0.0;
  }
Leavegd_SPLine:
  return 0;
}

//// Cubic spline difference function in grid interpolation
void gd_SPL(int n, float *x, float *y, int ni, float *xi, float *yi)
{
  float *b, *c, *d;
  int iflag=0, last=0, i=0;
  b = (float *)fdlib_mem_malloc_1d(n * sizeof(float), "gd_SPL");
  c = (float *)fdlib_mem_malloc_1d(n * sizeof(float), "gd_SPL");
  d = (float *)fdlib_mem_malloc_1d(n * sizeof(float), "gd_SPL");

  if (!d) { printf("no enough memory for b,c,d\n"); }
  else {
    gd_SPLine(n, 0, 0, 0.0, 0.0, x, y, b, c, d, &iflag);
    for (i = 0; i<ni; i++)
      yi[i] = gd_seval(ni, xi[i], n, x, y, b, c, d, &last);
    free(b);
    free(c);
    free(d);
  };
}

void
gd_curv_set_minmax(gd_t *gdcurv)
{
  float xmin = gdcurv->x3d[0], xmax = gdcurv->x3d[0];
  float ymin = gdcurv->y3d[0], ymax = gdcurv->y3d[0];
  float zmin = gdcurv->z3d[0], zmax = gdcurv->z3d[0];
  
  for (size_t i = 0; i < gdcurv->siz_icmp; i++){
      xmin = xmin < gdcurv->x3d[i] ? xmin : gdcurv->x3d[i];
      xmax = xmax > gdcurv->x3d[i] ? xmax : gdcurv->x3d[i];
      ymin = ymin < gdcurv->y3d[i] ? ymin : gdcurv->y3d[i];
      ymax = ymax > gdcurv->y3d[i] ? ymax : gdcurv->y3d[i];
      zmin = zmin < gdcurv->z3d[i] ? zmin : gdcurv->z3d[i];
      zmax = zmax > gdcurv->z3d[i] ? zmax : gdcurv->z3d[i];
  }

  gdcurv->xmin = xmin;
  gdcurv->xmax = xmax;
  gdcurv->ymin = ymin;
  gdcurv->ymax = ymax;
  gdcurv->zmin = zmin;
  gdcurv->zmax = zmax;

  return;
}

/*
 * convert cart coord to global index
 */

int
gd_cart_coord_to_glob_indx(gdinfo_t *gdinfo,
                           gd_t *gdcart,
                           float sx,
                           float sy,
                           float sz,
                           MPI_Comm comm,
                           int myid,
                           int   *ou_si, int *ou_sj, int *ou_sk,
                           float *ou_sx_inc, float *ou_sy_inc, float *ou_sz_inc)
{
  int ierr = 0;

  int si_glob = (int)( (sx - gdcart->x0_glob) / gdcart->dx + 0.5 );
  int sj_glob = (int)( (sy - gdcart->y0_glob) / gdcart->dy + 0.5 );
  int sk_glob = (int)( (sz - gdcart->z0_glob) / gdcart->dz + 0.5 );
  float sx_inc = si_glob * gdcart->dx + gdcart->x0_glob - sx;
  float sy_inc = sj_glob * gdcart->dy + gdcart->y0_glob - sy;
  float sz_inc = sk_glob * gdcart->dz + gdcart->z0_glob - sz;

  *ou_si = si_glob;
  *ou_sj = sj_glob;
  *ou_sk = sk_glob;
  *ou_sx_inc = sx_inc;
  *ou_sy_inc = sy_inc;
  *ou_sz_inc = sz_inc;

  return ierr; 
}

/*
 * convert curv coord to global index using MPI
 */

int
gd_curv_coord_to_glob_indx(gdinfo_t *gdinfo,
                           gd_t *gdcurv,
                           float sx,
                           float sy,
                           float sz,
                           MPI_Comm comm,
                           int myid,
                           int   *ou_si, int *ou_sj, int *ou_sk,
                           float *ou_sx_inc, float *ou_sy_inc, float *ou_sz_inc,
                           float *restrict wrk3d)
{
  int is_here = 0;

  int si_glob = 0;
  int sj_glob = 0;
  int sk_glob = 0;
  float sx_inc = 0.0;
  float sy_inc = 0.0;
  float sz_inc = 0.0;
  int si = 0;
  int sj = 0;
  int sk = 0;

  // if located in this thread
  is_here = gd_curv_coord_to_local_indx(gdinfo,gdcurv,sx,sy,sz,
                                    &si, &sj, &sk, &sx_inc, &sy_inc, &sz_inc,
                                    wrk3d);

  // if in this thread
  if ( is_here == 1)
  {
    // conver to global index
    si_glob = gd_info_ind_lcext2glphy_i(si, gdinfo);
    sj_glob = gd_info_ind_lcext2glphy_j(sj, gdinfo);
    sk_glob = gd_info_ind_lcext2glphy_k(sk, gdinfo);
    fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
    fprintf(stdout," -- located to global index = %d %d %d\n", 
                          si_glob, sj_glob, sk_glob);
    fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
  } else {
    //fprintf(stdout," -- not in this thread %d\n", myid);
  }

  // reduce global index and shift values
  int sendbufi = si_glob;
  MPI_Allreduce(&sendbufi, &si_glob, 1, MPI_INT, MPI_MAX, comm);

  sendbufi = sj_glob;
  MPI_Allreduce(&sendbufi, &sj_glob, 1, MPI_INT, MPI_MAX, comm);

  sendbufi = sk_glob;
  MPI_Allreduce(&sendbufi, &sk_glob, 1, MPI_INT, MPI_MAX, comm);

  float sendbuf = sx_inc;
  MPI_Allreduce(&sendbuf, &sx_inc, 1, MPI_INT, MPI_SUM, comm);

  sendbuf = sy_inc;
  MPI_Allreduce(&sendbuf, &sy_inc, 1, MPI_INT, MPI_SUM, comm);

  sendbuf = sz_inc;
  MPI_Allreduce(&sendbuf, &sz_inc, 1, MPI_INT, MPI_SUM, comm);

  fprintf(stdout," --myid=%d,index=%d %d %d,shift = %f %f %f\n",
      myid,si_glob,sj_glob,sk_glob, sx_inc,sy_inc,sz_inc);

  *ou_si = si_glob;
  *ou_sj = sj_glob;
  *ou_sk = sk_glob;
  *ou_sx_inc = sx_inc;
  *ou_sy_inc = sy_inc;
  *ou_sz_inc = sz_inc;

  return is_here; 
}

/* 
 * if the nearest point in this thread then search its grid index
 *   return value:
 *      1 - in this thread
 *      0 - not in this thread
 */

int
gd_curv_coord_to_local_indx(gdinfo_t *gdinfo,
                        gd_t *gd,
                        float sx, float sy, float sz,
                        int *si, int *sj, int *sk,
                        float *sx_inc, float *sy_inc, float *sz_inc,
                        float *restrict wrk3d)
{
  int is_here = 0; // default outside

  int nx = gdinfo->nx;
  int ny = gdinfo->ny;
  int nz = gdinfo->nz;
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line  = gdinfo->siz_iy;
  size_t siz_slice = gdinfo->siz_iz;

  float *restrict x3d = gd->x3d;
  float *restrict y3d = gd->y3d;
  float *restrict z3d = gd->z3d;

  // outside coord range
  if ( sx < gd->xmin || sx > gd->xmax ||
       sy < gd->ymin || sy > gd->ymax ||
       sz < gd->zmin || sz > gd->zmax)
  {
    is_here = 0;
    return is_here;
  }

  // init closest point
  float min_dist = sqrtf(  (sx - x3d[0]) * (sx - x3d[0])
      + (sy - y3d[0]) * (sy - y3d[0])
      + (sz - z3d[0]) * (sz - z3d[0]) );
  int min_dist_i = 0 ;
  int min_dist_j = 0 ;
  int min_dist_k = 0 ;

  // compute distance to each point
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;

        float x = x3d[iptr];
        float y = y3d[iptr];
        float z = z3d[iptr];

        float DistInt = sqrtf(  (sx - x) * (sx - x)
            + (sy - y) * (sy - y)
            + (sz - z) * (sz - z) );
        wrk3d[iptr] =  DistInt;

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = i;
          min_dist_j = j;
          min_dist_k = k;
        }
      }
    }
  }

  // if nearest index is outside phys region, not here
  if ( min_dist_i < ni1 || min_dist_i > ni2 ||
      min_dist_j < nj1 || min_dist_j > nj2 ||
      min_dist_k < nk1 || min_dist_k > nk2 )
  {
    is_here = 0;
    return is_here;
  }

  // in this thread
  is_here = 1;

  float points_x[8];
  float points_y[8];
  float points_z[8];
  float points_i[8];
  float points_j[8];
  float points_k[8];

  for (int kk=0; kk<2; kk++)
  {
    for (int jj=0; jj<2; jj++)
    {
      for (int ii=0; ii<2; ii++)
      {
        int cur_i = min_dist_i-1+ii;
        int cur_j = min_dist_j-1+jj;
        int cur_k = min_dist_k-1+kk;

        for (int n3=0; n3<2; n3++) {
          for (int n2=0; n2<2; n2++) {
            for (int n1=0; n1<2; n1++) {
              int iptr_cube = n1 + n2 * 2 + n3 * 4;
              int iptr = (cur_i+n1) + (cur_j+n2) * siz_line +
                (cur_k+n3) * siz_slice;
              points_x[iptr_cube] = x3d[iptr];
              points_y[iptr_cube] = y3d[iptr];
              points_z[iptr_cube] = z3d[iptr];
              points_i[iptr_cube] = cur_i+n1;
              points_j[iptr_cube] = cur_j+n2;
              points_k[iptr_cube] = cur_k+n3;
            }
          }
        }

        if (isPointInHexahedron(sx,sy,sz,points_x,points_y,points_z) == true)
        {
          float si_curv, sj_curv, sk_curv;
          //gd_curv_coord2index_rdinterp(sx,sy,sz,
          //              8,
          //              points_x,points_y,points_z,
          //              points_i,points_j,points_k,
          //              &si_curv, &sj_curv, &sk_curv);

          gd_curv_coord2index_sample(sx,sy,sz,
              8,
              points_x,points_y,points_z,
              points_i,points_j,points_k,
              100,100,100,
              &si_curv, &sj_curv, &sk_curv);

          // convert to return values
          *si = min_dist_i;
          *sj = min_dist_j;
          *sk = min_dist_k;
          *sx_inc = si_curv - min_dist_i;
          *sy_inc = sj_curv - min_dist_j;
          *sz_inc = sk_curv - min_dist_k;

          return is_here;
        }
      }
    }
  }

  // if not in any cube due to bug, set default value
  //    if everything is right, should be return 10 line before
  *si = min_dist_i;
  *sj = min_dist_j;
  *sk = min_dist_k;
  *sx_inc = 0.0;
  *sy_inc = 0.0;
  *sz_inc = 0.0;

  return is_here;
}

/* 
 * find curv index using sampling
 */

  int
gd_curv_coord2index_sample(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_j,
    float *points_k,
    int    nx_sample,
    int    ny_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sj_curv,
    float *sk_curv)
{
  float Lx[2], Ly[2], Lz[2];

  // init closest point
  float min_dist = sqrtf(  (sx - points_x[0]) * (sx - points_x[0])
      + (sy - points_y[0]) * (sy - points_y[0])
      + (sz - points_z[0]) * (sz - points_z[0]) );
  int min_dist_i = 0 ;
  int min_dist_j = 0 ;
  int min_dist_k = 0 ;

  // linear interp for all sample
  for (int n3=0; n3<nz_sample+1; n3++)
  {
    Lz[1] = (float)(n3) / (float)(nz_sample);
    Lz[0] = 1.0 - Lz[1];
    for (int n2=0; n2<ny_sample+1; n2++)
    {
      Ly[1] = (float)(n2) / (float)(ny_sample);
      Ly[0] = 1.0 - Ly[1];
      for (int n1=0; n1<nx_sample+1; n1++)
      {
        Lx[1] = (float)(n1) / (float)(nx_sample);
        Lx[0] = 1.0 - Lx[1];

        // interp
        float x_pt=0;
        float y_pt=0;
        float z_pt=0;
        for (int kk=0; kk<2; kk++) {
          for (int jj=0; jj<2; jj++) {
            for (int ii=0; ii<2; ii++)
            {
              int iptr_cube = ii + jj * 2 + kk * 4;
              x_pt += Lx[ii]*Ly[jj]*Lz[kk] * points_x[iptr_cube];
              y_pt += Lx[ii]*Ly[jj]*Lz[kk] * points_y[iptr_cube];
              z_pt += Lx[ii]*Ly[jj]*Lz[kk] * points_z[iptr_cube];
            }
          }
        }

        // find min dist
        float DistInt = sqrtf(  (sx - x_pt) * (sx - x_pt)
            + (sy - y_pt) * (sy - y_pt)
            + (sz - z_pt) * (sz - z_pt) );

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = n1;
          min_dist_j = n2;
          min_dist_k = n3;
        }
      } // n1
    } // n2
  } // n3

  *si_curv = points_i[0] + (float)min_dist_i / (float)nx_sample;
  *sj_curv = points_j[0] + (float)min_dist_j / (float)ny_sample;
  *sk_curv = points_k[0] + (float)min_dist_k / (float)nz_sample;

  return 0;
}

/* 
 * interp curv coord using inverse distance interp
 */

  int
gd_curv_coord2index_rdinterp(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_j,
    float *points_k,
    float *si_curv, // interped curv coord
    float *sj_curv,
    float *sk_curv)
{
  float weight[num_points];
  float total_weight = 0.0 ;

  // cal weight
  int at_point_indx = -1;
  for (int i=0; i<num_points; i++)
  {
    float dist = sqrtf ((sx - points_x[i]) * (sx - points_x[i])
        + (sy - points_y[i]) * (sy - points_y[i])
        + (sz - points_z[i]) * (sz - points_z[i])
        );
    if (dist < 1e-9) {
      at_point_indx = i;
    } else {
      weight[i]   = 1.0 / dist;
      total_weight += weight[i];
    }
  }
  // if at a point
  if (at_point_indx > 0) {
    total_weight = 1.0;
    // other weight 0
    for (int i=0; i<num_points; i++) {
      weight[i] = 0.0;
    }
    // point weight 1
    weight[at_point_indx] = 1.0;
  }

  // interp

  *si_curv = 0.0;
  *sj_curv = 0.0;
  *sk_curv = 0.0;

  for (int i=0; i<num_points; i++)
  {
    weight[i] *= 1.0 / total_weight ;

    (*si_curv) += weight[i] * points_i[i];
    (*sj_curv) += weight[i] * points_j[i]; 
    (*sk_curv) += weight[i] * points_k[i];  

    fprintf(stdout,"---- i=%d,weight=%f,points_i=%f,points_j=%f,points_k=%f\n",
        i,weight[i],points_i[i],points_j[i],points_k[i]);
  }

  return 0;
}

float
gd_coord_get_x(gd_t *gd, int i, int j, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->x1d[i];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + j * gd->siz_iy + k * gd->siz_iz;
    var = gd->x3d[iptr];
  }

  return var;
}

float
gd_coord_get_y(gd_t *gd, int i, int j, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->y1d[j];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + j * gd->siz_iy + k * gd->siz_iz;
    var = gd->y3d[iptr];
  }

  return var;
}

float
gd_coord_get_z(gd_t *gd, int i, int j, int k)
{
  float var = 0.0;

  if (gd->type == GD_TYPE_CART)
  {
    var = gd->z1d[k];
  }
  else if (gd->type == GD_TYPE_CURV)
  {
    size_t iptr = i + j * gd->siz_iy + k * gd->siz_iz;
    var = gd->z3d[iptr];
  }

  return var;
}
