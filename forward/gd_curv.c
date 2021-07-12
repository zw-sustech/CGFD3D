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
#include "gd_curv.h"

void 
gd_curv_init_c3d(
    size_t siz_volume,
    int *c3d_num_of_vars,
    float  **p_c3d,
    size_t **p_c3d_pos,
    char  ***p_c3d_name,
    char  ***p_coord_name)
{
  const int num_grid_vars = 3;
  /*
   * 0-2: x3d, y3d, z3d
   */
  
  // vars
  float *c3d = (float *) fdlib_mem_calloc_1d_float(siz_volume * num_grid_vars,
                                                   0.0,
                                                   "gd_curv_init_c3d");
  if (c3d == NULL) {
      fprintf(stderr,"Error: failed to alloc coord vars\n");
      fflush(stderr);
  }
  
  // position of each var
  size_t *c3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet(num_grid_vars,
                                                         0,
                                                         "gd_curv_init_c3d");
  
  // name of each var
  char **c3d_name = (char **) fdlib_mem_malloc_2l_char(num_grid_vars,
                                                       FD_MAX_STRLEN,
                                                       "gd_curv_init_c3d");
  
  // name of mapped coordinate
  char **coord_name = (char **) fdlib_mem_malloc_2l_char(FD_NDIM,
                                                       FD_MAX_STRLEN,
                                                       "gd_curv_init_c3d");
  
  // set value
  int ivar = GD_CURV_SEQ_X3D;
  c3d_pos[ivar] = ivar * siz_volume;
  sprintf(c3d_name[ivar],"%s","x");
  
  ivar = GD_CURV_SEQ_Y3D;
  c3d_pos[ivar] = ivar * siz_volume;
  sprintf(c3d_name[ivar],"%s","y");
  
  ivar = GD_CURV_SEQ_Z3D;
  c3d_pos[ivar] = ivar * siz_volume;
  sprintf(c3d_name[ivar],"%s","z");

  // coord name
  sprintf(coord_name[0],"%s","i");
  sprintf(coord_name[1],"%s","j");
  sprintf(coord_name[2],"%s","k");
  
  // set return values
  *p_c3d = c3d;
  *p_c3d_pos = c3d_pos;
  *p_c3d_name = c3d_name;
  *p_coord_name = coord_name;
  *c3d_num_of_vars = num_grid_vars;
}

void 
gd_curv_init_g3d(
    size_t siz_volume,
    int *g3d_num_of_vars,
    float  **p_g3d,
    size_t **p_g3d_pos,
    char  ***p_g3d_name)
{
  const int num_grid_vars = 10;
  /*
   * 0: jac
   * 1-3: xi_x, xi_y, xi_z
   * 4-6: eta_x, eta_y, eta_z
   * 7-9: zeta_x, zeta_y, zeta_z
   */
  
  // vars
  float *g3d = (float *) fdlib_mem_calloc_1d_float(siz_volume * num_grid_vars,
                                                   0.0,
                                                   "grid_curv_init_g3d");
  if (g3d == NULL) {
      fprintf(stderr,"Error: failed to alloc metric vars\n");
      fflush(stderr);
  }
  
  // position of each var
  size_t *g3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet(num_grid_vars,
                                                         0, 
                                                         "grid_curv_init_g3d");
  
  // name of each var
  char **g3d_name = (char **) fdlib_mem_malloc_2l_char(num_grid_vars,
                                                       FD_MAX_STRLEN,
                                                       "grid_curv_init_g3d");
  
  // set value
  int ivar; 
  ivar = GD_CURV_SEQ_JAC;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","jac");
  
  ivar = GD_CURV_SEQ_XIX;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","xi_x");
  
  ivar = GD_CURV_SEQ_XIY;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","xi_y");
  
  ivar = GD_CURV_SEQ_XIZ;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","xi_z");
  
  ivar = GD_CURV_SEQ_ETX;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","eta_x");
  
  ivar = GD_CURV_SEQ_ETY;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","eta_y");
  
  ivar = GD_CURV_SEQ_ETZ;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","eta_z");
  
  ivar = GD_CURV_SEQ_ZTX;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","zeta_x");
  
  ivar = GD_CURV_SEQ_ZTY;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","zeta_y");
  
  ivar = GD_CURV_SEQ_ZTZ;
  g3d_pos[ivar] = ivar * siz_volume;
  sprintf(g3d_name[ivar],"%s","zeta_z");
  
  // set return values
  *p_g3d = g3d;
  *p_g3d_pos = g3d_pos;
  *p_g3d_name = g3d_name;
  *g3d_num_of_vars = num_grid_vars;
}

//
// need to change to use fdlib_math.c
//
void
gd_curv_cal_metric(
    float *restrict c3d,
    float *restrict g3d,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    int nx, int ny, int nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    int fd_len, int *restrict fd_indx, float *restrict fd_coef)
{
  float x_xi, x_et, x_zt;
  float y_xi, y_et, y_zt;
  float z_xi, z_et, z_zt;
  float jac;
  float vec1[3], vec2[3], vec3[3], vecg[3];

  int n_fd;

  // point to each var
  float *restrict x3d  = c3d + GD_CURV_SEQ_X3D * siz_volume;
  float *restrict y3d  = c3d + GD_CURV_SEQ_Y3D * siz_volume;
  float *restrict z3d  = c3d + GD_CURV_SEQ_Z3D * siz_volume;
  float *restrict jac3d= g3d + GD_CURV_SEQ_JAC * siz_volume;
  float *restrict xi_x = g3d + GD_CURV_SEQ_XIX * siz_volume;
  float *restrict xi_y = g3d + GD_CURV_SEQ_XIY * siz_volume;
  float *restrict xi_z = g3d + GD_CURV_SEQ_XIZ * siz_volume;
  float *restrict et_x = g3d + GD_CURV_SEQ_ETX * siz_volume;
  float *restrict et_y = g3d + GD_CURV_SEQ_ETY * siz_volume;
  float *restrict et_z = g3d + GD_CURV_SEQ_ETZ * siz_volume;
  float *restrict zt_x = g3d + GD_CURV_SEQ_ZTX * siz_volume;
  float *restrict zt_y = g3d + GD_CURV_SEQ_ZTY * siz_volume;
  float *restrict zt_z = g3d + GD_CURV_SEQ_ZTZ * siz_volume;

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
    float *restrict c3d,
    size_t siz_volume,
    int nx, float dx, float x0,
    int ny, float dy, float y0,
    int nz, float dz, float z0)
{
  float *x3d = c3d + GD_CURV_SEQ_X3D * siz_volume;
  float *y3d = c3d + GD_CURV_SEQ_Y3D * siz_volume;
  float *z3d = c3d + GD_CURV_SEQ_Z3D * siz_volume;

  size_t iptr = 0;
  for (size_t k=0; k<nz; k++)
  {
    for (size_t j=0; j<ny; j++)
    {
      for (size_t i=0; i<nx; i++)
      {
        x3d[iptr] = x0 + i * dx;
        y3d[iptr] = y0 + j * dy;
        z3d[iptr] = z0 + k * dz;

        iptr++;
      }
    }
  }
}

//
// input/output
//
void
gd_curv_metric_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords)
{
  // construct file name
  char in_file[FD_MAX_STRLEN];
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
  char in_file[FD_MAX_STRLEN];
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
gd_curv_coord_export(float  *restrict c3d,
                     size_t *restrict c3d_pos,
                     char  **restrict c3d_name,
                     int number_of_vars,
                     int  nx,
                     int  ny,
                     int  nz,
                     int  ni1,
                     int  nj1,
                     int  nk1,
                     int  ni,
                     int  nj,
                     int  nk,
                     int  gni1,
                     int  gnj1,
                     int  gnk1,
                     char *fname_coords,
                     char *output_dir)
{
  // construct file name
  char ou_file[FD_MAX_STRLEN];
  sprintf(ou_file, "%s/coord_%s.nc", output_dir, fname_coords);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[FD_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i"  , nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "j" , ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, c3d_name[ivar], NC_FLOAT, FD_NDIM, dimid, &varid[ivar]);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nj1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,FD_NDIM,l_start);

  int g_start[] = { gni1, gnj1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,FD_NDIM,g_start);

  int l_count[] = { ni, nj, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,FD_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = c3d + c3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
}

void
gd_curv_metric_export(float  *restrict g3d,
                      size_t *restrict g3d_pos,
                      char  **restrict g3d_name,
                      int number_of_vars,
                      int  nx,
                      int  ny,
                      int  nz,
                      int  ni1,
                      int  nj1,
                      int  nk1,
                      int  ni,
                      int  nj,
                      int  nk,
                      int  gni1,
                      int  gnj1,
                      int  gnk1,
                      char *fname_coords,
                      char *output_dir)
{
  // construct file name
  char ou_file[FD_MAX_STRLEN];
  sprintf(ou_file, "%s/metric_%s.nc", output_dir, fname_coords);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[FD_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i"  , nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "j" , ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, g3d_name[ivar], NC_FLOAT, FD_NDIM, dimid, &varid[ivar]);
  }

  // attribute: index in output snapshot, index w ghost in thread
  int l_start[] = { ni1, nj1, nk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                   NC_INT,FD_NDIM,l_start);

  int g_start[] = { gni1, gnj1, gnk1 };
  nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                   NC_INT,FD_NDIM,g_start);

  int l_count[] = { ni, nj, nk };
  nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                   NC_INT,FD_NDIM,l_count);

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = g3d + g3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
}



// code by SunWenliang
int gd_curv_gen_layer(float*restrict c3d, int nLayers,
                      int* NCellPerlay, int* VmapSpacingIsequal,
                      int nx, int ni, int gni1, int fdx_nghosts, int n_total_grid_x,
                      int ny, int nj, int gnj1, int fdy_nghosts, int n_total_grid_y,
                      int nz, int nk, int gnk1, int fdz_nghosts, int n_total_grid_z)
{
  size_t siz_volume       = nx  * ny  * nz ; 
  size_t siz_volume_layer = n_total_grid_x  * n_total_grid_y  * (nLayers+1) ; 
  size_t siz_volume_layer_ch = ( n_total_grid_x + 2*fdx_nghosts )  * ( n_total_grid_y + 2*fdy_nghosts ) * (nLayers+1) ; 

  float * layer3d      = NULL;
  float * layer3d_load = NULL;
  float * layer3d_load_ch = NULL;
  layer3d         = ( float * ) malloc( sizeof( float ) * nx * ny * (nLayers+1) * 3 );
  layer3d_load    = ( float * ) malloc( sizeof( float ) * siz_volume_layer      * 3 );
  layer3d_load_ch = ( float * ) malloc( sizeof( float ) * siz_volume_layer_ch   * 3 );

  float dd;
  float *x3d = c3d + GD_CURV_SEQ_X3D * nx*ny*nz + nx*ny*fdz_nghosts;
  float *y3d = c3d + GD_CURV_SEQ_Y3D * nx*ny*nz + nx*ny*fdz_nghosts;
  float *z3d = c3d + GD_CURV_SEQ_Z3D * nx*ny*nz + nx*ny*fdz_nghosts;
  float *x3d_ch = c3d + GD_CURV_SEQ_X3D * nx*ny*nz;
  float *y3d_ch = c3d + GD_CURV_SEQ_Y3D * nx*ny*nz;
  float *z3d_ch = c3d + GD_CURV_SEQ_Z3D * nx*ny*nz;


  float *xlayer3d = layer3d + GD_CURV_SEQ_X3D * nx*ny*(nLayers+1);
  float *ylayer3d = layer3d + GD_CURV_SEQ_Y3D * nx*ny*(nLayers+1);
  float *zlayer3d = layer3d + GD_CURV_SEQ_Z3D * nx*ny*(nLayers+1);
  float *xlayer3d_load = layer3d_load + GD_CURV_SEQ_X3D * siz_volume_layer;
  float *ylayer3d_load = layer3d_load + GD_CURV_SEQ_Y3D * siz_volume_layer;
  float *zlayer3d_load = layer3d_load + GD_CURV_SEQ_Z3D * siz_volume_layer;
  float *xlayer3d_load_ch = layer3d_load_ch + GD_CURV_SEQ_X3D * siz_volume_layer_ch;
  float *ylayer3d_load_ch = layer3d_load_ch + GD_CURV_SEQ_Y3D * siz_volume_layer_ch;
  float *zlayer3d_load_ch = layer3d_load_ch + GD_CURV_SEQ_Z3D * siz_volume_layer_ch;
  
  float * zlayerpart  = NULL;
  float * ylayerpart  = NULL;
  float * xlayerpart  = NULL;
  zlayerpart = ( float * ) malloc( sizeof( float ) * (nLayers+1) );
  ylayerpart = ( float * ) malloc( sizeof( float ) * (nLayers+1) );
  xlayerpart = ( float * ) malloc( sizeof( float ) * (nLayers+1) );
  float * z3dpart  = NULL;
  float * y3dpart  = NULL;
  float * x3dpart  = NULL;
  z3dpart = ( float * ) malloc( sizeof( float ) * nk );
  y3dpart = ( float * ) malloc( sizeof( float ) * nk );
  x3dpart = ( float * ) malloc( sizeof( float ) * nk );
////////////////////////////////////////////////////////////////////////////////////////
  FILE * fp;
  fp = fopen("./layer3d_hill.dat","rb");
  if (!fp){
    printf("Failed to open input file of interface!!");
  }
  fread( layer3d_load, sizeof(float),  siz_volume_layer * 3, fp );
  fclose( fp );
///////////////////////////////////////////////////////////////////////////////////////
 
  size_t iptr1 = 0;
  size_t iptr2 = 0;
  size_t iptr3 = 0;
  float dd_ch;

  for ( int k=0; k<nLayers+1; k++ )
  {
    for (int j=0; j<n_total_grid_y; j++)
    {
      for (int i=0; i<n_total_grid_x; i++)
      {
        iptr1 = INDEX( i + fdx_nghosts, j + fdy_nghosts, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr2 = INDEX( i              , j              , k, n_total_grid_x                , n_total_grid_y                 );
        xlayer3d_load_ch[ iptr1 ] = xlayer3d_load[ iptr2 ];
        ylayer3d_load_ch[ iptr1 ] = ylayer3d_load[ iptr2 ];
        zlayer3d_load_ch[ iptr1 ] = zlayer3d_load[ iptr2 ];
      }
    }
  }

  for ( int k=0; k<nLayers+1; k++ )
  {
    for (int j=0; j<n_total_grid_y+ 2*fdy_nghosts; j++)
    {
      for (int i=0; i<fdx_nghosts; i++)
      {
        iptr1 = INDEX( i            , j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr2 = INDEX( fdx_nghosts  , j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr3 = INDEX( fdx_nghosts+1, j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        dd_ch = xlayer3d_load_ch[ iptr3 ] - xlayer3d_load_ch[ iptr2 ];
        xlayer3d_load_ch[ iptr1 ]  = xlayer3d_load_ch[ iptr2 ] - dd_ch * ( fdx_nghosts - i);
        ylayer3d_load_ch[ iptr1 ]  = ylayer3d_load_ch[ iptr2 ];
        zlayer3d_load_ch[ iptr1 ]  = zlayer3d_load_ch[ iptr2 ];

        iptr1 = INDEX( n_total_grid_x + 2*fdx_nghosts -1 - i , j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr2 = INDEX( n_total_grid_x + 1*fdx_nghosts -1     , j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr3 = INDEX( n_total_grid_x + 1*fdx_nghosts -2     , j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        dd_ch = xlayer3d_load_ch[ iptr2 ] - xlayer3d_load_ch[ iptr3 ];
        xlayer3d_load_ch[ iptr1 ]  = xlayer3d_load_ch[ iptr2 ] + dd_ch * ( fdx_nghosts - i);
        ylayer3d_load_ch[ iptr1 ]  = ylayer3d_load_ch[ iptr2 ];
        zlayer3d_load_ch[ iptr1 ]  = zlayer3d_load_ch[ iptr2 ];
      }
    }
  }

  for ( int k=0; k<nLayers+1; k++ )
  {
    for (int j=0; j<fdy_nghosts; j++)
    {
      for ( int i=0; i<n_total_grid_x + 2*fdx_nghosts; i++)
      {
        iptr1 = INDEX( i, j          , k, n_total_grid_x + 2*fdx_nghosts   , n_total_grid_y + 2*fdy_nghosts );
        iptr2 = INDEX( i, fdy_nghosts, k, n_total_grid_x + 2*fdx_nghosts   , n_total_grid_y + 2*fdy_nghosts );
        iptr3 = INDEX( i, fdy_nghosts, k, n_total_grid_x + 2*fdx_nghosts + 1, n_total_grid_y + 2*fdy_nghosts );
        dd_ch = ylayer3d_load_ch[ iptr3 ] - ylayer3d_load_ch[ iptr2 ];
        xlayer3d_load_ch[ iptr1 ]  = xlayer3d_load_ch[ iptr2 ];
        ylayer3d_load_ch[ iptr1 ]  = ylayer3d_load_ch[ iptr2 ] - dd_ch * ( fdy_nghosts - j);
        zlayer3d_load_ch[ iptr1 ]  = zlayer3d_load_ch[ iptr2 ];

        iptr1 = INDEX( i, n_total_grid_y + 2*fdy_nghosts - 1 - j, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr2 = INDEX( i, n_total_grid_y + 1*fdy_nghosts - 1    , k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        iptr3 = INDEX( i, n_total_grid_y + 1*fdy_nghosts - 2    , k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        dd_ch = ylayer3d_load_ch[ iptr2 ] - ylayer3d_load_ch[ iptr3 ];
        xlayer3d_load_ch[ iptr1 ]  = xlayer3d_load_ch[ iptr2 ];
        ylayer3d_load_ch[ iptr1 ]  = ylayer3d_load_ch[ iptr2 ] + dd_ch * ( fdy_nghosts - j);
        zlayer3d_load_ch[ iptr1 ]  = zlayer3d_load_ch[ iptr2 ];
      }
    }
  }

  for (int k=0; k<( nLayers + 1 ); k++)
  {
    for (int j=0; j<ny; j++)
    {
      for (int i=0; i<nx; i++)
      {
        iptr1 = INDEX( i       , j       , k, nx                            , ny                             );
        iptr2 = INDEX( i + gni1, j + gnj1, k, n_total_grid_x + 2*fdx_nghosts, n_total_grid_y + 2*fdy_nghosts );
        xlayer3d[ iptr1 ]  = xlayer3d_load_ch[ iptr2 ];
        ylayer3d[ iptr1 ]  = ylayer3d_load_ch[ iptr2 ];
        zlayer3d[ iptr1 ]  = zlayer3d_load_ch[ iptr2 ];
      }
    }
  }

  dd = zlayer3d[INDEX( 0, 0, 1, nx, ny )] - zlayer3d[INDEX( 0, 0, 0, nx, ny )];
  for ( int i = 0; i < nx; i ++ )
  {
    for ( int j = 0; j < ny; j ++ )
    {
      gd_grid_z_interp( i, j, z3d, zlayer3d, NCellPerlay, VmapSpacingIsequal, nLayers, nx, ny );
      float abszfirstpts = fabs( z3d[ INDEX( i, j, 0, nx, ny ) ] );
      if ( dd > 0 )
      {
        for ( int k = 0; k < nk; k ++ )
        {
          z3dpart[ k ] = abszfirstpts + z3d[ INDEX( i, j, k, nx, ny  ) ] + 1;
        }
        for ( int k = 0; k < nLayers+1; k ++ )
        {
          xlayerpart[ k ] = xlayer3d[ INDEX( i, j, k, nx, ny ) ] ;
          ylayerpart[ k ] = ylayer3d[ INDEX( i, j, k, nx, ny ) ] ;
          zlayerpart[ k ] = abszfirstpts + zlayer3d[ INDEX( i, j, k, nx, ny ) ] + 1;
        }
      }
      else
      {
        for ( int k = 0; k < nk; k ++ )
        {
          z3dpart[ k ] = abszfirstpts - z3d[ INDEX( i, j, k, nx, ny ) ] + 1;
        }
        for ( int k = 0; k < nLayers+1; k ++ )
        {
          xlayerpart[ k ] = xlayer3d[ INDEX( i, j, k, nx, ny ) ] ;
          ylayerpart[ k ] = ylayer3d[ INDEX( i, j, k, nx, ny ) ] ;
          zlayerpart[ k ] =  abszfirstpts - zlayer3d[ INDEX( i, j, k, nx, ny ) ] + 1;
        }
      }
      if (i*j==1)
      {
        for (int ii = 0; ii<60; ii++ )
        {
           printf( "%d: %f -> %f\n",  ii, z3d[ INDEX( i, j, ii, nx, ny )], xlayerpart[ ii ]);
        }
      }
      
      gd_SPL( nLayers+1, zlayerpart, xlayerpart, nk, z3dpart, x3dpart);
      gd_SPL( nLayers+1, zlayerpart, ylayerpart, nk, z3dpart, y3dpart);
      for ( int k = 0; k < nk; k ++)
      {
        x3d[ INDEX( i, j, k, nx, ny ) ] = x3dpart[ k ];
        y3d[ INDEX( i, j, k, nx, ny ) ] = y3dpart[ k ];
      }
    }
  }

  for ( int k=0; k<fdz_nghosts; k++ )
  {
    for (int j=0; j<ny; j++)
    {
      for ( int i=0; i<nx; i++)
      {
        iptr1 = INDEX( i, j, k              , nx, ny );
        iptr2 = INDEX( i, j, fdz_nghosts    , nx, ny );
        iptr3 = INDEX( i, j, fdz_nghosts + 1, nx, ny );
        dd_ch = z3d_ch[ iptr3 ] - z3d_ch[ iptr2 ];
        x3d_ch[ iptr1 ]  = x3d_ch[ iptr2 ];
        y3d_ch[ iptr1 ]  = y3d_ch[ iptr2 ];
        z3d_ch[ iptr1 ]  = z3d_ch[ iptr2 ] - dd_ch * ( fdz_nghosts - k ) ;

        iptr1 = INDEX( i, j, nz - 1 - k              , nx, ny );
        iptr2 = INDEX( i, j, nz - 1 -fdz_nghosts     , nx, ny );
        iptr3 = INDEX( i, j, nz - 1 -fdz_nghosts - 1 , nx, ny );
        dd_ch = z3d_ch[ iptr2 ] - z3d_ch[ iptr3 ];
        x3d_ch[ iptr1 ]  = x3d_ch[ iptr2 ];
        y3d_ch[ iptr1 ]  = y3d_ch[ iptr2 ];
        z3d_ch[ iptr1 ]  = z3d_ch[ iptr2 ] + dd_ch * ( fdz_nghosts - k );
      }
    }
  }
  free( layer3d );
  free( layer3d_load );
  free( layer3d_load_ch );
  free( zlayerpart );
  free( ylayerpart );
  free( xlayerpart );
  free( z3dpart );
  free( y3dpart );
  free( x3dpart );
  return 0;
}

int gd_grid_z_interp(int xi, int yi, float* z3d, float* zlayer3d, int* NCellPerlay, int* VmapSpacingIsequal, int nLayers, int nx, int ny )
{
  int i = 0, ii = 0, j = 0, N1 = 0, N2 = 0;
  float distance = 0.0, dmidpoint = 0.0;
  int sumNCellPerlay = 0;
  float * LayerDz = NULL;
  float * range1  = NULL;
  float * range2  = NULL;
  LayerDz = ( float * ) malloc( sizeof( float ) * nLayers );
  for( i = 0; i < nLayers; i ++ ){
    LayerDz[i] = (zlayer3d[INDEX( xi, yi, i+1, nx, ny )] - zlayer3d[INDEX( xi, yi, i, nx, ny )]) / NCellPerlay[i];
  }
  for( i = 0; i < nLayers; i ++ )
  {
    for (ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i]; ii ++ )
    {
      z3d[INDEX( xi, yi, ii, nx, ny )] = zlayer3d[INDEX( xi, yi, i, nx, ny )] + LayerDz[i]/2 + ( ii - sumNCellPerlay ) * LayerDz[i] ;
    }
    if (VmapSpacingIsequal[i] < 1 && i > 0 && i < nLayers-1) // The grid spacing is equal)
    {
      range1 = ( float * ) malloc( sizeof( float ) * (NCellPerlay[i] +1) );
      range2 = ( float * ) malloc( sizeof( float ) * (NCellPerlay[i] +1) );
      range1[0] = zlayer3d[INDEX( xi, yi, i, nx, ny )];
      range2[0] = zlayer3d[INDEX( xi, yi, i, nx, ny )];
      if ( (LayerDz[i] - LayerDz[i-1]) * (LayerDz[i+1] - LayerDz[i]) > 0 )
      {
        N1 = floor(NCellPerlay[i]/2);
        N2 = NCellPerlay[i] - N1;
        distance = zlayer3d[INDEX( xi, yi, i+1, nx, ny )] - zlayer3d[INDEX( xi, yi, i, nx, ny )];
        dmidpoint = ( distance*2 - N1 * LayerDz[i-1] - N2*LayerDz[i+1] ) / NCellPerlay[i];
        
        for ( j = 0; j < N1; j ++)
        {
          range1[j+1] = range1[j] + LayerDz[i-1] + j*(dmidpoint - LayerDz[i-1])/(N1-1);
        }
        for ( j = N1; j < N1+N2; j ++)
        {
          range1[j+1] = range1[j] + dmidpoint + (j-N1)*(LayerDz[i+1] - dmidpoint)/(N2-1);
        }
        // Vertical smooth
        range2[N1+N2] = range1[N1+N2];
        for (j = 1; j < N1+N2; j ++ )
        {
          range2[j] = ( range1[j-1] + range1[j] + range1[j+1] ) / 3;
        }
        for (j = 1; j < N1+N2; j ++ )
        {
          range1[j] = ( range2[j-1] + range2[j] + range2[j+1] ) / 3;
        }
        for (ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i]; ii ++ )
        {
          z3d[INDEX( xi, yi, ii, nx, ny )] = ( range1[ii-sumNCellPerlay] + range1[ii+1-sumNCellPerlay] ) / 2;
        }
      }
    }
    sumNCellPerlay += NCellPerlay[i];
  }
  free( LayerDz );
  free( range1 );
  free( range2 );
  return 0;
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
    /* Forward elimination */
    for (i = 1; i < n; ++i)
    {
      t = d[i - 1] / b[i - 1];
      b[i] = b[i] - t * d[i - 1];
      c[i] = c[i] - t * c[i - 1];
    }
    /* Back substitution */
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

void gd_SPL(int n, float *x, float *y, int ni, float *xi, float *yi)
{
  float *b, *c, *d;
  int iflag=0, last=0, i=0;
  b = (float *)malloc(sizeof(float) * n);
  c = (float *)malloc(sizeof(float) * n);
  d = (float *)malloc(sizeof(float) * n);
  if (!d) { printf("no enough memory for b,c,d\n"); }
  else {
    gd_SPLine(n, 0, 0, 0, 0, x, y, b, c, d, &iflag);
    for (i = 0; i<ni; i++)
      yi[i] = gd_seval(ni, xi[i], n, x, y, b, c, d, &last);
    free(b);
    free(c);
    free(d);
  };
}

void getMinMaxCoor(float *x3d, float *y3d, float *z3d,
                   size_t siz_volume,
                   float *xmin, float *ymin, float *zmin,
                   float *xmax, float *ymax, float *zmax) 
{
    for (size_t i = 0; i < siz_volume; i++){
        *xmin = *xmin < x3d[i] ? *xmin : x3d[i];
        *xmax = *xmax > x3d[i] ? *xmax : x3d[i];
        *ymin = *ymin < y3d[i] ? *ymin : y3d[i];
        *ymax = *ymax > y3d[i] ? *ymax : y3d[i];
        *zmin = *zmin < z3d[i] ? *zmin : z3d[i];
        *zmax = *zmax > z3d[i] ? *zmax : z3d[i];
    }
}



