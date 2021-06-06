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

/*
 * generate grid from a input grid layer file
 * code by SunWenliang  2021/06/02
 */
int gd_curv_gen_layer(char *in_grid_layer_file,
                      int *grid_layer_interp_factor,
                      int *grid_layer_startend,
                      int n_total_grid_x,
                      int n_total_grid_y,
                      int n_total_grid_z,
                      float *restrict c3d,
                      int nx, int ni, int gni1, int fdx_nghosts, 
                      int ny, int nj, int gnj1, int fdy_nghosts, 
                      int nz, int nk, int gnk1, int fdz_nghosts)
{
  int nLayers;
  int n_Interfaces;
  int nx_layers;
  int ny_layers;
  int nx_first;
  int ny_first;
  int nx_interp;
  int ny_interp;
  size_t iptr1 ;
  size_t iptr2 ;
  size_t iptr3 ;

//
// read data from layer file
//

  FILE * fp;
  fp = fopen( in_grid_layer_file, "r");
  if (!fp){
    printf("Failed to open input file of interface!!");
  }

  // read number of interface
  fscanf(fp,"%d", &n_Interfaces);

  nLayers = n_Interfaces - 1;
  int NCellPerlay[nLayers];
  int VmapSpacingIsequal[nLayers];

  // read cells and is_equal of each layer
  for ( int i=0; i<nLayers; i++)
  {
    fscanf(fp,"%d",&NCellPerlay[i]);
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
  layer3d_In   = ( float * ) malloc( sizeof( float ) * siz_volume_layerIn * 3 );
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

  int x_interp_factor = abs(grid_layer_interp_factor[0]);
  int y_interp_factor = abs(grid_layer_interp_factor[1]);
  int z_interp_factor = abs(grid_layer_interp_factor[2]);

  // effective layer nx after downsampling
  if ( grid_layer_interp_factor[0] < 0 ) {
    nx_interp = floor( (nx_layers               +1) / x_interp_factor );
    nx_first  = floor( (grid_layer_startend[0]+1) / x_interp_factor );
  }
  // effective layer nx after upsampling
  else
  { 
    nx_interp = (nx_layers               -1) * x_interp_factor + 1;
    nx_first  = (grid_layer_startend[0]-1) * x_interp_factor + 1;
  }

  if ( grid_layer_interp_factor[1] < 0 ) {
    ny_interp = floor( (ny_layers               +1) / y_interp_factor );
    ny_first  = floor( (grid_layer_startend[1]+1) / y_interp_factor );
  }else{ 
    ny_interp = (ny_layers               -1) * y_interp_factor + 1;
    ny_first  = (grid_layer_startend[1]-1) * y_interp_factor + 1;
  }

// resample based on interp_factor on input grid layer
  size_t siz_volume_layer_interp = nx_interp * ny_interp * (nLayers + 1);
  float *layer3d_interp = NULL;
  //-- attention: check null layer
  layer3d_interp = (float *)malloc(sizeof(float) * siz_volume_layer_interp * 3);
  
  // 
  float   *xx1_lenX = NULL;
  float   *xx2_lenX = NULL;
  float *yy1_x_lenX = NULL;
  float *yy2_x_lenX = NULL;
  float *yy1_y_lenX = NULL;
  float *yy2_y_lenX = NULL;
  float *yy1_z_lenX = NULL;
  float *yy2_z_lenX = NULL;
  float   *xx1_lenY = NULL;
  float   *xx2_lenY = NULL;
  float *yy1_x_lenY = NULL;
  float *yy2_x_lenY = NULL;
  float *yy1_y_lenY = NULL;
  float *yy2_y_lenY = NULL;
  float *yy1_z_lenY = NULL;
  float *yy2_z_lenY = NULL;
    xx1_lenX = (float *)malloc(sizeof(float) * nx_interp);
    xx2_lenX = (float *)malloc(sizeof(float) * nx_layers);
  yy1_x_lenX = (float *)malloc(sizeof(float) * nx_interp);
  yy2_x_lenX = (float *)malloc(sizeof(float) * nx_layers);
  yy1_y_lenX = (float *)malloc(sizeof(float) * nx_interp);
  yy2_y_lenX = (float *)malloc(sizeof(float) * nx_layers);
  yy1_z_lenX = (float *)malloc(sizeof(float) * nx_interp);
  yy2_z_lenX = (float *)malloc(sizeof(float) * nx_layers);
    xx1_lenY = (float *)malloc(sizeof(float) * ny_interp);
    xx2_lenY = (float *)malloc(sizeof(float) * ny_layers);
  yy1_x_lenY = (float *)malloc(sizeof(float) * ny_interp);
  yy2_x_lenY = (float *)malloc(sizeof(float) * ny_layers);
  yy1_y_lenY = (float *)malloc(sizeof(float) * ny_interp);
  yy2_y_lenY = (float *)malloc(sizeof(float) * ny_layers);
  yy1_z_lenY = (float *)malloc(sizeof(float) * ny_interp);
  yy2_z_lenY = (float *)malloc(sizeof(float) * ny_layers);

  // Downsampling in X and Y
  if (grid_layer_interp_factor[0] < 2 && grid_layer_interp_factor[1] < 2 ) 
  {
    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int jj = 0; jj < ny_interp; jj++) {
        for (int ii = 0; ii < nx_interp; ii++) {
          iptr1 = INDEX(ii, jj, kk, nx_interp, ny_interp);
          iptr2 = INDEX(ii * x_interp_factor, jj * y_interp_factor, kk, nx_layers,
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
  else if( grid_layer_interp_factor[0] < 2 && grid_layer_interp_factor[1] >= 2  )  
  {
    for ( int jj1=0; jj1<ny_interp; jj1++) xx1_lenY[jj1] = jj1;
    for ( int jj1=0; jj1<ny_layers; jj1++) xx2_lenY[jj1] = jj1*y_interp_factor;

    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int ii = 0; ii < nx_interp; ii++) {
        for ( int jj1=0; jj1<ny_layers; jj1++)
        {
          yy2_x_lenY[jj1] = layer3d_In[INDEX(ii*x_interp_factor, jj1, kk, nx_layers,
                                       ny_layers)                        ];
          yy2_y_lenY[jj1] = layer3d_In[INDEX(ii*x_interp_factor, jj1, kk, nx_layers,
                                             ny_layers) + siz_volume_layerIn   ];
          yy2_z_lenY[jj1] = layer3d_In[INDEX(ii*x_interp_factor, jj1, kk, nx_layers,
                                             ny_layers) + siz_volume_layerIn*2 ];
        }

        gd_SPL(ny_layers, xx2_lenY, yy2_x_lenY, ny_interp, xx1_lenY, yy1_x_lenY);
        gd_SPL(ny_layers, xx2_lenY, yy2_y_lenY, ny_interp, xx1_lenY, yy1_y_lenY);
        gd_SPL(ny_layers, xx2_lenY, yy2_z_lenY, ny_interp, xx1_lenY, yy1_z_lenY);

        for (int jj = 0; jj < ny_interp; jj++)
        {
          iptr1 = INDEX(ii, jj, kk, nx_interp, ny_interp);
          layer3d_interp[iptr1] = yy1_x_lenY[jj];
          layer3d_interp[iptr1 + siz_volume_layer_interp] = yy1_y_lenY[jj];
          layer3d_interp[iptr1 + siz_volume_layer_interp * 2] = yy1_z_lenY[jj];
        }
      }
    }
  }
  //Interpolating in in X, Downsampling in Y
  else if( grid_layer_interp_factor[0] >=2 && grid_layer_interp_factor[1] < 2  )  
  {
    for ( int ii1=0; ii1<nx_interp; ii1++) xx1_lenX[ii1] = ii1;
    for ( int ii1=0; ii1<nx_layers; ii1++) xx2_lenX[ii1] = ii1*x_interp_factor;

    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int jj = 0; jj < ny_interp; jj++) {
        for ( int ii1=0; ii1<nx_layers; ii1++)
        {
          yy2_x_lenX[ii1] = layer3d_In[INDEX(ii1, jj*y_interp_factor, kk, nx_layers,
                                            ny_layers)                        ];
          yy2_y_lenX[ii1] = layer3d_In[INDEX(ii1, jj*y_interp_factor, kk, nx_layers,
                                            ny_layers) + siz_volume_layerIn   ];
          yy2_z_lenX[ii1] = layer3d_In[INDEX(ii1, jj*y_interp_factor, kk, nx_layers,
                                            ny_layers) + siz_volume_layerIn*2 ];
        }

        gd_SPL(nx_layers, xx2_lenX, yy2_x_lenX, nx_interp, xx1_lenX, yy1_x_lenX);
        gd_SPL(nx_layers, xx2_lenX, yy2_y_lenX, nx_interp, xx1_lenX, yy1_y_lenX);
        gd_SPL(nx_layers, xx2_lenX, yy2_z_lenX, nx_interp, xx1_lenX, yy1_z_lenX);

        for (int ii = 0; ii < nx_interp; ii++)
        {
            iptr1 = INDEX(ii, jj, kk, nx_interp, ny_interp);
            layer3d_interp[iptr1] = yy1_x_lenX[ii];
            layer3d_interp[iptr1 + siz_volume_layer_interp] = yy1_y_lenX[ii];
            layer3d_interp[iptr1 + siz_volume_layer_interp * 2] = yy1_z_lenX[ii];
        }
      }
    }
  }
  // Interpolating in in X and Y
  else if( grid_layer_interp_factor[0] >=2 && grid_layer_interp_factor[1] >= 2  ) 
  {
    // interp x
    for ( int ii1=0; ii1<nx_interp; ii1++) xx1_lenX[ii1] = ii1;
    for ( int ii1=0; ii1<nx_layers; ii1++) xx2_lenX[ii1] = ii1*x_interp_factor;
    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int jj = 0; jj < ny_layers; jj++) {
        for ( int ii1=0; ii1<nx_layers; ii1++)
        {
          yy2_x_lenX[ii1] = layer3d_In[INDEX(ii1, jj, kk, nx_layers, ny_layers)];
          yy2_y_lenX[ii1] = layer3d_In[INDEX(ii1, jj, kk, nx_layers, ny_layers) 
                                       + siz_volume_layerIn   ];
          yy2_z_lenX[ii1] = layer3d_In[INDEX(ii1, jj, kk, nx_layers, ny_layers) 
                                       + siz_volume_layerIn*2 ];
        }

        gd_SPL( nx_layers, xx2_lenX, yy2_x_lenX, nx_interp, xx1_lenX, yy1_x_lenX);
        gd_SPL( nx_layers, xx2_lenX, yy2_y_lenX, nx_interp, xx1_lenX, yy1_y_lenX);
        gd_SPL( nx_layers, xx2_lenX, yy2_z_lenX, nx_interp, xx1_lenX, yy1_z_lenX);

        for (int ii = 0; ii < nx_interp; ii++)
        {
          iptr1 = INDEX( ii, jj*y_interp_factor, kk, nx_interp, ny_interp );
          layer3d_interp[ iptr1                           ] = yy1_x_lenX[ ii ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp  ] = yy1_y_lenX[ ii ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp*2] = yy1_z_lenX[ ii ];
        }
      }
    }

    // interp y
    for ( int jj1=0; jj1<ny_interp; jj1++) xx1_lenY[jj1] = jj1;
    for ( int jj1=0; jj1<ny_layers; jj1++) xx2_lenY[jj1] = jj1*y_interp_factor;
    for (int kk = 0; kk < nLayers+1; kk++) {
      for (int ii = 0; ii < nx_interp; ii++) {
        for ( int jj1=0; jj1<ny_layers; jj1++) {
          yy2_x_lenY[jj1] = layer3d_interp[INDEX(ii, jj1 * y_interp_factor, kk,
                               nx_interp, ny_interp)];
          yy2_y_lenY[jj1] = layer3d_interp[INDEX(ii, jj1 * y_interp_factor, kk,
                               nx_interp, ny_interp) + siz_volume_layer_interp];
          yy2_z_lenY[jj1] = layer3d_interp[INDEX(ii, jj1 * y_interp_factor, kk,
                               nx_interp, ny_interp) + siz_volume_layer_interp * 2];
          }
          gd_SPL( ny_layers, xx2_lenY, yy2_x_lenY, ny_interp, xx1_lenY, yy1_x_lenY);
          gd_SPL( ny_layers, xx2_lenY, yy2_y_lenY, ny_interp, xx1_lenY, yy1_y_lenY);
          gd_SPL( ny_layers, xx2_lenY, yy2_z_lenY, ny_interp, xx1_lenY, yy1_z_lenY);
        for (int jj = 0; jj < ny_interp; jj++) {
          iptr1 = INDEX( ii, jj, kk, nx_interp, ny_interp );
          layer3d_interp[ iptr1                           ] = yy1_x_lenY[ jj ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp  ] = yy1_y_lenY[ jj ];
          layer3d_interp[ iptr1 +siz_volume_layer_interp*2] = yy1_z_lenY[ jj ];
        }
      }
    }  
  }

  free(xx1_lenX);
  free(xx2_lenX);
  free(yy1_x_lenX);
  free(yy2_x_lenX);
  free(yy1_y_lenX);
  free(yy2_y_lenX);
  free(yy1_z_lenX);
  free(yy2_z_lenX);
  free(xx1_lenY);
  free(xx2_lenY);
  free(yy1_x_lenY);
  free(yy2_x_lenY);
  free(yy1_y_lenY);
  free(yy2_y_lenY);
  free(yy1_z_lenY);
  free(yy2_z_lenY);

//
// generate FD grid
//

  size_t siz_volume  = nx  * ny  * nz ; 
  float * layer3d      = NULL;
  layer3d = ( float * ) malloc( sizeof( float ) * nx * ny * (nLayers+1) * 3 );

  // suffix ch means:
  float zdiff_two_interface; // thickness between interfaces
  float zdiff_two_interface_ch;

  float *x3d = c3d + GD_CURV_SEQ_X3D * nx*ny*nz + nx*ny*fdz_nghosts;
  float *y3d = c3d + GD_CURV_SEQ_Y3D * nx*ny*nz + nx*ny*fdz_nghosts;
  float *z3d = c3d + GD_CURV_SEQ_Z3D * nx*ny*nz + nx*ny*fdz_nghosts;
  float *x3d_ch = c3d + GD_CURV_SEQ_X3D * nx*ny*nz;
  float *y3d_ch = c3d + GD_CURV_SEQ_Y3D * nx*ny*nz;
  float *z3d_ch = c3d + GD_CURV_SEQ_Z3D * nx*ny*nz;

  float *xlayer3d = layer3d + GD_CURV_SEQ_X3D * nx*ny*(nLayers+1);
  float *ylayer3d = layer3d + GD_CURV_SEQ_Y3D * nx*ny*(nLayers+1);
  float *zlayer3d = layer3d + GD_CURV_SEQ_Z3D * nx*ny*(nLayers+1);
  
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

  // layer regions of input model covered by this thread
  size_t x_gd_first = nx_first + gni1 - fdx_nghosts;
  size_t y_gd_first = ny_first + gnj1 - fdy_nghosts;
  for (int k = 0; k < (nLayers + 1); k++) {
    for (int j = 0; j < ny; j++) {
      for (int i = 0; i < nx; i++)
      {
        iptr1 = INDEX(i, j, k, nx, ny);
        iptr2 = INDEX(i + x_gd_first, j + y_gd_first, k, nx_interp, ny_interp);
        xlayer3d[iptr1] = layer3d_interp[iptr2];
        ylayer3d[iptr1] = layer3d_interp[iptr2 + siz_volume_layer_interp];
        zlayer3d[iptr1] = layer3d_interp[iptr2 + siz_volume_layer_interp * 2];
      }
    }
  }

  // interp both z and x,y

  zdiff_two_interface = zlayer3d[INDEX( 0, 0, 1, nx, ny )] - zlayer3d[INDEX( 0, 0, 0, nx, ny )];

  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Interpolating the Z coordinate of the grid by interface control point.
      gd_grid_z_interp(i, j, z3d, zlayer3d, NCellPerlay, VmapSpacingIsequal,
                       nLayers, nx, ny);

      // interp x and y
      float abszfirstpts = fabs(z3d[INDEX(i, j, 0, nx, ny)]);

      if ( zdiff_two_interface > 0 )
      {
        for ( int k = 0; k < nk; k ++ )
        {
          z3dpart[ k ] = abszfirstpts + z3d[ INDEX( i, j, k, nx, ny  ) ] + 1;
        }
        for ( int k = 0; k < nLayers+1; k ++ )
        {
          xlayerpart[k] = xlayer3d[INDEX(i, j, k, nx, ny)];
          ylayerpart[k] = ylayer3d[INDEX(i, j, k, nx, ny)];
          zlayerpart[k] = abszfirstpts + zlayer3d[INDEX(i, j, k, nx, ny)] + 1;
        }
      }
      else
      {
        for ( int k = 0; k < nk; k ++ )
        {
          z3dpart[ k ] = abszfirstpts - z3d[ INDEX( i, j, k, nx, ny ) ] + 1;
        }
        for (int k = 0; k < nLayers + 1; k++)
        {
          xlayerpart[ k ] = xlayer3d[ INDEX( i, j, k, nx, ny ) ] ;
          ylayerpart[k] = ylayer3d[INDEX(i, j, k, nx, ny)];
          zlayerpart[k] = abszfirstpts - zlayer3d[INDEX(i, j, k, nx, ny)] + 1;
        }
      }

      // Interpolating the X and Y coordinate of the grid ...
      // by cubic spline interpolation method. 
      gd_SPL( nLayers+1, zlayerpart, xlayerpart, nk, z3dpart, x3dpart);
      gd_SPL( nLayers+1, zlayerpart, ylayerpart, nk, z3dpart, y3dpart);
      for ( int k = 0; k < nk; k ++)
      {
        x3d[ INDEX( i, j, k, nx, ny ) ] = x3dpart[ k ];
        y3d[ INDEX( i, j, k, nx, ny ) ] = y3dpart[ k ];
      }
    }
  }

  // Grids outside the boundary of the upper and lower interfaces. 
  //   
  for ( int k=0; k<fdz_nghosts; k++ )
  {
    for (int j=0; j<ny; j++)
    {
      for ( int i=0; i<nx; i++)
      {
        iptr1 = INDEX(i, j, k, nx, ny);
        iptr2 = INDEX(i, j, fdz_nghosts, nx, ny);
        iptr3 = INDEX(i, j, fdz_nghosts + 1, nx, ny);
        zdiff_two_interface_ch = z3d_ch[iptr3] - z3d_ch[iptr2];
        x3d_ch[iptr1] = x3d_ch[iptr2];
        y3d_ch[iptr1] = y3d_ch[iptr2];
        z3d_ch[iptr1] = z3d_ch[iptr2] - zdiff_two_interface_ch * (fdz_nghosts - k);

        iptr1 = INDEX(i, j, nz - 1 - k, nx, ny);
        iptr2 = INDEX(i, j, nz - 1 - fdz_nghosts, nx, ny);
        iptr3 = INDEX(i, j, nz - 1 - fdz_nghosts - 1, nx, ny);
        zdiff_two_interface_ch = z3d_ch[iptr2] - z3d_ch[iptr3];
        x3d_ch[iptr1] = x3d_ch[iptr2];
        y3d_ch[iptr1] = y3d_ch[iptr2];
        z3d_ch[iptr1] = z3d_ch[iptr2] + zdiff_two_interface_ch * (fdz_nghosts - k);
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
  return 0;
}

//  Interpolating the Z coordinate
int gd_grid_z_interp(int xi, int yi, float* z3d, float* zlayer3d, int* NCellPerlay,
                     int* VmapSpacingIsequal, int nLayers, int nx, int ny )
{
  int i = 0, ii = 0, j = 0, N1 = 0, N2 = 0;
  float distance = 0.0, dmidpoint = 0.0;
  int sumNCellPerlay = 0;
  float * LayerDz = NULL;
  float * range1  = NULL;
  float * range2  = NULL;

  LayerDz = ( float * ) malloc( sizeof( float ) * nLayers );

  for( i = 0; i < nLayers; i ++ ){
    LayerDz[i] = (zlayer3d[INDEX(xi, yi, i + 1, nx, ny)] - 
                  zlayer3d[INDEX(xi, yi, i, nx, ny)]) / NCellPerlay[i];
  }

  for( i = 0; i < nLayers; i ++ )
  {
    for (ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i]; ii ++ )
    {
      z3d[INDEX( xi, yi, ii, nx, ny )] = zlayer3d[INDEX( xi, yi, i, nx, ny )] 
                        + LayerDz[i]/2 + ( ii - sumNCellPerlay ) * LayerDz[i] ;
    }
    // The grid spacing is equal)
    if (VmapSpacingIsequal[i] < 1 && i > 0 && i < nLayers-1) 
    {
      range1 = ( float * ) malloc( sizeof( float ) * (NCellPerlay[i] +1) );
      range2 = ( float * ) malloc( sizeof( float ) * (NCellPerlay[i] +1) );
      range1[0] = zlayer3d[INDEX( xi, yi, i, nx, ny )];
      range2[0] = zlayer3d[INDEX( xi, yi, i, nx, ny )];
      if ( (LayerDz[i] - LayerDz[i-1]) * (LayerDz[i+1] - LayerDz[i]) > 0 )
      {
        N1 = floor(NCellPerlay[i] / 2);
        N2 = NCellPerlay[i] - N1;
        distance = zlayer3d[INDEX(xi, yi, i + 1, nx, ny)] 
                            - zlayer3d[INDEX(xi, yi, i, nx, ny)];
        dmidpoint = (distance * 2 - N1 * LayerDz[i - 1] 
                     - N2 * LayerDz[i + 1]) / NCellPerlay[i];

        for ( j = 0; j < N1; j ++)
        {
          range1[j+1] = range1[j] + LayerDz[i-1] 
                        + j*(dmidpoint - LayerDz[i-1])/(N1-1);
        }
        for ( j = N1; j < N1+N2; j ++)
        {
          range1[j + 1] = range1[j] + dmidpoint 
                          + (j - N1) * (LayerDz[i + 1] - dmidpoint) / (N2 - 1);
        }
        // Vertical smooth
        range2[N1+N2] = range1[N1+N2];
        for (j = 1; j < N1+N2; j ++ ) {
          range2[j] = ( range1[j-1] + range1[j] + range1[j+1] ) / 3;
        }
        for (j = 1; j < N1+N2; j ++ ) {
          range1[j] = ( range2[j-1] + range2[j] + range2[j+1] ) / 3;
        }
        for (ii = sumNCellPerlay; ii < sumNCellPerlay + NCellPerlay[i]; ii ++ ) {
          z3d[INDEX( xi, yi, ii, nx, ny )] = ( range1[ii-sumNCellPerlay] + 
                                              range1[ii+1-sumNCellPerlay] ) / 2;
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
  b = (float *)malloc(sizeof(float) * n);
  c = (float *)malloc(sizeof(float) * n);
  d = (float *)malloc(sizeof(float) * n);
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



