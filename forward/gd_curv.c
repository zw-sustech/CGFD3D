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
    char  ***p_c3d_name)
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
  
  // set value
  int ivar = GD_CURV_SEQ_X3D;
  c3d_pos[ivar] = ivar * siz_volume;
  strcpy(c3d_name[ivar],"x");
  
  ivar = GD_CURV_SEQ_Y3D;
  c3d_pos[ivar] = ivar * siz_volume;
  strcpy(c3d_name[ivar],"y");
  
  ivar = GD_CURV_SEQ_Z3D;
  c3d_pos[ivar] = ivar * siz_volume;
  strcpy(c3d_name[ivar],"z");
  
  // set return values
  *p_c3d = c3d;
  *p_c3d_pos = c3d_pos;
  *p_c3d_name = c3d_name;
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
  strcpy(g3d_name[ivar],"jac");
  
  ivar = GD_CURV_SEQ_XIX;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"xi_x");
  
  ivar = GD_CURV_SEQ_XIY;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"xi_y");
  
  ivar = GD_CURV_SEQ_XIZ;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"xi_z");
  
  ivar = GD_CURV_SEQ_ETX;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"eta_x");
  
  ivar = GD_CURV_SEQ_ETY;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"eta_y");
  
  ivar = GD_CURV_SEQ_ETZ;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"eta_z");
  
  ivar = GD_CURV_SEQ_ZTX;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"zeta_x");
  
  ivar = GD_CURV_SEQ_ZTY;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"zeta_y");
  
  ivar = GD_CURV_SEQ_ZTZ;
  g3d_pos[ivar] = ivar * siz_volume;
  strcpy(g3d_name[ivar],"zeta_z");
  
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
        int number_of_vars, size_t siz_volume, char *in_dir, int *myid2)
{
  // construct file name
  char in_file[FD_MAX_STRLEN];
  sprintf(in_file, "%s/metric_mpi%02d%02d00.nc", in_dir, myid2[0], myid2[1]);
  
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
        int number_of_vars, size_t siz_volume, char *in_dir, int *myid2)
{
  // construct file name
  char in_file[FD_MAX_STRLEN];
  sprintf(in_file, "%s/coord_mpi%02d%02d00.nc", in_dir, myid2[0], myid2[1]);
  
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
                     int *myid2,
                     char *output_dir)
{
  // construct file name
  char ou_file[FD_MAX_STRLEN];
  sprintf(ou_file, "%s/coord_mpi%02d%02d00.nc", output_dir, myid2[0], myid2[1]);
  
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
  ierr = nc_def_dim(ncid, "xi"  , nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "eta" , ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "zeta", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, c3d_name[ivar], NC_FLOAT, FD_NDIM, dimid, &varid[ivar]);
  }

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
                      int *myid2,
                      char *output_dir)
{
  // construct file name
  char ou_file[FD_MAX_STRLEN];
  sprintf(ou_file, "%s/metric_mpi%02d%02d00.nc", output_dir, myid2[0], myid2[1]);
  
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
  ierr = nc_def_dim(ncid, "xi"  , nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "eta" , ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "zeta", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, g3d_name[ivar], NC_FLOAT, FD_NDIM, dimid, &varid[ivar]);
  }

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
