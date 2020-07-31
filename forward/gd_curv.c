/*
********************************************************************************
* Curve grid metric calculation using MacCormack scheme                        *
********************************************************************************
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "netcdf.h"
#include "gd_curv.h"
#include "fdlib_math.h"

//int grid_curv_init_vars(size_t siz_volume, int *number_of_vars,
int gd_curv_init_vars(size_t siz_volume, int *number_of_vars,
    float **p_g3d, size_t **p_g3d_pos, char ***p_g3d_name)
{
    int ivar;

    *number_of_vars = 12;
    /*
     * 0-2: x3d, y3d, z3d
     * 3: jac
     * 4-6: xi_x, xi_y, xi_z
     * 7-9: eta_x, eta_y, eta_z
     * 10-12: zeta_x, zeta_y, zeta_z
     */

    // vars
    *p_g3d = (float *) fdlib_mem_calloc_1d_float( 
                 siz_volume * (*number_of_vars), 0.0, "grid_curv_init_vars");

    if (*p_g3d == NULL) {
        fprintf(stderr,"Error: failed to alloc grid vars\n");
        fflush(stderr);
        ierr = -1;
    }

    // position of each var
    *p_g3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet( 
                 *number_of_vars, 0, "grid_curv_init_vars");

    // name of each var
    *p_g3d_name = (char **) fdlib_mem_malloc_2d_char( 
                 *number_of_vars, FDSYS_MAX_STR_LEN, "grid_curv_init_vars");

    // set value
    ivar = GRID_CURV_SEQ_X3D;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"x");

    ivar = GRID_CURV_SEQ_Y3D;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"y");

    ivar = GRID_CURV_SEQ_Z3D;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"z");

    ivar = GRID_CURV_SEQ_JAC;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"jac");

    ivar = GRID_CURV_SEQ_XIX;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"xi_x");

    ivar = GRID_CURV_SEQ_XIY;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"xi_y");

    ivar = GRID_CURV_SEQ_XIZ;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"xi_z");

    ivar = GRID_CURV_SEQ_ETX;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"eta_x");

    ivar = GRID_CURV_SEQ_ETY;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"eta_y");

    ivar = GRID_CURV_SEQ_ETZ;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"eta_z");

    ivar = GRID_CURV_SEQ_ZTX;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"zeta_x");

    ivar = GRID_CURV_SEQ_ZTY;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"zeta_y");

    ivar = GRID_CURV_SEQ_ZTZ;
    (*p_g3d)[ivar] = ivar * siz_volume;
    strcpy((*p_g3d_name)[ivar],"zeta_z");

    return 0;
}

//
//
//
int gd_curv_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_metric_dir, int *myid3)
{
    int ierr = 0;
    char in_file[FDSYS_MAX_STR_LEN];

    int ncid, varid;

    // construct file name
    sprintf(in_file, "%s/metric_mpi%02d%02d%02d.nc", in_metric_dir, myid3[0], myid3[1], myid3[2]);

    // read in nc
    ierr = nc_open(in_file, NC_NOWRITE, &ncid);
    if (ierr != NC_NOERR) errh(ierr);

    for (int ivar=0; ivar<number_of_vars; ivar++) {
        ierr = nc_inq_varid(ncid, g3d_name[ivar], &varid);
        if (ierr != NC_NOERR) errh(ierr);

        ierr = nc_get_vara_float(ncid,varid,g3d+g3d_pos[ivar]);
        if (ierr != NC_NOERR) errh(ierr);
    }
    
    // close file
    ierr = nc_close(ncid);
    if (ierr != NC_NOERR) errh(ierr);

    return ierr;
}

//
// need to change to use fdlib_math.c
//
int gd_curv_cal_metric(float *restrict g3d, size_t nx, size_t ny, size_t nz, int op_half_len,
        int fd_length, 
        size_t *restrict fdx_shift, float *restrict fdx_coef,
        size_t *restrict fdy_shift, float *restrict fdy_coef,
        size_t *restrict fdz_shift, float *restrict fdz_coef)
{
  size_t ni1, ni2, nj1, nj2, nk1, nk2, siz_volume;
  size_t i, j, k, iptr;
  int ivar;

  float x_xi, x_et, x_zt;
  float y_xi, y_et, y_zt;
  float z_xi, z_et, z_zt;
  float jac;
  float vec1[3], vec2[3], vec3[3], vecg[3];

  float *restrict jac; 
  float *restrict xi_x, float *restrict xi_y, float *restrict xi_z;
  float *restrict eta_x, float *restrict eta_y, float *restrict eta_z;
  float *restrict zeta_x, float *restrict zeta_y, float *restrict zeta_z;

  siz_volume = nx * ny * nz;

  // point to each var
  jac  = g3d + GRID_CURV_SEQ_JAC * siz_volume;
  xi_x = g3d + GRID_CURV_SEQ_XIX * siz_volume;
  xi_y = g3d + GRID_CURV_SEQ_XIY * siz_volume;
  xi_z = g3d + GRID_CURV_SEQ_XIZ * siz_volume;
  et_x = g3d + GRID_CURV_SEQ_ETX * siz_volume;
  et_y = g3d + GRID_CURV_SEQ_ETY * siz_volume;
  et_z = g3d + GRID_CURV_SEQ_ETZ * siz_volume;
  zt_x = g3d + GRID_CURV_SEQ_ZTX * siz_volume;
  zt_y = g3d + GRID_CURV_SEQ_ZTY * siz_volume;
  zt_z = g3d + GRID_CURV_SEQ_ZTZ * siz_volume;

  ni1 =  0 + op_half_len;
  ni2 = nx - op_half_len;
  nj1 =  0 + op_half_len;
  nj2 = ny - op_half_len;
  nk1 =  0 + op_half_len;
  nk2 = nz - op_half_len;

  for (k = nk1; k <= nk2; k++){
    for (j = nj1; j <= nj2; j++) {
      for (i = ni1; i <= ni2; i++)
      {
        iptr = i + j * nx + k * nx * ny;

        x_xi = 0.0; x_et = 0.0; x_zt = 0.0;
        y_xi = 0.0; y_et = 0.0; y_zt = 0.0;
        z_xi = 0.0; z_et = 0.0; z_zt = 0.0;

        M_CAL_FD_VAL(x_xi, x3d, iptr, fd_length, fdx_shift, fdx_coef, n);
        M_CAL_FD_VAL(y_xi, y3d, iptr, fd_length, fdx_shift, fdx_coef, n);
        M_CAL_FD_VAL(z_xi, z3d, iptr, fd_length, fdx_shift, fdx_coef, n);

        M_CAL_FD_VAL(x_ei, x3d, iptr, fd_length, fdy_shift, fdy_coef, n);
        M_CAL_FD_VAL(y_ei, y3d, iptr, fd_length, fdy_shift, fdy_coef, n);
        M_CAL_FD_VAL(z_ei, z3d, iptr, fd_length, fdy_shift, fdy_coef, n);

        M_CAL_FD_VAL(x_zt, x3d, iptr, fd_length, fdz_shift, fdz_coef, n);
        M_CAL_FD_VAL(y_zt, y3d, iptr, fd_length, fdz_shift, fdz_coef, n);
        M_CAL_FD_VAL(z_zt, z3d, iptr, fd_length, fdz_shift, fdz_coef, n);

        vec1[0] = x_xi; vec1[1] = y_xi; vec1[2] = z_xi;
        vec2[0] = x_et; vec2[1] = y_et; vec2[2] = z_et;
        vec3[0] = x_zt; vec3[1] = y_zt; vec3[2] = z_zt;

        cross_product(vec1, vec2, vecg);
        jac_val = dot_product(vecg, vec3);
        jac[iptr]  = jac_val;

        cross_product(vec2, vec3, vecg);
        xi_x[iptr] = vecg[0] / jac_val;
        xi_y[iptr] = vecg[1] / jac_val;
        xi_z[iptr] = vecg[2] / jac_val;

        cross_product(vec3, vec1, vecg);
        eta_x[iptr] = vecg[0] / jac_val;
        eta_y[iptr] = vecg[1] / jac_val;
        eta_z[iptr] = vecg[2] / jac_val;

        cross_product(vec1, vec2, vecg);
        zeta_x[iptr] = vecg[0] / jac_val;
        zeta_y[iptr] = vecg[1] / jac_val;
        zeta_z[iptr] = vecg[2] / jac_val;
      }
    }
  }

  return;
}

