/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "fd_t.h"
#include "md_el_iso.h"

void
md_el_iso_init_vars(
    size_t siz_volume,
    int *m3d_num_of_vars, 
    float **p_m3d,
    size_t **p_m3d_pos,
    char ***p_m3d_name)
{
  const int num_medium_vars = 3;
  /*
   * 0: rho
   * 1: lambda
   * 2: mu
   */

  // vars
  float *m3d = (float *) fdlib_mem_calloc_1d_float(siz_volume * num_medium_vars,
                                                   0.0,
                                                   "md_el_iso_init_vars");
  if (m3d == NULL) {
      fprintf(stderr,"Error: failed to alloc medium_el_iso\n");
      fflush(stderr);
  }

  // position of each var
  size_t *m3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet(num_medium_vars,
                                                         0,
                                                         "medium_el_iso_init_vars");

  // name of each var
  char **m3d_name = (char **) fdlib_mem_malloc_2l_char(num_medium_vars,
                                                       FD_MAX_STRLEN,
                                                       "medium_el_iso_init_vars");

  // init
  int ivar = MD_EL_ISO_SEQ_RHO;
  m3d_pos[ivar] = ivar * siz_volume;
  sprintf(m3d_name[ivar],"%s","rho");

  ivar = MD_EL_ISO_SEQ_LAMBDA;
  m3d_pos[ivar] = ivar * siz_volume;
  sprintf(m3d_name[ivar],"%s","lambda");

  ivar = MD_EL_ISO_SEQ_MU;
  m3d_pos[ivar] = ivar * siz_volume;
  sprintf(m3d_name[ivar],"%s","mu");

  // set return values
  *p_m3d = m3d;
  *p_m3d_pos = m3d_pos;
  *p_m3d_name = m3d_name;
  *m3d_num_of_vars = num_medium_vars;
}

//
//
//
void
md_el_iso_import(float *restrict m3d, size_t *restrict m3d_pos, char **restrict m3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords)
{
  char in_file[FD_MAX_STRLEN];
  
  int ncid, varid;
  
  // construct file name
  sprintf(in_file, "%s/media_%s.nc", in_dir, fname_coords);
  
  // read in nc
  int ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
  
  for (int ivar=0; ivar<number_of_vars; ivar++) {
      ierr = nc_inq_varid(ncid, m3d_name[ivar], &varid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  
      ierr = nc_get_var_float(ncid,varid,m3d+m3d_pos[ivar]);
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
md_el_iso_export(float  *restrict m3d,
                 size_t *restrict m3d_pos,
                 char  **restrict m3d_name,
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
  sprintf(ou_file, "%s/media_%s.nc", output_dir, fname_coords);
  
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
  ierr = nc_def_dim(ncid, "i", nx, &dimid[2]);
  ierr = nc_def_dim(ncid, "j", ny, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, m3d_name[ivar], NC_FLOAT, FD_NDIM, dimid, &varid[ivar]);
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
    float *ptr = m3d + m3d_pos[ivar];
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
 * test
 */

void
md_el_iso_gen_test(
    float *restrict m3d,
    float *restrict x3d,
    float *restrict y3d,
    float *restrict z3d,
    int nx,
    int ny,
    int nz,
    size_t siz_line,
    size_t siz_slice,
    size_t siz_volume)
{
  float *lam3d = m3d + MD_EL_ISO_SEQ_LAMBDA * siz_volume;
  float  *mu3d = m3d + MD_EL_ISO_SEQ_MU     * siz_volume;
  float *rho3d = m3d + MD_EL_ISO_SEQ_RHO    * siz_volume;

  for (size_t k=0; k<nz; k++)
  {
    for (size_t j=0; j<ny; j++)
    {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        float Vp=3000.0;
        float Vs=2000.0;
        float rho=1500.0;
        float mu = Vs*Vs*rho;
        float lam = Vp*Vp*rho - 2.0*mu;
        lam3d[iptr] = lam;
         mu3d[iptr] = mu;
        rho3d[iptr] = rho;
      }
    }
  }
}

/*
 * convert rho to slowness to reduce number of arithmetic cal
 */

void
md_el_iso_rho_to_slow(float *restrict m3d, size_t siz_volume)
{
  float *rho = m3d + MD_EL_ISO_SEQ_RHO * siz_volume;

  /*
  for (size_t k=0; k<nx; k++) {
    for (size_t j=0; j<ny; j++) {
      for (size_t i=0; i<nx; i++) {
      }
    }
  }
  */
  for (size_t iptr=0; iptr<siz_volume; iptr++) {
    if (rho[iptr] > 1e-10) {
      rho[iptr] = 1.0 / rho[iptr];
    } else {
      rho[iptr] = 0.0;
    }
  }
}
