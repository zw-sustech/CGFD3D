/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "md_el_iso.h"

int
md_el_iso_init(gdinfo_t *gdinfo, mdeliso_t *mdeliso)
{
  int ierr = 0;

  mdeliso->nx   = gdinfo->nx;
  mdeliso->ny   = gdinfo->ny;
  mdeliso->nz   = gdinfo->nz;
  mdeliso->ncmp = 3;
  /*
   * 0: rho
   * 1: lambda
   * 2: mu
   */

  mdeliso->siz_iy   = mdeliso->nx;
  mdeliso->siz_iz   = mdeliso->nx * mdeliso->ny;
  mdeliso->siz_icmp = mdeliso->nx * mdeliso->ny * mdeliso->nz;
  
  // vars
  mdeliso->v4d = (float *) fdlib_mem_calloc_1d_float(
                          mdeliso->siz_icmp * mdeliso->ncmp,
                          0.0, "md_el_iso_init");
  if (mdeliso->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc medium_el_iso\n");
      fflush(stderr);
  }

  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(mdeliso->ncmp,
                                                         0,
                                                         "medium_el_iso_init");

  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(mdeliso->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "medium_el_iso_init");

  // set pos
  for (int icmp=0; icmp < mdeliso->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * mdeliso->siz_icmp;
  }

  // init
  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","rho");
  mdeliso->rho = mdeliso->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","lambda");
  mdeliso->lambda = mdeliso->v4d + cmp_pos[icmp];

  icmp += 1;
  sprintf(cmp_name[icmp],"%s","mu");
  mdeliso->mu = mdeliso->v4d + cmp_pos[icmp];
  
  // set pointer
  mdeliso->cmp_pos  = cmp_pos;
  mdeliso->cmp_name = cmp_name;

  return ierr;
}

//
//
//

int
md_el_iso_import(float *restrict m3d, size_t *restrict m3d_pos, char **restrict m3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords)
{
  int ierr = 0;

  char in_file[CONST_MAX_STRLEN];
  
  int ncid, varid;
  
  // construct file name
  sprintf(in_file, "%s/media_%s.nc", in_dir, fname_coords);
  
  // read in nc
  ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }
  
  for (int icmp=0; icmp<number_of_vars; icmp++) {
      ierr = nc_inq_varid(ncid, m3d_name[icmp], &varid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  
      ierr = nc_get_var_float(ncid,varid,m3d+m3d_pos[icmp]);
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

  return ierr;
}

int
md_el_iso_export(gdinfo_t  *gdinfo,
                 mdeliso_t *mdeliso,
                 char *fname_coords,
                 char *output_dir)
{
  int ierr = 0;

  size_t *restrict m3d_pos   = mdeliso->cmp_pos;
  char  **restrict m3d_name  = mdeliso->cmp_name;
  int  number_of_vars = mdeliso->ncmp;
  int  nx = mdeliso->nx;
  int  ny = mdeliso->ny;
  int  nz = mdeliso->nz;
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
  sprintf(ou_file, "%s/media_%s.nc", output_dir, fname_coords);
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
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
    ierr = nc_def_var(ncid, m3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
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
    float *ptr = mdeliso->v4d + m3d_pos[ivar];
    ierr = nc_put_var_float(ncid, varid[ivar],ptr);
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return ierr;
}

/*
 * test
 */

int
md_el_iso_gen_test(mdeliso_t *mdeliso)
{
  int ierr = 0;

  int nx = mdeliso->nx;
  int ny = mdeliso->ny;
  int nz = mdeliso->nz;
  int siz_line  = mdeliso->siz_iy;
  int siz_slice = mdeliso->siz_iz;

  float *lam3d = mdeliso->lambda;
  float  *mu3d = mdeliso->mu;
  float *rho3d = mdeliso->rho;

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

  return ierr;
}

/*
 * convert rho to slowness to reduce number of arithmetic cal
 */

int
md_el_iso_rho_to_slow(float *restrict rho, size_t siz_volume)
{
  int ierr = 0;

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

  return ierr;
}
