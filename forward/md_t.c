/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"

#include "constants.h"
#include "fdlib_mem.h"
#include "md_t.h"

int
md_init(gdinfo_t *gdinfo, md_t *md, int media_type, int visco_type)
{
  int ierr = 0;

  md->nx   = gdinfo->nx;
  md->ny   = gdinfo->ny;
  md->nz   = gdinfo->nz;

  md->siz_iy   = md->nx;
  md->siz_iz   = md->nx * md->ny;
  md->siz_icmp = md->nx * md->ny * md->nz;

  // media type
  md->medium_type = media_type;
  if (media_type == CONST_MEDIUM_ACOUSTIC_ISO) {
    md->ncmp = 2;
  } else if (media_type == CONST_MEDIUM_ELASTIC_ISO) {
    md->ncmp = 3;
  } else if (media_type == CONST_MEDIUM_ELASTIC_ISO) {
    md->ncmp = 6; // 5 + rho
  } else {
    md->ncmp = 22; // 21 + rho
  }

  // visco
  md->visco_type = visco_type;
  if (visco_type == CONST_VISCO_GRAVES_QS) {
   md->ncmp += 1;
  }

  /*
   * 0: rho
   * 1: lambda
   * 2: mu
   */
  
  // vars
  md->v4d = (float *) fdlib_mem_calloc_1d_float(
                          md->siz_icmp * md->ncmp,
                          0.0, "md_init");
  if (md->v4d == NULL) {
      fprintf(stderr,"Error: failed to alloc medium_el_iso\n");
      fflush(stderr);
  }

  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(md->ncmp,
                                                         0,
                                                         "medium_el_iso_init");

  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(md->ncmp,
                                                       CONST_MAX_STRLEN,
                                                       "medium_el_iso_init");

  // set pos
  for (int icmp=0; icmp < md->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * md->siz_icmp;
  }

  // init
  int icmp = 0;
  sprintf(cmp_name[icmp],"%s","rho");
  md->rho = md->v4d + cmp_pos[icmp];

  // acoustic iso
  if (media_type == CONST_MEDIUM_ACOUSTIC_ISO) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","kappa");
    md->kappa = md->v4d + cmp_pos[icmp];
  }

  // iso
  if (media_type == CONST_MEDIUM_ELASTIC_ISO) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","lambda");
    md->lambda = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","mu");
    md->mu = md->v4d + cmp_pos[icmp];
  }

  // aniso
  if (media_type == CONST_MEDIUM_ELASTIC_ANISO) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c11");
    md->c11 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c12");
    md->c12 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c13");
    md->c13 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c14");
    md->c14 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c15");
    md->c15 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c16");
    md->c16 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c22");
    md->c22 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c23");
    md->c23 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c24");
    md->c24 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c25");
    md->c25 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c26");
    md->c26 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c33");
    md->c33 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c34");
    md->c34 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c35");
    md->c35 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c36");
    md->c36 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c44");
    md->c44 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c45");
    md->c45 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c46");
    md->c46 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c55");
    md->c55 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c56");
    md->c56 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c66");
    md->c66 = md->v4d + cmp_pos[icmp];
  }

  // plus Qs
  if (visco_type == CONST_VISCO_GRAVES_QS) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","Qs");
    md->Qs = md->v4d + cmp_pos[icmp];
  }
  
  // set pointer
  md->cmp_pos  = cmp_pos;
  md->cmp_name = cmp_name;

  return ierr;
}

//
//
//

int
md_import(md_t *md, char *fname_coords, char *in_dir)
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
  
  for (int icmp=0; icmp < md->ncmp; icmp++) {
      ierr = nc_inq_varid(ncid, md->cmp_name[icmp], &varid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }
  
      ierr = nc_get_var_float(ncid,varid,md->v4d + md->cmp_pos[icmp]);
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
md_export(gdinfo_t  *gdinfo,
                 md_t *md,
                 char *fname_coords,
                 char *output_dir)
{
  int ierr = 0;

  size_t *restrict m3d_pos   = md->cmp_pos;
  char  **restrict m3d_name  = md->cmp_name;
  int  number_of_vars = md->ncmp;
  int  nx = md->nx;
  int  ny = md->ny;
  int  nz = md->nz;
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
    float *ptr = md->v4d + m3d_pos[ivar];
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
md_gen_test_ac_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int ny = md->ny;
  int nz = md->nz;
  int siz_line  = md->siz_iy;
  int siz_slice = md->siz_iz;

  float *kappa3d = md->kappa;
  float *rho3d = md->rho;

  for (size_t k=0; k<nz; k++)
  {
    for (size_t j=0; j<ny; j++)
    {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        float Vp=3000.0;
        float rho=1500.0;
        float kappa = Vp*Vp*rho;
        kappa3d[iptr] = kappa;
        rho3d[iptr] = rho;
      }
    }
  }

  return ierr;
}

int
md_gen_test_el_iso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int ny = md->ny;
  int nz = md->nz;
  int siz_line  = md->siz_iy;
  int siz_slice = md->siz_iz;

  float *lam3d = md->lambda;
  float  *mu3d = md->mu;
  float *rho3d = md->rho;

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

int
md_gen_test_Qs(md_t *md, float Qs_freq)
{
  int ierr = 0;

  int nx = md->nx;
  int ny = md->ny;
  int nz = md->nz;
  int siz_line  = md->siz_iy;
  int siz_slice = md->siz_iz;

  md->visco_Qs_freq = Qs_freq;

  float *Qs = md->Qs;

  for (size_t k=0; k<nz; k++)
  {
    for (size_t j=0; j<ny; j++)
    {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        Qs[iptr] = 20;
      }
    }
  }

  return ierr;
}

int
md_gen_test_el_aniso(md_t *md)
{
  int ierr = 0;

  int nx = md->nx;
  int ny = md->ny;
  int nz = md->nz;
  int siz_line  = md->siz_iy;
  int siz_slice = md->siz_iz;

  for (size_t k=0; k<nz; k++)
  {
    for (size_t j=0; j<ny; j++)
    {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;

        float rho=1500.0;

        md->rho[iptr] = rho;

	      md->c11[iptr] = 25.2*1e9;//lam + 2.0f*mu;
	      md->c12[iptr] = 0.0;//lam;
	      md->c13[iptr] = 10.9620*1e9;//lam;
	      md->c14[iptr] = 0.0;
	      md->c15[iptr] = 0.0;
	      md->c16[iptr] = 0.0;
	      md->c22[iptr] = 0.0;//lam + 2.0f*mu;
	      md->c23[iptr] = 0.0;//lam;
	      md->c24[iptr] = 0.0;
	      md->c25[iptr] = 0.0;
	      md->c26[iptr] = 0.0;
	      md->c33[iptr] = 18.0*1e9;//lam + 2.0f*mu;
	      md->c34[iptr] = 0.0;
	      md->c35[iptr] = 0.0;
	      md->c36[iptr] = 0.0;
	      md->c44[iptr] = 0.0;//mu;
	      md->c45[iptr] = 0.0;
	      md->c46[iptr] = 0.0;
	      md->c55[iptr] = 5.12*1e9;//mu;
	      md->c56[iptr] = 0.0;
        md->c66[iptr] = 7.1680*1e9;//mu;

        //-- Vp ~ sqrt(c11/rho) = 4098

        // convert to VTI media
        md->c12[iptr] = md->c11[iptr] - 2.0*md->c66[iptr]; 
	      md->c22[iptr] = md->c11[iptr];
        md->c23[iptr] = md->c13[iptr];
	      md->c44[iptr] = md->c55[iptr]; 
      }
    }
  }

  return ierr;
}

/*
 * convert rho to slowness to reduce number of arithmetic cal
 */

int
md_rho_to_slow(float *restrict rho, size_t siz_volume)
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
