/*
 *
 */

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"
#include "netcdf_par.h"

#include "constants.h"
#include "fdlib_mem.h"
#include "md_t.h"
#include "par_t.h"

#ifndef M_NCRET
#define M_NCRET(ierr) {fprintf(stderr,"md nc error: %s\n", nc_strerror(ierr)); exit(1);}
#endif

#ifndef M_NCERR
#define M_NCERR {fprintf(stderr,"md nc error\n"); exit(1);}
#endif

int
md_init(gd_t *gdinfo, md_t *md, int media_type, int visco_type, int nmaxwell)
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
  } else if (media_type == CONST_MEDIUM_ELASTIC_VTI) {
    md->ncmp = 6; // 5 + rho
  } else if (media_type == CONST_MEDIUM_VISCOELASTIC_ISO) {
    md->nmaxwell = nmaxwell;
    md->ncmp = 3 + 2*md->nmaxwell+2;
  } else{
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

  if (visco_type == CONST_VISCO_GMB) 
  {
      md->Ylam = (float **) malloc(md->nmaxwell * sizeof(float *)); 
      md->Ymu  = (float **) malloc(md->nmaxwell * sizeof(float *)); 
      md->wl   = (float *)  malloc(md->nmaxwell * sizeof(float *));
  }

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

  // vti
  if (media_type == CONST_MEDIUM_ELASTIC_VTI) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c11");
    md->c11 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c13");
    md->c13 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c33");
    md->c33 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c55");
    md->c55 = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","c66");
    md->c66 = md->v4d + cmp_pos[icmp];
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

  // vis_iso
  if (media_type == CONST_MEDIUM_VISCOELASTIC_ISO) {
    icmp += 1;
    sprintf(cmp_name[icmp],"%s","lambda");
    md->lambda = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","mu");
    md->mu = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","Qp");
    md->Qp = md->v4d + cmp_pos[icmp];

    icmp += 1;
    sprintf(cmp_name[icmp],"%s","Qs");
    md->Qs = md->v4d + cmp_pos[icmp];

    for(int i=0; i < md->nmaxwell; i++)
    { 
      icmp += 1;
      sprintf(cmp_name[icmp],"%s%d","Ylam",i+1);
      md->Ylam[i] = md->v4d + cmp_pos[icmp];

      icmp += 1;
      sprintf(cmp_name[icmp],"%s%d","Ymu",i+1);
      md->Ymu[i] = md->v4d + cmp_pos[icmp];
    }
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
md_import(gd_t *gd,
          md_t *md,
          int is_parallel_netcdf,
          MPI_Comm comm, 
          char *fname_coords, char *in_dir)
{
  int ierr = 0;

  int  nx = gd->nx;
  int  ny = gd->ny;
  int  nz = gd->nz;

  // construct file name
  char in_file[CONST_MAX_STRLEN];
  
  int ncid, varid;

  // default for seperated nc
  int start_i  = 0;
  int start_j  = 0;
  int start_k  = 0;

  if (is_parallel_netcdf == 1)
  {
    sprintf(in_file, "%s/media.nc", in_dir);

    if (ierr=nc_open_par(in_file, NC_NOWRITE | NC_NETCDF4, 
                      comm, MPI_INFO_NULL, &ncid)) M_NCRET(ierr);

    start_i  = gd->nx1_to_glob_halo0;
    start_j  = gd->ny1_to_glob_halo0;
    start_k  = gd->nz1_to_glob_halo0;
  }
  else
  {
    sprintf(in_file, "%s/media_%s.nc", in_dir, fname_coords);

    if (ierr = nc_open(in_file, NC_NOWRITE, &ncid)) M_NCRET(ierr);
  }
  
  for (int icmp=0; icmp < md->ncmp; icmp++)
  {
    float *ptr = md->v4d + md->cmp_pos[icmp];
    size_t startp[] = { start_k, start_j, start_i };
    size_t countp[] = { nz, ny, nx };

    ierr = nc_inq_varid(ncid, md->cmp_name[icmp], &varid);
    if (ierr != NC_NOERR){
      fprintf(stderr,"md nc inq error: %s\n", nc_strerror(ierr));
      exit(-1);
    }
  
    ierr = nc_get_vara_float(ncid, varid,startp,countp,ptr);
    if (ierr != NC_NOERR){
      fprintf(stderr,"md nc get var error: %s\n", nc_strerror(ierr));
      exit(-1);
    }
  }
  
  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"md nc close error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  return ierr;
}

int
md_export(gd_t  *gdinfo,
          md_t *md,
          int is_parallel_netcdf,
          MPI_Comm comm, 
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
  if(md->medium_type == CONST_MEDIUM_VISCOELASTIC_ISO){
    number_of_vars = 5; //Qp,Qs,not output coef
  }

  // construct file name
  char ou_file[CONST_MAX_STRLEN];

  // default dim size as seperated nc
  int dimx_siz = nx;
  int dimy_siz = ny;
  int dimz_siz = nz;
  int start_i  = 0;
  int start_j  = 0;
  int start_k  = 0;
  
  // read in nc
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  if (is_parallel_netcdf == 1)
  {
    sprintf(ou_file, "%s/media.nc", output_dir);

    if (ierr=nc_create_par(ou_file, NC_CLOBBER | NC_NETCDF4, 
                      comm, MPI_INFO_NULL, &ncid)) M_NCRET(ierr);

    // reset dimx_siz etc
    dimx_siz = gdinfo->glob_nx;
    dimy_siz = gdinfo->glob_ny;
    dimz_siz = gdinfo->glob_nz;
    start_i  = gdinfo->nx1_to_glob_halo0;
    start_j  = gdinfo->ny1_to_glob_halo0;
    start_k  = gdinfo->nz1_to_glob_halo0;
  }
  else
  {
    sprintf(ou_file, "%s/media_%s.nc", output_dir, fname_coords);

    ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
    if (ierr != NC_NOERR){
      fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }
  }

  // define dimension
  ierr = nc_def_dim(ncid, "i", dimx_siz, &dimid[2]);
  ierr = nc_def_dim(ncid, "j", dimy_siz, &dimid[1]);
  ierr = nc_def_dim(ncid, "k", dimz_siz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, m3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
  }

  // attribute: index in output snapshot, index w ghost in thread
  if (is_parallel_netcdf == 0)
  {
    int l_start[] = { ni1, nj1, nk1 };
    nc_put_att_int(ncid,NC_GLOBAL,"local_index_of_first_physical_points",
                     NC_INT,CONST_NDIM,l_start);

    int g_start[] = { gni1, gnj1, gnk1 };
    nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                     NC_INT,CONST_NDIM,g_start);

    int l_count[] = { ni, nj, nk };
    nc_put_att_int(ncid,NC_GLOBAL,"count_of_physical_points",
                     NC_INT,CONST_NDIM,l_count);
  }

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++)
  {
    float *ptr = md->v4d + m3d_pos[ivar];
    size_t startp[] = { start_k, start_j, start_i };
    size_t countp[] = { nz, ny, nx };
    ierr = nc_put_vara_float(ncid, varid[ivar],startp,countp,ptr);
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
md_gen_test_el_vti(md_t *md)
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
	      md->c13[iptr] = 10.9620*1e9;//lam;
	      md->c33[iptr] = 18.0*1e9;//lam + 2.0f*mu;
	      md->c55[iptr] = 5.12*1e9;//mu;
        md->c66[iptr] = 7.1680*1e9;//mu;
        //-- Vp ~ sqrt(c11/rho) = 4098
        
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

int
md_ac_Vp_to_kappa(float *restrict rho, float *restrict kappa, size_t siz_volume)
{
  int ierr = 0;

  for (size_t iptr=0; iptr<siz_volume; iptr++) {
    if (rho[iptr] > 1e-10) {
      float Vp = kappa[iptr];
      kappa[iptr] = Vp * Vp * rho[iptr];
    } else {
      kappa[iptr] = 0.0;
    }
  }

  return ierr;
}

int 
md_vis_GMB_cal_Y(md_t *md, float freq, float fmin, float fmax)
{
  int ierr = 0;

  md->visco_GMB_freq = freq;
  md->visco_GMB_fmin = fmin;
  md->visco_GMB_fmax = fmax;

  int kmax = 2*md->nmaxwell-1;
  float wr = 2*PI*md->visco_GMB_freq;
  float wmin = 2*PI*md->visco_GMB_fmin;
  float wmax = 2*PI*md->visco_GMB_fmax;
  float wratio = wmax/wmin;

  int nmaxwell = md->nmaxwell;
  int nx = md->nx;
  int ny = md->ny;
  int nz = md->nz;
  int siz_line  = md->siz_iy;
  int siz_slice = md->siz_iz;
  

  float *wk = (float *) fdlib_mem_calloc_1d_float(kmax,0,
                                                  "medium_visco_iso_cal");
  float *YP = (float *) fdlib_mem_calloc_1d_float(nmaxwell,0,
                                                  "medium_visco_iso_cal");
  float *YS = (float *) fdlib_mem_calloc_1d_float(nmaxwell,0,
                                                  "medium_visco_iso_cal");
  float **GP = (float **) fdlib_mem_calloc_2l_float(kmax, nmaxwell, 0,
                                                   "medium_visco_iso_cal");
  float **GS = (float **) fdlib_mem_calloc_2l_float(kmax, nmaxwell, 0,
                                                   "medium_visco_iso_cal");

  for(int k=0; k<kmax; k++)
  {
      wk[k] = wmin*pow(wratio,k/(float)(kmax-1));
  }

  float *wl = md->wl;
  float *Qp = md->Qp;
  float *Qs = md->Qs; 
  float *lambda = md->lambda;
  float *mu = md->mu;
  float QP1,QS1,theta1,theta2,thetatmp,lam,muu;
  float R, kappa,muunrelax;


  for(int n=0; n<nmaxwell; n++)
  {
      wl[n] = wk[2*n];
  }

  for (size_t k=0; k<nz; k++)
  {
    for(size_t j=0 ;j<ny; j++)
    {
      for(size_t i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        QP1 = 1.0/round(Qp[iptr]);
        QS1 = 1.0/round(Qs[iptr]);
        lam = lambda[iptr];
        muu = mu[iptr];
        for(size_t m=0; m<kmax; m++)
          {
            for(size_t n=0; n<nmaxwell; n++) 
              {
                GP[m][n] = (wl[n]*wk[m]+pow(wl[n],2)*QP1)/(pow(wl[n],2)+pow(wk[m],2));
                GS[m][n] = (wl[n]*wk[m]+pow(wl[n],2)*QS1)/(pow(wl[n],2)+pow(wk[m],2));
              }
          }

        md_visco_LS(GP,YP,QP1,kmax,nmaxwell);
        md_visco_LS(GS,YS,QS1,kmax,nmaxwell);

        //P
        theta1 = 0.0;
        theta2 = 0.0;
        thetatmp = 0.0;

        for(int n=0; n<nmaxwell; n++){
            thetatmp = thetatmp+YP[n]/(1+pow(wr/wl[n],2));
        }
        theta1 = 1-thetatmp;

        thetatmp = 0.0;
        for(int n=0; n<nmaxwell; n++){
            thetatmp = thetatmp+YP[n]*wr/wl[n]/(1+pow(wr/wl[n],2));
        }
        theta2 = thetatmp;

        R = sqrt(pow(theta1,2)+pow(theta2,2));
        kappa = (lam+2*muu)*(R+theta1)/(2*pow(R,2));

        //S 
        theta1 = 0.0;
        theta2 = 0.0;
        thetatmp = 0.0;
    
        for(int n=0; n<nmaxwell; n++){
            thetatmp = thetatmp+YS[n]/(1+pow(wr/wl[n],2));
        }
        theta1 = 1-thetatmp;
    
        thetatmp = 0.0;
        for(int n=0; n<nmaxwell; n++){
            thetatmp = thetatmp+YS[n]*wr/wl[n]/(1+pow(wr/wl[n],2));
        }
        theta2 = thetatmp;
    
        R = sqrt(pow(theta1,2)+pow(theta2,2));
        muunrelax = muu*(R+theta1)/(2*pow(R,2));
    
        md->mu[iptr] = muunrelax;
    
        md->lambda[iptr] = kappa-2*muunrelax;
    
        for(int n=0; n<nmaxwell; n++){
          md->Ylam[n][iptr] = (1+2*md->mu[iptr]/md->lambda[iptr])*YP[n]
                              -2*md->mu[iptr]/md->lambda[iptr]*YS[n];
          md->Ymu[n][iptr]  = YS[n];
	  if(md->Ylam[n][iptr] > 1 || md->Ymu[n][iptr] > 1){
	    fprintf(stdout,"attention, the coef over the normal range!");
            fprintf(stdout,"Ylam[%d][%d][%d][%d]=%f,Ymu[%d][%d][%d][%d]=%f\n",
	      	      n,i,j,k,md->Ylam[n][iptr],n,i,j,k,md->Ymu[n][iptr]);
	  }
        }
      }
    }
  }
  
  //printf("lam=%f\n",lam);
  //printf("muu=%f\n",muu);
  //printf("kappa=%f\n",kappa);
  //printf("lambda=%f\n",md->lambda[150]);
  //printf("mu=%f\n",md->mu[150]);
  //printf("Ylam1=%f\n",md->Ylam[0][150]);
  //printf("Ylam2=%f\n",md->Ylam[1][150]);
  //printf("Ymu1=%f\n",md->Ymu[0][150]);
  //printf("Ymu2=%f\n",md->Ymu[1][150]);
  
  free(YP);
  free(YS);

  fdlib_mem_free_2l_float(GP, kmax, "visco_GP");
  fdlib_mem_free_2l_float(GS, kmax, "visco_GS");
  
  return ierr;

}

int 
md_visco_LS(float **restrict input, float *restrict output, float d, int m, int n)
{

  // G=md,m=(G^TG)^(-1)G^Td

  int ierr = 0;

  float trans[VISCO_LS_MAXSIZE][VISCO_LS_MAXSIZE];
  float multi[VISCO_LS_MAXSIZE][VISCO_LS_MAXSIZE];
  float inver[VISCO_LS_MAXSIZE][VISCO_LS_MAXSIZE];

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
        if(i==j)
          inver[i][j] = 1;
        else
          inver[i][j] = 0;
    }
  }
  
  // transposition
  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
        trans[i][j] = input[j][i];
    }
  }

  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      multi[i][j] = 0;
      for(int k=0;k<m;k++){
        multi[i][j] = multi[i][j]+trans[i][k]*input[k][j];
      }
    }
  }

 md_visco_LS_mat_inv(multi,inver,n);

  for(int i=0; i<n; i++){
    for(int j=0; j<m; j++){
      multi[i][j] = 0;
      for(int k=0; k<n; k++){
        multi[i][j] = multi[i][j]+inver[i][k]*trans[k][j];
      }
    }
  }

  float sum;
  for(int i=0; i<n; i++){
    sum = 0.0;
    for(int j=0; j<m; j++){
      sum = sum+multi[i][j]*d;
    }
    output[i] = sum;
  }

    return ierr;
}

int 
md_visco_LS_mat_inv(float matrix[][VISCO_LS_MAXSIZE], float inverse[][VISCO_LS_MAXSIZE], int n)
{
  int ierr = 0;
  
  float tmp;
  int j=0;

  for(int k=0; k<n; k++)
  {
    if(matrix[k][k] == 0)
    {
      for(int jj=k+1; jj<n; jj++)
      {
        j = jj;
        if(matrix[j][k]!=0)
            break;
      }
      if(j == n)
      {
        fprintf(stderr,"Error: Matrix is not inversible(Cal of visco coef)\n");
        fflush(stderr);
        exit(-1);
      }
      for(int i=0; i<n; i++)
      {
        tmp = matrix[k][i];
        matrix[k][i] = matrix[j][i];
        matrix[j][i] = tmp;
        tmp = inverse[k][i];
        inverse[k][i] = inverse[j][i];
        inverse[j][i] = tmp;
      }
    }
          
    tmp = matrix[k][k];
    for(int j=0; j<n; j++)
    {
        matrix[k][j] = matrix[k][j]/tmp;
        inverse[k][j] = inverse[k][j]/tmp;
    }
    
    for(int i=0; i<n; i++)
    {
      tmp = matrix[i][k];
      for(j=0;j<n;j++)
      {
        if(i==k)
          break;
        matrix[i][j] = matrix[i][j]-matrix[k][j]*tmp;
        inverse[i][j] = inverse[i][j]-inverse[k][j]*tmp;
      }
    }
  }

  return ierr;
}

int
md_gen_test_vis_iso(md_t *md)
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
  float *Qs = md->Qs;
  float *Qp = md->Qp;

  float Vp=3000.0;
  float Vs=2000.0;
  float rho=1500.0;
  float mu = Vs*Vs*rho;
  float lam = Vp*Vp*rho - 2.0*mu;

  for (size_t k=0; k<nz; k++)
  {
    for (size_t j=0; j<ny; j++)
    {
      for (size_t i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        lam3d[iptr] = lam;
         mu3d[iptr] = mu;
        rho3d[iptr] = rho;
        Qs[iptr] = 40;
        Qp[iptr] = 80;
      }
    }
  }

  return ierr;
}

// stress 2 strain for single trace, the input is also output

int
md_stress2strain_trace(float *tij, // time fastest, then cmp
                       int nt, md_t *md, size_t iptr)
{
  int ierr = 0;

  if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO)
  {
    float lam = md->lambda[iptr];
    float mu  = md->mu    [iptr];

    md_stress2strain_trace_el_iso(tij, nt, lam, mu);

  } else {

    fprintf(stdout, "--- warning: stress2strain for medium_type=%d has not been implemented\n",
           md->medium_type); fflush(stdout);
    // temp set to 0 with unsafe memset
    memset(tij, 0, nt * CONST_NDIM_2 * sizeof(float));
  }

  return ierr;
}

int
md_stress2strain_trace_el_iso(float *tij, // time fastest, then cmp
                              int nt, float lam, float mu)
{
  int ierr = 0;

  // use fake evt_x etc. since did not implement gather evt_x by mpi
  float evt_x = 0.0;
  float evt_y = 0.0;
  float evt_z = 0.0;
  float evt_d = 0.0;

  // cmp seq hard-coded, need to revise in the future
  float *Txx = tij + 0 * nt;
  float *Tyy = tij + 1 * nt;
  float *Tzz = tij + 2 * nt;
  float *Tyz = tij + 3 * nt;
  float *Txz = tij + 4 * nt;
  float *Txy = tij + 5 * nt;

  float E1, E2, E3;

  if (mu < 1e-16) {
    // need to check and revise
    E1 = 0.0;
    E2 = 1.0 / lam;
    E3 = 0.0;
  } else {
    E1 = (lam + mu) / (mu * ( 3.0 * lam + 2.0 * mu));
    E2 = - lam / ( 2.0 * mu * (3.0 * lam + 2.0 * mu));
    E3 = 1.0 / mu;
  }

  // conver to strain per time step
  for (int it = 0; it < nt; it++)
  {
    float E0 = E2 * (Txx[it] + Tyy[it] + Tzz[it]);

    Txx[it] = E0 - (E2 - E1) * Txx[it];
    Tyy[it] = E0 - (E2 - E1) * Tyy[it];
    Tzz[it] = E0 - (E2 - E1) * Tzz[it];
    Tyz[it] = 0.5 * E3 * Tyz[it];
    Txz[it] = 0.5 * E3 * Txz[it];
    Txy[it] = 0.5 * E3 * Txy[it];
  }

  return ierr;
}

int
md_stress2strain_trace_el_vti(float *tij, int nt,
                        float *restrict c11, float *restrict c13,
                        float *restrict c33, float *restrict c55,
                        float *restrict c66)
{
  //not implement

  return 0;
}

int
md_stress2strain_trace_el_aniso(float *tij, int nt,
                        float *restrict c11d, float *restrict c12d,
                        float *restrict c13d, float *restrict c14d,
                        float *restrict c15d, float *restrict c16d,
                        float *restrict c22d, float *restrict c23d,
                        float *restrict c24d, float *restrict c25d,
                        float *restrict c26d, float *restrict c33d,
                        float *restrict c34d, float *restrict c35d,
                        float *restrict c36d, float *restrict c44d,
                        float *restrict c45d, float *restrict c46d,
                        float *restrict c55d, float *restrict c56d,
                        float *restrict c66d)
{
  //not implement

  return 0;
}
