/*
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"
#include "netcdf_par.h"
#include "sacLib.h"

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "constants.h"
#include "fd_t.h"
#include "gd_t.h"
#include "io_funcs.h"

#ifndef M_NCRET
#define M_NCRET(ierr) { \
  fprintf(stderr,"io nc error: %s [%s:%d]\n", \
    nc_strerror(ierr),  __FILE__, __LINE__); exit(1);}
#endif

#ifndef M_NCERR
#define M_NCERR {fprintf(stderr,"io nc error\n"); exit(1);}
#endif

/*
 * export to a single file
 */

void
io_build_fname(char *out_dir,
               char *prefix,
               char *sufix,
               int *topoid,
               char *ou_fname)
{
  sprintf(ou_fname,"%s/%s_%d_%d%s",out_dir,prefix,topoid[0],topoid[1],sufix);

  //fprintf(stdout,"out_dir=%s\n",out_dir);
  //fprintf(stdout,"prefix=%s\n",prefix);
  //fprintf(stdout,"sufix=%s\n",sufix);
  //fprintf(stdout,"ou_fname=%s\n",ou_fname);
  //fflush(stdout);
}

void
io_build_fname_time(char *out_dir,
                    char *prefix,
                    char *sufix,
                    int *topoid,
                    int  it,
                    char *ou_fname)
{
  sprintf(ou_fname,"%s/%s_%d_%d_it%d%s",out_dir,prefix,topoid[0],topoid[1],it,sufix);
}

/*
 * append to a single file
 */

/*
void
io_snapshot_append(FILE *fp,
                   float *restrict var,
                   size_t ni,
                   size_t nj,
                   size_t nk,
                   int verbose)
{

}
*/

/*
 *  export all var3ds to netcdf file
 */
void
io_var3d_export_nc(char   *ou_file,
                   float  *restrict v3d,
                   size_t *restrict v3d_pos,
                   char  **restrict v3d_name,
                   int   number_of_vars,
                   char  **restrict coord_name,
                   int  nx,
                   int  ny,
                   int  nz)
{
  int ncid;
  int varid[number_of_vars];
  int dimid[CONST_NDIM];

  int ierr = nc_create(ou_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"ou_file=%s\n",ou_file);
    fprintf(stderr,"creat coord nc error: %s\n", nc_strerror(ierr));
    exit(-1);
  }

  // define dimension: last dim varis fastest for c nc file
  ierr = nc_def_dim(ncid, coord_name[0], nx, &dimid[2]);
  ierr = nc_def_dim(ncid, coord_name[1], ny, &dimid[1]);
  ierr = nc_def_dim(ncid, coord_name[2], nz, &dimid[0]);

  // define vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    ierr = nc_def_var(ncid, v3d_name[ivar], NC_FLOAT, CONST_NDIM, dimid, &varid[ivar]);
  }

  // end def
  ierr = nc_enddef(ncid);

  // add vars
  for (int ivar=0; ivar<number_of_vars; ivar++) {
    float *ptr = v3d + v3d_pos[ivar];
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
 * get next non-comment line
 */

int
io_get_nextline(FILE *fp, char *str, int length)
{
  int ierr = 0;

  do
  {
    if (fgets(str, length, fp) == NULL)
    {
      ierr = 1;
      return ierr;
    }
  } while (str[0] == '#' || str[0] == '\n');

  // remove newline char
  int len = strlen(str);

  if (len > 0 && str[len-1] == '\n') {
    str[len-1] = '\0';
  }

  // for debug:
  //fprintf(stdout," --return: %s\n", str);

  return ierr;
}

/*******************************************************************************
 * snapshot output
 *******************************************************************************/

void
io_snapshot_export_binary(char *fname,
                   float *restrict var,
                   int nx,
                   int ny,
                   int nz,
                   int *snap_indx,
                   int verbose)
{
  FILE *fp=fopen(fname,"wb");
  if (fp == NULL) {
    fprintf(stderr,"Error: can't create : %s\n", fname);
    exit(1);
  }

  // index triple
  int i1 = snap_indx[0];
  int j1 = snap_indx[1];
  int k1 = snap_indx[2];
  int ic = snap_indx[3];
  int jc = snap_indx[4];
  int kc = snap_indx[5];
  int di = snap_indx[6];
  int dj = snap_indx[7];
  int dk = snap_indx[8];

  for (int n3=0; n3<kc; n3++)
  {
    int k = k1 + n3 * dk;
    for (int n2=0; n2<jc; n2++)
    {
      int j = j1 + n2 * dj;
      for (int n1=0; n1<ic; n1++)
      {
        int i = i1 + n1 * di;
        int iptr = i + j * nx + k * nx * ny;
        fwrite(var+iptr, 1, sizeof(float), fp);
      }
    }
  }

  fclose(fp);
}

void
io_snapshot_locate(gd_t *gdinfo,
                   iosnap_t *iosnap,
                    int  number_of_snapshot,
                    char **snapshot_name,
                    int *snapshot_index_start,
                    int *snapshot_index_count,
                    int *snapshot_index_incre,
                    int *snapshot_time_start,
                    int *snapshot_time_incre,
                    int *snapshot_save_velocity,
                    int *snapshot_save_stress,
                    int *snapshot_save_strain,
                    int *snapshot_save_coord,
                    char *output_dir)
{
  // keep total number for parallel netcdf
  iosnap->num_snap_total  = number_of_snapshot;

  // malloc to max, num of snap will not be large
  if (number_of_snapshot > 0)
  {
    iosnap->fname_main = (char **) fdlib_mem_malloc_2l_char(number_of_snapshot,
                                    CONST_MAX_STRLEN,"snap_fname_main");
    iosnap->i1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->j1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->k1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->ni = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->nj = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->nk = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->di = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->dj = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->dk = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->it1 = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->dit = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_vel    = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_stress = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_strain = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->out_coord  = (int *) malloc(number_of_snapshot * sizeof(int));

    iosnap->i1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->j1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->k1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));

    iosnap->in_this_proc = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->ni_total = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->nj_total = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->nk_total = (int *) malloc(number_of_snapshot * sizeof(int));
  }

  // init

  iosnap->siz_max_wrk = 0;

  int isnap = 0;

  for (int n=0; n < number_of_snapshot; n++)
  {
    int iptr0 = n*CONST_NDIM;

    iosnap->ni_total[n] = snapshot_index_count[iptr0+0];
    iosnap->nj_total[n] = snapshot_index_count[iptr0+1];
    iosnap->nk_total[n] = snapshot_index_count[iptr0+2];

    // scan output k-index in this proc
    int gk1 = -1; int ngk =  0; int k_in_nc = 0;
    for (int n3=0; n3<snapshot_index_count[iptr0+2]; n3++)
    {
      int gk = snapshot_index_start[iptr0+2] + n3 * snapshot_index_incre[iptr0+2];
      if (gd_gindx_is_inner_k(gk,gdinfo) == 1)
      {
        if (gk1 == -1) {
          gk1 = gk;
          k_in_nc = n3;
        }
        ngk++;
      }
      if (gk > gdinfo->gnk2) break; // no need to larger k
    }

    // scan output j-index in this proc
    int gj1 = -1; int ngj =  0; int j_in_nc = 0;
    for (int n2=0; n2<snapshot_index_count[iptr0+1]; n2++)
    {
      int gj = snapshot_index_start[iptr0+1] + n2 * snapshot_index_incre[iptr0+1];
      if (gd_gindx_is_inner_j(gj,gdinfo) == 1)
      {
        if (gj1 == -1) {
          gj1 = gj;
          j_in_nc = n2;
        }
        ngj++;
      }
      if (gj > gdinfo->gnj2) break;
    }

    // scan output i-index in this proc
    int gi1 = -1; int ngi =  0; int i_in_nc = 0;
    for (int n1=0; n1<snapshot_index_count[iptr0+0]; n1++)
    {
      int gi = snapshot_index_start[iptr0+0] + n1 * snapshot_index_incre[iptr0+0];
      if (gd_gindx_is_inner_i(gi,gdinfo) == 1)
      {
        if (gi1 == -1) {
          gi1 = gi;
          i_in_nc = n1;
        }
        ngi++;
      }
      if (gi > gdinfo->gni2) break;
    }

    // if in this proc
    if (ngi>0 && ngj>0 && ngk>0)
    {
      iosnap->in_this_proc[n] = 1;

      iosnap->i1[n]  = gd_indx_glphy2lcext_i(gi1, gdinfo);
      iosnap->j1[n]  = gd_indx_glphy2lcext_j(gj1, gdinfo);
      iosnap->k1[n]  = gd_indx_glphy2lcext_k(gk1, gdinfo);
      iosnap->ni[n]  = ngi;
      iosnap->nj[n]  = ngj;
      iosnap->nk[n]  = ngk;
      iosnap->di[n]  = snapshot_index_incre[iptr0+0];
      iosnap->dj[n]  = snapshot_index_incre[iptr0+1];
      iosnap->dk[n]  = snapshot_index_incre[iptr0+2];

      iosnap->it1[n]  = snapshot_time_start[n];
      iosnap->dit[n]  = snapshot_time_incre[n];

      iosnap->out_vel   [n] = snapshot_save_velocity[n];
      iosnap->out_stress[n] = snapshot_save_stress[n];
      iosnap->out_strain[n] = snapshot_save_strain[n];
      iosnap->out_coord [n] = snapshot_save_coord [n];

      iosnap->i1_to_glob[n] = i_in_nc;
      iosnap->j1_to_glob[n] = j_in_nc;
      iosnap->k1_to_glob[n] = k_in_nc;

      sprintf(iosnap->fname_main[n],"%s/%s",output_dir,
                                                 snapshot_name[n]);
                                                 

      // for max wrk
      size_t snap_siz =  ngi * ngj * ngk;
      iosnap->siz_max_wrk = snap_siz > iosnap->siz_max_wrk ? 
                            snap_siz : iosnap->siz_max_wrk;

      isnap += 1;
    } // if in this
    else
    {
      iosnap->in_this_proc[n] = 0;

      iosnap->i1[n]  = -1;
      iosnap->j1[n]  = -1;
      iosnap->k1[n]  = -1;
      iosnap->ni[n]  = 0;
      iosnap->nj[n]  = 0;
      iosnap->nk[n]  = 0;
      iosnap->di[n]  = 1;
      iosnap->dj[n]  = 1;
      iosnap->dk[n]  = 1;

      iosnap->it1[n]  = snapshot_time_start[n];
      iosnap->dit[n]  = snapshot_time_incre[n];

      iosnap->out_vel   [n] = snapshot_save_velocity[n];
      iosnap->out_stress[n] = snapshot_save_stress[n];
      iosnap->out_strain[n] = snapshot_save_strain[n];
      iosnap->out_coord [n] = snapshot_save_coord [n];

      iosnap->i1_to_glob[n] = -1;
      iosnap->j1_to_glob[n] = -1;
      iosnap->k1_to_glob[n] = -1;

      sprintf(iosnap->fname_main[n],"%s/%s",output_dir,
                                                 snapshot_name[n]);
    }
  } // loop all snap

  iosnap->num_of_snap = isnap;
}

/*
 * combine init and creat to reduce funcs call
 */

int
io_snap_nc_create(iosnap_t    *iosnap,
                  iosnap_nc_t *iosnap_nc,
                  gd_t        *gd,
                  char *output_fname_part,
                  float *buff,
                  int is_parallel_netcdf,
                  MPI_Comm comm, 
                  int *topoid)
{
  int ierr = 0;

  size_t siz_line  = gd->siz_iy;
  size_t siz_slice = gd->siz_iz;

  float *x3d = gd->x3d;
  float *y3d = gd->y3d;
  float *z3d = gd->z3d;

  iosnap_nc->num_of_snap    = iosnap->num_of_snap;
  iosnap_nc->num_snap_total = iosnap->num_snap_total; 

  int num_snap_total = iosnap_nc->num_snap_total; 

  iosnap_nc->is_parallel = is_parallel_netcdf;

  iosnap_nc->in_this_proc = (int *)malloc(num_snap_total*sizeof(int));

  iosnap_nc->fname = (char **) fdlib_mem_malloc_2l_char(num_snap_total,
                                    CONST_MAX_STRLEN,"snap_fname");

  iosnap_nc->ncid = (int *)malloc(num_snap_total*sizeof(int));
  iosnap_nc->timeid = (int *)malloc(num_snap_total*sizeof(int));

  iosnap_nc->varid_V = (int *)malloc(num_snap_total*CONST_NDIM*sizeof(int));
  iosnap_nc->varid_T = (int *)malloc(num_snap_total*CONST_NDIM_2*sizeof(int));
  iosnap_nc->varid_E = (int *)malloc(num_snap_total*CONST_NDIM_2*sizeof(int));

  // will be used in put step
  iosnap_nc->cur_it = (int *)malloc(num_snap_total*sizeof(int));
  for (int n=0; n<num_snap_total; n++) {
    iosnap_nc->cur_it[n] = 0;
  }

  // for parallel netcdf
  iosnap_nc->start_i = (int *)malloc(num_snap_total*sizeof(int));
  iosnap_nc->start_j = (int *)malloc(num_snap_total*sizeof(int));
  iosnap_nc->start_k = (int *)malloc(num_snap_total*sizeof(int));

  int *in_this_proc = iosnap_nc->in_this_proc;

  int *ncid    = iosnap_nc->ncid;
  int *timeid  = iosnap_nc->timeid;
  int *varid_V = iosnap_nc->varid_V;
  int *varid_T = iosnap_nc->varid_T;
  int *varid_E = iosnap_nc->varid_E;
  int *start_i   = iosnap_nc->start_i;
  int *start_j   = iosnap_nc->start_j;
  int *start_k   = iosnap_nc->start_k;
  int varid_x, varid_y, varid_z;

  for (int n=0; n<num_snap_total; n++)
  {
    in_this_proc[n] = iosnap->in_this_proc[n];

    int dimid[4];
    int snap_i1  = iosnap->i1[n];
    int snap_j1  = iosnap->j1[n];
    int snap_k1  = iosnap->k1[n];
    int snap_ni  = iosnap->ni[n];
    int snap_nj  = iosnap->nj[n];
    int snap_nk  = iosnap->nk[n];
    int snap_di  = iosnap->di[n];
    int snap_dj  = iosnap->dj[n];
    int snap_dk  = iosnap->dk[n];

    int snap_out_V = iosnap->out_vel[n];
    int snap_out_T = iosnap->out_stress[n];
    int snap_out_E = iosnap->out_strain[n];
    int snap_out_coord = iosnap->out_coord[n];

    // default for out per proc
    int snap_dim_ni = snap_ni;
    int snap_dim_nj = snap_nj;
    int snap_dim_nk = snap_nk;
    start_i[n] = 0;
    start_j[n] = 0;
    start_k[n] = 0;

    // if parallel nc
    if (is_parallel_netcdf == 1)
    {
      sprintf(iosnap_nc->fname[n],"%s.nc", iosnap->fname_main[n]);
      if (nc_create_par(iosnap_nc->fname[n], NC_CLOBBER | NC_NETCDF4, 
                        comm, MPI_INFO_NULL, &ncid[n])) M_NCERR;
      // reset ni,nj,nk to total
      snap_dim_ni = iosnap->ni_total[n];
      snap_dim_nj = iosnap->nj_total[n];
      snap_dim_nk = iosnap->nk_total[n];
      start_i[n] = iosnap->i1_to_glob[n];
      start_j[n] = iosnap->j1_to_glob[n];
      start_k[n] = iosnap->k1_to_glob[n];
    } 
    // if output one file per proc
    else if (in_this_proc[n] == 1)
    {
      sprintf(iosnap_nc->fname[n],"%s_%s.nc", iosnap->fname_main[n], output_fname_part);
      if (nc_create(iosnap_nc->fname[n], NC_CLOBBER, &ncid[n])) M_NCERR;
    }
    // output per proc and not in this proc, continue to next
    else
    {
      continue;
    }

    if (nc_def_dim(ncid[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid[n], "k", snap_dim_nk    , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid[n], "j", snap_dim_nj    , &dimid[2])) M_NCERR;
    if (nc_def_dim(ncid[n], "i", snap_dim_ni    , &dimid[3])) M_NCERR;
    // time var
    if (nc_def_var(ncid[n], "time", NC_FLOAT, 1, dimid+0, &timeid[n])) M_NCERR;
    if (is_parallel_netcdf==1) {
      if (ierr=nc_var_par_access(ncid[n], timeid[n], NC_COLLECTIVE)) M_NCRET(ierr);
    }
    // if coord
    if (snap_out_coord==1) {
       if (nc_def_var(ncid[n],"x",NC_FLOAT,CONST_NDIM,dimid+1,&varid_x)) M_NCERR;
       if (nc_def_var(ncid[n],"y",NC_FLOAT,CONST_NDIM,dimid+1,&varid_y)) M_NCERR;
       if (nc_def_var(ncid[n],"z",NC_FLOAT,CONST_NDIM,dimid+1,&varid_z)) M_NCERR;
       if (is_parallel_netcdf==1) {
         if (ierr=nc_var_par_access(ncid[n], varid_x, NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n], varid_y, NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n], varid_z, NC_COLLECTIVE)) M_NCRET(ierr);
       }
    }
    // other vars
    if (snap_out_V==1) {
       if (nc_def_var(ncid[n],"Vx",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Vy",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Vz",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+2])) M_NCERR;
       if (is_parallel_netcdf==1) {
         if (ierr=nc_var_par_access(ncid[n],varid_V[n*CONST_NDIM+0], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_V[n*CONST_NDIM+1], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_V[n*CONST_NDIM+2], NC_COLLECTIVE)) M_NCRET(ierr);
       }
    }
    if (snap_out_T==1) {
       if (nc_def_var(ncid[n],"Txx",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Tyy",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Tzz",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+2])) M_NCERR;
       if (nc_def_var(ncid[n],"Txz",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+3])) M_NCERR;
       if (nc_def_var(ncid[n],"Tyz",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+4])) M_NCERR;
       if (nc_def_var(ncid[n],"Txy",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+5])) M_NCERR;
       if (is_parallel_netcdf==1) {
         if (ierr=nc_var_par_access(ncid[n],varid_T[n*CONST_NDIM_2+0], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_T[n*CONST_NDIM_2+1], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_T[n*CONST_NDIM_2+2], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_T[n*CONST_NDIM_2+3], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_T[n*CONST_NDIM_2+4], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_T[n*CONST_NDIM_2+5], NC_COLLECTIVE)) M_NCRET(ierr);
       }
    }
    if (snap_out_E==1) {
       if (nc_def_var(ncid[n],"Exx",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Eyy",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Ezz",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+2])) M_NCERR;
       if (nc_def_var(ncid[n],"Exz",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+3])) M_NCERR;
       if (nc_def_var(ncid[n],"Eyz",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+4])) M_NCERR;
       if (nc_def_var(ncid[n],"Exy",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+5])) M_NCERR;
       if (is_parallel_netcdf==1) {
         if (ierr=nc_var_par_access(ncid[n],varid_E[n*CONST_NDIM_2+0], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_E[n*CONST_NDIM_2+1], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_E[n*CONST_NDIM_2+2], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_E[n*CONST_NDIM_2+3], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_E[n*CONST_NDIM_2+4], NC_COLLECTIVE)) M_NCRET(ierr);
         if (ierr=nc_var_par_access(ncid[n],varid_E[n*CONST_NDIM_2+5], NC_COLLECTIVE)) M_NCRET(ierr);
       }
    }

    // attribute: index in output snapshot, index w ghost in thread
    if (is_parallel_netcdf == 0)
    {
      int g_start[] = { iosnap->i1_to_glob[n],
                        iosnap->j1_to_glob[n],
                        iosnap->k1_to_glob[n] };
      nc_put_att_int(ncid[n],NC_GLOBAL,"first_index_to_snapshot_output",
                     NC_INT,CONST_NDIM,g_start);

      int l_start[] = { snap_i1, snap_j1, snap_k1 };
      nc_put_att_int(ncid[n],NC_GLOBAL,"first_index_in_this_thread_with_ghosts",
                     NC_INT,CONST_NDIM,l_start);

      int l_count[] = { snap_di, snap_dj, snap_dk };
      nc_put_att_int(ncid[n],NC_GLOBAL,"index_stride_in_this_thread",
                     NC_INT,CONST_NDIM,l_count);
      nc_put_att_int(ncid[n],NC_GLOBAL,"coords_of_mpi_topo",
                     NC_INT,2,topoid);
    }
    else
    {
      // set collective I/O globally (for all variables)
      //if (ierr=nc_var_par_access(ncid[n], NC_GLOBAL, NC_COLLECTIVE)) M_NCRET(ierr);
    }

    if (nc_enddef(ncid[n])) M_NCERR;

    // write out coord if out coord
    if (snap_out_coord==1)
    {
      size_t startp[] = { start_k[n], start_j[n], start_i[n] };
      size_t countp[] = { iosnap->nk[n], iosnap->nj[n], iosnap->ni[n] };

      io_snap_pack_buff(x3d,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
      //fprintf(stdout,"-- 31\n"); fflush(stdout);
      if (nc_put_vara_float(iosnap_nc->ncid[n],varid_x,startp,countp,buff)) M_NCERR;
      //fprintf(stdout,"-- 32\n"); fflush(stdout);

      io_snap_pack_buff(y3d,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
      nc_put_vara_float(iosnap_nc->ncid[n],varid_y,startp,countp,buff);

      io_snap_pack_buff(z3d,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
      nc_put_vara_float(iosnap_nc->ncid[n],varid_z,startp,countp,buff);
    }
  } // loop snap

  return ierr;
}

/*
 * 
 */

int
io_snap_nc_put(iosnap_t *iosnap,
               iosnap_nc_t *iosnap_nc,
               gd_t    *gdinfo,
               md_t    *md,
               wav_t   *wav,
               float *restrict w4d,
               float *restrict buff,
               int   nt_total,
               int   it,
               float time,
               int is_run_out_vel,     // for stg, out vel and stress at sep call
               int is_run_out_stress,  // 
               int is_incr_cur_it)     // for stg, should output cur_it once
{
  int ierr = 0;

  //int num_of_snap = iosnap->num_of_snap;
  int num_snap_total = iosnap->num_snap_total;
  size_t siz_line  = gdinfo->siz_iy;
  size_t siz_slice = gdinfo->siz_iz;

  int *start_i   = iosnap_nc->start_i;
  int *start_j   = iosnap_nc->start_j;
  int *start_k   = iosnap_nc->start_k;

  for (int n=0; n<num_snap_total; n++)
  {
    // skip if not in this proc
    if (iosnap->in_this_proc[n] == 0) {
      //fprintf(stdout,"-- 11,n=%d, in_this_proc=%d\n",n,iosnap->in_this_proc[n]); fflush(stdout);
      continue;
    }

    //fprintf(stdout,"-- 12,n=%d, in_this_proc=%d\n",n,iosnap->in_this_proc[n]); fflush(stdout);

    int snap_i1  = iosnap->i1[n];
    int snap_j1  = iosnap->j1[n];
    int snap_k1  = iosnap->k1[n];
    int snap_ni  = iosnap->ni[n];
    int snap_nj  = iosnap->nj[n];
    int snap_nk  = iosnap->nk[n];
    int snap_di  = iosnap->di[n];
    int snap_dj  = iosnap->dj[n];
    int snap_dk  = iosnap->dk[n];

    int snap_it1 = iosnap->it1[n];
    int snap_dit = iosnap->dit[n];

    int snap_out_V = iosnap->out_vel[n];
    int snap_out_T = iosnap->out_stress[n];
    int snap_out_E = iosnap->out_strain[n];

    int snap_it_mod = (it - snap_it1) % snap_dit;
    int snap_it_num = (it - snap_it1) / snap_dit;
    int snap_nt_total = (nt_total - snap_it1) / snap_dit;

    int snap_max_num = snap_ni * snap_nj * snap_nk;

    if (it>=snap_it1 && snap_it_num<=snap_nt_total && snap_it_mod==0)
    {
      //size_t startp[] = { iosnap_nc->cur_it[n], 0, 0, 0 };
      size_t startp[] = { iosnap_nc->cur_it[n], start_k[n], start_j[n], start_i[n] };
      size_t countp[] = { 1, snap_nk, snap_nj, snap_ni };
      size_t start_tdim = iosnap_nc->cur_it[n];

      //fprintf(stdout,"-- 20,n=%d, startp=%d,%d,%d,%d count=%d,%d,%d,%d\n",n,
      //        startp[0],startp[1],startp[2],startp[3],
      //        countp[0],countp[1],countp[2],countp[3]); fflush(stdout);

      // put time var
      //fprintf(stdout,"-- 20b,n=%d, ncid=%d,tid=%d,vxid=%d,vyid=%d\n",
      //          n,iosnap_nc->ncid[n],iosnap_nc->timeid[n],
      //          iosnap_nc->varid_V[n*CONST_NDIM+0],
      //          iosnap_nc->varid_V[n*CONST_NDIM+1]
      //          );
      //fflush(stdout);
      if (ierr=nc_put_var1_float(iosnap_nc->ncid[n],iosnap_nc->timeid[n],&start_tdim,&time)) M_NCRET(ierr);
      //fprintf(stdout,"-- 21,n=%d, start_tdim=%d\n",start_tdim); fflush(stdout);

      // vel
      if (is_run_out_vel == 1 && snap_out_V==1)
      {
        io_snap_pack_buff(w4d + wav->Vx_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        if (nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+0],
              startp,countp,buff)) M_NCERR;

        io_snap_pack_buff(w4d + wav->Vy_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+1],
              startp,countp,buff);

        io_snap_pack_buff(w4d + wav->Vz_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+2],
              startp,countp,buff);
      }
      if (is_run_out_stress==1 && snap_out_T==1)
      {
        io_snap_pack_buff(w4d + wav->Txx_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*CONST_NDIM_2+0],
              startp,countp,buff);

        io_snap_pack_buff(w4d + wav->Tyy_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*CONST_NDIM_2+1],
              startp,countp,buff);
        
        io_snap_pack_buff(w4d + wav->Tzz_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*CONST_NDIM_2+2],
              startp,countp,buff);
        
        io_snap_pack_buff(w4d + wav->Txz_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*CONST_NDIM_2+3],
              startp,countp,buff);

        io_snap_pack_buff(w4d + wav->Tyz_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*CONST_NDIM_2+4],
              startp,countp,buff);

        io_snap_pack_buff(w4d + wav->Txy_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n*CONST_NDIM_2+5],
              startp,countp,buff);
      }
      if (is_run_out_stress==1 && snap_out_E==1)
      {
        // convert to strain
        md_stress2strain_snap_pack(md,
                                   w4d + wav->Txx_pos,
                                   w4d + wav->Tyy_pos,
                                   w4d + wav->Tzz_pos,
                                   w4d + wav->Tyz_pos,
                                   w4d + wav->Txz_pos,
                                   w4d + wav->Txy_pos,
                                   buff + wav->Txx_pos,
                                   buff + wav->Tyy_pos,
                                   buff + wav->Tzz_pos,
                                   buff + wav->Tyz_pos,
                                   buff + wav->Txz_pos,
                                   buff + wav->Txy_pos,
                                   siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
                                   snap_dj,snap_k1,snap_nk,snap_dk);
        // export
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*CONST_NDIM_2+0],
              startp,countp,buff+wav->Txx_pos);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*CONST_NDIM_2+1],
              startp,countp,buff+wav->Tyy_pos);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*CONST_NDIM_2+2],
              startp,countp,buff+wav->Tzz_pos);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*CONST_NDIM_2+3],
              startp,countp,buff+wav->Txz_pos);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*CONST_NDIM_2+4],
              startp,countp,buff+wav->Tyz_pos);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_E[n*CONST_NDIM_2+5],
              startp,countp,buff+wav->Txy_pos);
      }

      if (is_incr_cur_it == 1) {
        iosnap_nc->cur_it[n] += 1;
      }

    } // if it
  } // loop snap

  return ierr;
}


int
io_snap_stress_to_strain_eliso(float *restrict lam3d,
                               float *restrict mu3d,
                               float *restrict Txx,
                               float *restrict Tyy,
                               float *restrict Tzz,
                               float *restrict Tyz,
                               float *restrict Txz,
                               float *restrict Txy,
                               float *restrict Exx,
                               float *restrict Eyy,
                               float *restrict Ezz,
                               float *restrict Eyz,
                               float *restrict Exz,
                               float *restrict Exy,
                               size_t siz_line,
                               size_t siz_slice,
                               int starti,
                               int counti,
                               int increi,
                               int startj,
                               int countj,
                               int increj,
                               int startk,
                               int countk,
                               int increk)
{
  size_t iptr_snap=0;
  size_t i,j,k,iptr,iptr_j,iptr_k;
  float lam,mu,E1,E2,E3,E0;

  for (int n3=0; n3<countk; n3++)
  {
    k = startk + n3 * increk;
    iptr_k = k * siz_slice;
    for (int n2=0; n2<countj; n2++)
    {
      j = startj + n2 * increj;
      iptr_j = j * siz_line + iptr_k;

      for (int n1=0; n1<counti; n1++)
      {
        i = starti + n1 * increi;
        iptr = i + iptr_j;

        lam = lam3d[iptr];
        mu  =  mu3d[iptr];

        E1 = (lam + mu) / (mu * ( 3.0 * lam + 2.0 * mu));
        E2 = - lam / ( 2.0 * mu * (3.0 * lam + 2.0 * mu));
        E3 = 1.0 / mu;

        E0 = E2 * (Txx[iptr] + Tyy[iptr] + Tzz[iptr]);

        Exx[iptr_snap] = E0 - (E2 - E1) * Txx[iptr];
        Eyy[iptr_snap] = E0 - (E2 - E1) * Tyy[iptr];
        Ezz[iptr_snap] = E0 - (E2 - E1) * Tzz[iptr];
        Eyz[iptr_snap] = 0.5 * E3 * Tyz[iptr];
        Exz[iptr_snap] = 0.5 * E3 * Txz[iptr];
        Exy[iptr_snap] = 0.5 * E3 * Txy[iptr];

        iptr_snap++;
      } // i
    } //j
  } //k

  return 0;
}

int
io_snap_pack_buff(float *restrict var,
                  size_t siz_line,
                  size_t siz_slice,
                  int starti,
                  int counti,
                  int increi,
                  int startj,
                  int countj,
                  int increj,
                  int startk,
                  int countk,
                  int increk,
                  float *restrict buff)
{
  size_t iptr_snap=0;
  for (int n3=0; n3<countk; n3++)
  {
    size_t k = startk + n3 * increk;
    for (int n2=0; n2<countj; n2++)
    {
      size_t j = startj + n2 * increj;
      for (int n1=0; n1<counti; n1++)
      {
        size_t i = starti + n1 * increi;
        size_t iptr = i + j * siz_line + k * siz_slice;
        buff[iptr_snap] = var[iptr];
        iptr_snap++;
      }
    }
  }

  return 0;
}

int
io_snap_nc_close(iosnap_nc_t *iosnap_nc)
{
  int ierr = 0;

  for (int n=0; n < iosnap_nc->num_snap_total; n++)
  {
    if (iosnap_nc->is_parallel==1 || iosnap_nc->in_this_proc[n]==1) {
      if (ierr=nc_close(iosnap_nc->ncid[n])) M_NCRET(ierr);
    }
  }
  return(ierr);
}

int
iosnap_print(iosnap_t *iosnap)
{    
  fprintf(stdout, "--> num_of_snap = %d\n", iosnap->num_of_snap);
  fprintf(stdout, "--> num_snap_total = %d\n", iosnap->num_snap_total);
  fprintf(stdout, "#  ishere i0 j0 k0 ni nj nk di dj dk it0 dit vel stress strain gi1 gj1 gk1\n");
  for (int n=0; n < iosnap->num_snap_total; n++)
  {
    fprintf(stdout, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
              n,
              iosnap->in_this_proc[n],
              iosnap->i1[n], iosnap->j1[n], iosnap->k1[n],
              iosnap->ni[n], iosnap->nj[n], iosnap->nk[n],
              iosnap->di[n], iosnap->dj[n], iosnap->dk[n],
              iosnap->it1[n], iosnap->dit[n], 
              iosnap->out_vel[n],
              iosnap->out_stress[n],
              iosnap->out_strain[n],
              iosnap->i1_to_glob[n],
              iosnap->j1_to_glob[n],
              iosnap->k1_to_glob[n]);
  }

  return(0);
}

/*******************************************************************************
 * station output
 *******************************************************************************/

/*
 * read in station list file, prototype, not finished
 */

//int
//io_read_station_list(char *in_filenm, int *sta_num, char ***p_sta_name, float **p_sta_xyz)
//{
//  FILE *fp;
//  char line[500];
//
//  if (!(fp = fopen (in_filenm, "rt")))
//	{
//	    fprintf (stdout, "Cannot open input file %s\n", in_filenm);
//	    fflush (stdout);
//	    return(1);
//	}
//
//  // first round: get valid num
//  int sta_num = 0;
//  // scan line
//  while ( fgets(line,500,fp) )
//  {
//    // skip comment
//    if (line[0]=='#') continue;
//    sta_num += 1;
//  }
//
//  // alloc
//  char **sta_name;
//  float *sta_xyz;
//
//  sta_name = (char **) fdlib_mem_malloc_2l_char(sta_num,CONST_MAX_STRLEN,"io_read_station_list");
//  sta_xyz  = (float *) fdlib_mem_malloc_1d(sta_num*CONST_NDIM*sizeof(float),"io_read_station_list");
//
//  // second round: read values
//
//  fseek(fp, 0, SEEK_SET);
//
//  int ir=0;
//  while ( fgets(line,500,fp) )
//  {
//    // skip comment
//    if (line[0]=='#') continue;
//
//    sscanf(line, "%s %f %f %f", sta_name[ir],
//                      sta_xzy+ir*CONST_NDIM,
//                      sta_xzy+ir*CONST_NDIM+1,
//                      sta_xzy+ir*CONST_NDIM+2);
//    ir += 1;
//  }
//
//  fclose(fp);
//
//  *p_sta_num  = sta_num;
//  *p_sta_xyz  = sta_xyz;
//  *p_sta_name = sta_name;
//
//  return(0);
//}

int
io_recv_read_locate(gd_t     *gd,
                    iorecv_t *iorecv,
                    int       nt_total,
                    int       nt_per_out,
                    int       file_type_sac,
                    int       save_velocity,
                    int       save_stress,
                    int       save_strain,
                    char     *in_filenm,
                    MPI_Comm  comm,
                    int       myid,
                    int       verbose)
{

  //int file_type_sac = 0; // for test purpose, hard coded to output nc now

  // keep par
  iorecv->file_type_sac = file_type_sac;
  iorecv->save_velocity = save_velocity;
  iorecv->save_stress   = save_stress;
  iorecv->save_strain   = save_strain;
  iorecv->max_nt        = nt_total;
  iorecv->nt_per_out    = nt_per_out;

  // init values
  iorecv->it_to_this   = -1; // for += in _keep func
  iorecv->it0_to_start =  0; // for nc put

  // reset nt_per_out if nt_total is very small
  if (nt_per_out < 0 || nt_per_out > nt_total) {
    iorecv->nt_per_out = nt_total;
  }

  // reset nt_per_out for sac as sac is written only once
  if (file_type_sac == 1) {
    iorecv->nt_per_out = nt_total;
  }

  // default value
  iorecv->nt_this_out = iorecv->nt_per_out;

  // open input file and read
  FILE *fp;
  if (!(fp = fopen (in_filenm, "rt")))
	{
	    fprintf (stdout, "Cannot open input file %s\n", in_filenm);
	    fflush (stdout);
	    return(1);
	}

  // number of station
  int nr;
  char line[500];

  io_get_nextline(fp, line, 500);
  sscanf(line, "%d", &nr);

  //fprintf(stdout, "-- nr=%d\n", nr);
  iorecv->total_number = nr;

  // alloc all, check fail in the future
  iorecv_one_t *recvone      = (iorecv_one_t *)malloc(nr * sizeof(iorecv_one_t));

  // read coord and locate

  int ir=0;
  int nr_this = 0; // in this thread
  int flag_coord;
  int flag_depth;

  for (ir=0; ir<nr; ir++)
  {
    float rx, ry, rz;
    float rx_inc, ry_inc, rz_inc; // shift in computational space grid
    int ix, iy, iz; // global index

    // read one line
    io_get_nextline(fp, line, 500);

    // get values
    sscanf(line, "%s %d %d %g %g %g", 
              recvone[ir].name, &flag_coord, &flag_depth, &rx, &ry, &rz);
    
    // by grid index
    if (flag_coord == 0)
    {
      // if sz is relative to surface, convert to normal index
      if (flag_depth == 1) {
        rz = gd->gnk2 - rz;
      }

      // do not take nearest value, but use smaller value
      ix = (int) (rx + 0.0);
      iy = (int) (ry + 0.0);
      iz = (int) (rz + 0.0);
      rx_inc = rx - ix;
      ry_inc = ry - iy;
      rz_inc = rz - iz;
    }
    // by axis
    else
    {
      // convert coord to global index
      if (gd->type == GD_TYPE_CURV)
      {
        // if rz is depth, convert to axis when it is in this thread
        if (flag_depth == 1) {
          gd_curv_depth_to_axis(gd,rx,ry,&rz,comm,myid);
        }
        gd_curv_coord_to_glob_indx(gd,rx,ry,rz,comm,myid,
                               &ix,&iy,&iz,&rx_inc,&ry_inc,&rz_inc);
      }
      else if (gd->type == GD_TYPE_CART)
      {
        // if sz is depth, convert to axis
        if (flag_depth == 1) {
          rz = gd->z1d[gd->nk2] - rz;
        }
        gd_cart_coord_to_glob_indx(gd,rx,ry,rz,comm,myid,
                               &ix,&iy,&iz,&rx_inc,&ry_inc,&rz_inc);
      }

      // conver minus shift to plus to use linear interp with all grid in this thread
      //    there may be problem if the receiver is located just bewteen two mpi block
      //    we should exchange data first then before save receiver waveform
      if (rx_inc < 0.0) {
        rx_inc = 1.0 + rx_inc;
        ix -= 1;
      }
      if (ry_inc < 0.0) {
        ry_inc = 1.0 + ry_inc;
        iy -= 1;
      }
      if (rz_inc < 0.0) {
        rz_inc = 1.0 + rz_inc;
        iz -= 1;
      }
    }

    if (gd_gindx_is_inner(ix,iy,iz,gd) == 1)
    {
      // convert to local index w ghost
      int i_local = gd_indx_glphy2lcext_i(ix,gd);
      int j_local = gd_indx_glphy2lcext_j(iy,gd);
      int k_local = gd_indx_glphy2lcext_k(iz,gd);

      // get coord
      if (flag_coord == 0)
      {
        rx = gd_coord_get_x(gd,i_local,j_local,k_local);
        ry = gd_coord_get_y(gd,i_local,j_local,k_local);
        rz = gd_coord_get_z(gd,i_local,j_local,k_local);
      }

      //-- same for coord?
      //int ptr_this = nr_this * CONST_NDIM;

      iorecv_one_t *this_recv = recvone + nr_this;

      sprintf(this_recv->name, "%s", recvone[ir].name);

      // set seq_id
      this_recv->seq_id = ir;
      // get coord
      this_recv->x = rx;
      this_recv->y = ry;
      this_recv->z = rz;
      // set point and shift
      this_recv->i=i_local;
      this_recv->j=j_local;
      this_recv->k=k_local;
      this_recv->di = rx_inc;
      this_recv->dj = ry_inc;
      this_recv->dk = rz_inc;

      this_recv->indx1d[0] = i_local   + j_local     * gd->siz_iy + k_local * gd->siz_iz;
      this_recv->indx1d[1] = i_local+1 + j_local     * gd->siz_iy + k_local * gd->siz_iz;
      this_recv->indx1d[2] = i_local   + (j_local+1) * gd->siz_iy + k_local * gd->siz_iz;
      this_recv->indx1d[3] = i_local+1 + (j_local+1) * gd->siz_iy + k_local * gd->siz_iz;
      this_recv->indx1d[4] = i_local   + j_local     * gd->siz_iy + (k_local+1) * gd->siz_iz;
      this_recv->indx1d[5] = i_local+1 + j_local     * gd->siz_iy + (k_local+1) * gd->siz_iz;
      this_recv->indx1d[6] = i_local   + (j_local+1) * gd->siz_iy + (k_local+1) * gd->siz_iz;
      this_recv->indx1d[7] = i_local+1 + (j_local+1) * gd->siz_iy + (k_local+1) * gd->siz_iz;

      //fprintf(stdout,"== ir_this=%d,name=%s,i=%d,j=%d,k=%d\n",
      //      nr_this,sta_name[nr_this],i_local,j_local,k_local); fflush(stdout);

      nr_this += 1;
    }
  }

  fclose(fp);

  iorecv->nr_here = nr_this;
  iorecv->recvone = recvone;

  // malloc seismo
  for (int ir=0; ir < iorecv->nr_here; ir++)
  {
    recvone = iorecv->recvone + ir;
    if (save_velocity==1) {
      recvone->vi = (float *) malloc(CONST_NDIM * iorecv->nt_per_out * sizeof(float));
    }
    if (save_stress==1 || save_strain==1) {
      recvone->tij = (float *) malloc(CONST_NDIM_2 * iorecv->nt_per_out * sizeof(float));
    }
  }

  return(0);
}

int
io_recv_nc_create(iorecv_t *iorecv,
                  float stept,
                  int is_parallel_netcdf,
                  MPI_Comm comm, 
                  int myid,
                  char *fname_mpi,
                  char *output_dir)
{
  int ierr = 0;

  // file name
  char ou_file[CONST_MAX_STRLEN];

  // dim size and index
  int ncid;
  int dimnr_siz;
  int start_ir;

  if (is_parallel_netcdf == 1)
  {
    sprintf(ou_file, "%s/station.nc", output_dir);

    if (ierr=nc_create_par(ou_file, NC_CLOBBER | NC_NETCDF4, 
                      comm, MPI_INFO_NULL, &ncid)) M_NCRET(ierr);

    // set dim and index
    dimnr_siz = iorecv->total_number;
  }
  else
  {
    sprintf(ou_file, "%s/station_%s.nc", output_dir, fname_mpi);

    ierr = nc_create(ou_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
    if (ierr != NC_NOERR){
      fprintf(stderr,"station creat nc error: %s\n", nc_strerror(ierr));
      exit(-1);
    }

    dimnr_siz = iorecv->nr_here;
  }

  // define dimension
  int dimid[2];
  int dimid_nmlen;
  ierr = nc_def_dim(ncid, "time", iorecv->max_nt, &dimid[1]); // time varies fastest
  ierr = nc_def_dim(ncid, "nr"  , dimnr_siz,      &dimid[0]);
  ierr = nc_def_dim(ncid, "name_length", CONST_MAX_NMLEN, &dimid_nmlen);

  // dim for station name as char array
  int dimid_stanm[] = { dimid[0], dimid_nmlen };

  // define const vars
  int varid_t;
  int varid_seq;
  int varid_stanm;
  if (ierr=nc_def_var(ncid,"time"         ,NC_FLOAT,1,dimid+1,&varid_t))   M_NCRET(ierr);
  if (ierr=nc_def_var(ncid,"seq_no"       ,NC_INT  ,1,dimid+0,&varid_seq)) M_NCRET(ierr);
  if (ierr=nc_def_var(ncid,"station_name" ,NC_CHAR ,2,dimid_stanm,&varid_stanm)) M_NCRET(ierr);

  // define vars
  if (iorecv->save_velocity == 1)
  {
    if (nc_def_var(ncid,"Vx",NC_FLOAT,2,dimid,iorecv->varid_vi+0)) M_NCERR;
    if (nc_def_var(ncid,"Vy",NC_FLOAT,2,dimid,iorecv->varid_vi+1)) M_NCERR;
    if (nc_def_var(ncid,"Vz",NC_FLOAT,2,dimid,iorecv->varid_vi+2)) M_NCERR;
  }

  if (iorecv->save_stress == 1)
  {
    if (nc_def_var(ncid,"Txx",NC_FLOAT,2,dimid,iorecv->varid_tij+0)) M_NCERR;
    if (nc_def_var(ncid,"Tyy",NC_FLOAT,2,dimid,iorecv->varid_tij+1)) M_NCERR;
    if (nc_def_var(ncid,"Tzz",NC_FLOAT,2,dimid,iorecv->varid_tij+2)) M_NCERR;
    if (nc_def_var(ncid,"Tyz",NC_FLOAT,2,dimid,iorecv->varid_tij+3)) M_NCERR;
    if (nc_def_var(ncid,"Txz",NC_FLOAT,2,dimid,iorecv->varid_tij+4)) M_NCERR;
    if (nc_def_var(ncid,"Txy",NC_FLOAT,2,dimid,iorecv->varid_tij+5)) M_NCERR;
  }

  if (iorecv->save_strain == 1)
  {
    if (nc_def_var(ncid,"Exx",NC_FLOAT,2,dimid,iorecv->varid_eij+0)) M_NCERR;
    if (nc_def_var(ncid,"Eyy",NC_FLOAT,2,dimid,iorecv->varid_eij+1)) M_NCERR;
    if (nc_def_var(ncid,"Ezz",NC_FLOAT,2,dimid,iorecv->varid_eij+2)) M_NCERR;
    if (nc_def_var(ncid,"Eyz",NC_FLOAT,2,dimid,iorecv->varid_eij+3)) M_NCERR;
    if (nc_def_var(ncid,"Exz",NC_FLOAT,2,dimid,iorecv->varid_eij+4)) M_NCERR;
    if (nc_def_var(ncid,"Exy",NC_FLOAT,2,dimid,iorecv->varid_eij+5)) M_NCERR;
  }

  if (nc_enddef(ncid)) M_NCERR;

  // write time var
  if (is_parallel_netcdf == 0 || myid==0)
  {
    size_t start_tdim;
    float  time;
    for (int it=0; it<iorecv->max_nt; it++) {
      time       = it * stept;
      start_tdim = it;
      if (ierr=nc_put_var1_float(ncid,varid_t,&start_tdim,&time)) M_NCRET(ierr);
    }
  }

  // write other const vars
  size_t startp[] = { 0, 0               };
  size_t countp[] = { 1, CONST_MAX_NMLEN };
  size_t start1d;
  char stanm[CONST_MAX_NMLEN];

  for (int ir=0; ir < iorecv->nr_here; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;

    if (is_parallel_netcdf == 1) {
      start1d = this_recv->seq_id;
      startp[0] = this_recv->seq_id;
    } else {
      start1d   = ir;
      startp[0] = ir;
    }

    // seq no
    if (ierr=nc_put_var1_int(ncid,varid_seq,&start1d,&(this_recv->seq_id))) M_NCRET(ierr);

    // station name
    sprintf(stanm,"%s",this_recv->name);
    for (int j=strlen(stanm); j<CONST_MAX_NMLEN; j++) {
      stanm[j] = '\0';
    }
    if (ierr=nc_put_vara_text(ncid,varid_stanm,startp,countp,stanm)) M_NCRET(ierr);
  }

  // save to iorecv
  iorecv->ncid = ncid;

  return(ierr);
}

int
io_recv_keep(iorecv_t *iorecv, float *restrict w4d,
             int it, int siz_icmp)
{
  int ierr = 0;

  float Lx1, Lx2, Ly1, Ly2, Lz1, Lz2;

  // increase it_to_this
  iorecv->it_to_this += 1;

  // reset it_to_this if overflow
  if (iorecv->it_to_this >= iorecv->nt_this_out)
  {
    iorecv->it_to_this   = 0;
    iorecv->it0_to_start = it;

    // reset nt_this_out if no eough it left
    if (iorecv->nt_per_out > iorecv->max_nt - it) {
      iorecv->nt_this_out = iorecv->max_nt - it;
    }
  }
  // if it_to_this == nt_per_out - 1, should output after this func

  for (int n=0; n < iorecv->nr_here; n++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + n;
    int *indx1d = this_recv->indx1d;

    // get coef of linear interp
    Lx2 = this_recv->di; Lx1 = 1.0 - Lx2;
    Ly2 = this_recv->dj; Ly1 = 1.0 - Ly2;
    Lz2 = this_recv->dk; Lz1 = 1.0 - Lz2;

    if (iorecv->save_velocity==1)
    {
      for (int icmp=0; icmp < CONST_NDIM; icmp++)
      {
        int iptr_sta = icmp * iorecv->nt_per_out + iorecv->it_to_this;
        size_t iptr_cmp = icmp * siz_icmp;
        this_recv->vi[iptr_sta] = 
            w4d[iptr_cmp + indx1d[0]] * Lx1 * Ly1 * Lz1
          + w4d[iptr_cmp + indx1d[1]] * Lx2 * Ly1 * Lz1
          + w4d[iptr_cmp + indx1d[2]] * Lx1 * Ly2 * Lz1
          + w4d[iptr_cmp + indx1d[3]] * Lx2 * Ly2 * Lz1
          + w4d[iptr_cmp + indx1d[5]] * Lx1 * Ly1 * Lz2
          + w4d[iptr_cmp + indx1d[5]] * Lx2 * Ly1 * Lz2
          + w4d[iptr_cmp + indx1d[6]] * Lx1 * Ly2 * Lz2
          + w4d[iptr_cmp + indx1d[7]] * Lx2 * Ly2 * Lz2;
      }
    }

    if (iorecv->save_stress==1 || iorecv->save_strain==1)
    {
      // save stress
      for (int icmp=0; icmp < CONST_NDIM_2; icmp++)
      {
        int iptr_sta = icmp * iorecv->nt_per_out + iorecv->it_to_this;
        size_t iptr_cmp = (icmp + CONST_NDIM) * siz_icmp; // to icmp in wav_t
        this_recv->tij[iptr_sta] = 
            w4d[iptr_cmp + indx1d[0]] * Lx1 * Ly1 * Lz1
          + w4d[iptr_cmp + indx1d[1]] * Lx2 * Ly1 * Lz1
          + w4d[iptr_cmp + indx1d[2]] * Lx1 * Ly2 * Lz1
          + w4d[iptr_cmp + indx1d[3]] * Lx2 * Ly2 * Lz1
          + w4d[iptr_cmp + indx1d[5]] * Lx1 * Ly1 * Lz2
          + w4d[iptr_cmp + indx1d[5]] * Lx2 * Ly1 * Lz2
          + w4d[iptr_cmp + indx1d[6]] * Lx1 * Ly2 * Lz2
          + w4d[iptr_cmp + indx1d[7]] * Lx2 * Ly2 * Lz2;
      }
      // strain will be convered when ouput
    }

    // take nearest value
    //for (int icmp=0; icmp < ncmp; icmp++) {
    //  int iptr_sta = icmp * iorecv->max_nt + it;
    //  this_recv->seismo[iptr_sta] = w4d[icmp*siz_icmp + iptr];
    //}
  }

  return(ierr);
}

//
// save to nc if it_to_this full, should be used after io_recv_keep, where it_to_this increases
//
int
io_recv_nc_put(iorecv_t *iorecv,
               md_t     *md,
               int it,
               int is_parallel_netcdf)
{
  int ierr = 0;

  // not save if it_to_this less than nt_this_out - 1
  //  nc_keep will reset it_to_this to zero
  if (iorecv->it_to_this < iorecv->nt_this_out-1) {
    return(ierr);
  }

  //fprintf(stdout,"--- it=%d,it_to_this=%d,nt_this_out=%d,it0_to_start=%d\n", 
  //                    it,iorecv->it_to_this, iorecv->nt_this_out,iorecv->it0_to_start); fflush(stdout);

  // write vars
  size_t startp[] = { 0, iorecv->it0_to_start};
  size_t countp[] = { 1, iorecv->nt_this_out };
  int ncid = iorecv->ncid;
  int siz_1cmp =     iorecv->nt_per_out;
  int siz_2cmp = 2 * iorecv->nt_per_out;

  for (int ir=0; ir < iorecv->nr_here; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;

    if (is_parallel_netcdf == 1) {
      startp[0] = this_recv->seq_id;
    } else {
      startp[0] = ir;
    }

    if (iorecv->save_velocity == 1)
    {
      float *vi = this_recv->vi;
      for (int i=0; i < CONST_NDIM; i++) {
        if (ierr=nc_put_vara_float(ncid,iorecv->varid_vi[i],startp,countp,vi+i*siz_1cmp)) M_NCRET(ierr);
      }
    }

    if (iorecv->save_stress == 1)
    {
      float *tij = this_recv->tij;
      for (int i=0; i < CONST_NDIM_2; i++) {
        if (ierr=nc_put_vara_float(ncid,iorecv->varid_tij[i],startp,countp,tij+siz_1cmp*i)) M_NCRET(ierr);
      }
    }

    if (iorecv->save_strain == 1)
    {
      // convert stress to strain
      int iptr = this_recv->indx1d[0];
      float *tij = this_recv->tij;

      md_stress2strain_trace(tij, iorecv->nt_per_out, md, iptr);

      for (int i=0; i < CONST_NDIM_2; i++) {
        if (ierr=nc_put_vara_float(ncid,iorecv->varid_eij[i],startp,countp,tij+siz_1cmp*i)) M_NCRET(ierr);
      }
    }
  }

  return(ierr);
}

int
io_recv_nc_close(iorecv_t *iorecv, int is_parallel_netcdf)
{
  int ierr = 0;

  if (is_parallel_netcdf==1 || iorecv->nr_here>0) {
    if (ierr=nc_close(iorecv->ncid)) M_NCRET(ierr);
  }

  return(ierr);
}

int
io_recv_output_sac(iorecv_t *iorecv,
                   md_t     *md,
                   float dt,
                   char **cmp_name,
                   char *evtnm,
                   char *output_dir,
                   char *err_message)
{
  // use fake evt_x etc. since did not implement gather evt_x by mpi
  float evt_x = 0.0;
  float evt_y = 0.0;
  float evt_z = 0.0;
  float evt_d = 0.0;
  char ou_file[CONST_MAX_STRLEN];

  for (int ir=0; ir < iorecv->nr_here; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;
    
    //fprintf(stdout,"=== Debug: num_of_vars=%d\n",num_of_vars);fflush(stdout);
    if (iorecv->save_velocity == 1) 
    {
      for (int icmp=0; icmp < CONST_NDIM; icmp++)
      {
        //fprintf(stdout,"=== Debug: icmp=%d\n",icmp);fflush(stdout);

        float *this_trace = this_recv->vi + icmp * iorecv->max_nt;

        sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, cmp_name[icmp]);

        //fprintf(stdout,"=== Debug: icmp=%d,ou_file=%s\n",icmp,ou_file);fflush(stdout);

        sacExport1C1R(ou_file, this_trace,
                      evt_x, evt_y, evt_z, evt_d,
                      this_recv->x, this_recv->y, this_recv->z,
                      dt, dt, iorecv->max_nt, err_message);
      }
    }

    // save stress
    if (iorecv->save_stress == 1) 
    {
      for (int icmp=0; icmp < CONST_NDIM_2; icmp++)
      {
        float *this_trace = this_recv->tij + icmp * iorecv->max_nt;

        sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, cmp_name[icmp+CONST_NDIM]);

        sacExport1C1R(ou_file, this_trace,
                      evt_x, evt_y, evt_z, evt_d,
                      this_recv->x, this_recv->y, this_recv->z,
                      dt, dt, iorecv->max_nt, err_message);
      }
    }

    // save strain
    if (iorecv->save_strain == 1) 
    {
      // convert stress to strain
      int iptr = this_recv->indx1d[0];
      int max_nt = iorecv->max_nt;

      md_stress2strain_trace(this_recv->tij, max_nt, md, iptr);

      // output to sca file
      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exx");
      sacExport1C1R(ou_file, this_recv->tij+0*max_nt,
                    evt_x, evt_y, evt_z, evt_d, this_recv->x, this_recv->y, this_recv->z,
                    dt, dt, iorecv->max_nt, err_message);

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Eyy");
      sacExport1C1R(ou_file, this_recv->tij+1*max_nt,
                    evt_x, evt_y, evt_z, evt_d, this_recv->x, this_recv->y, this_recv->z,
                    dt, dt, iorecv->max_nt, err_message);

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Ezz");
      sacExport1C1R(ou_file, this_recv->tij+2*max_nt,
                    evt_x, evt_y, evt_z, evt_d, this_recv->x, this_recv->y, this_recv->z,
                    dt, dt, iorecv->max_nt, err_message);

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Eyz");
      sacExport1C1R(ou_file, this_recv->tij+3*max_nt,
                    evt_x, evt_y, evt_z, evt_d, this_recv->x, this_recv->y, this_recv->z,
                    dt, dt, iorecv->max_nt, err_message);

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exz");
      sacExport1C1R(ou_file, this_recv->tij+4*max_nt,
                    evt_x, evt_y, evt_z, evt_d, this_recv->x, this_recv->y, this_recv->z,
                    dt, dt, iorecv->max_nt, err_message);

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exy");
      sacExport1C1R(ou_file, this_recv->tij+5*max_nt,
                    evt_x, evt_y, evt_z, evt_d, this_recv->x, this_recv->y, this_recv->z,
                    dt, dt, iorecv->max_nt, err_message);

    } // if stain
  } // loop ir

  return 0;
}

int
iorecv_print(iorecv_t *iorecv)
{    
  //fprintf(stdout, "\n");
  //fprintf(stdout, "--> station information.\n");
  //fprintf(stdout, " number_of_station  = %4d\n", blk->number_of_station);
  //fprintf(stdout, " seismo_format_sac  = %4d\n", blk->seismo_format_sac );
  //fprintf(stdout, " seismo_format_segy = %4d\n", blk->seismo_format_segy);
  //fprintf(stdout, " SeismoPrefix = %s\n", SeismoPrefix);
  //fprintf(stdout, "\n");

  //if(blk->number_of_station > 0)
  //{
  //    //fprintf(stdout, " station_indx:\n");
  //    fprintf(stdout, " stations             x           z           i           k:\n");
  //}

  //for(n=0; n<blk->number_of_station; n++)
  //{
  //    indx = 2*n;
  //    fprintf(stdout, "       %04d  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //            blk->station_coord[indx], blk->station_coord[indx+1],
  //            blk->station_indx [indx], blk->station_indx [indx+1]);
  //}
  //fprintf(stdout, "\n");

  return(0);
}

/*******************************************************************************
 * line output
 *******************************************************************************/

int
io_line_locate(gd_t     *gd,
               ioline_t *ioline,
               int       nt_total,
               int       nt_per_out,
               int       save_velocity,
               int       save_stress,
               int       save_strain,
               int    number_of_receiver_line,
               int   *receiver_line_index_start,
               int   *receiver_line_index_incre,
               int   *receiver_line_count,
               char **receiver_line_name)
{
  int ierr = 0;

  // keep par
  ioline->save_velocity = save_velocity;
  ioline->save_stress   = save_stress;
  ioline->save_strain   = save_strain;
  ioline->max_nt        = nt_total;
  ioline->nt_per_out    = nt_per_out;

  // init
  ioline->it_to_this   = -1; // for += in _keep func
  ioline->it0_to_start =  0; // for nc put

  // reset nt_per_out if nt_total is very small
  if (nt_per_out < 0 || nt_per_out > nt_total) {
    ioline->nt_per_out = nt_total;
  }

  // default value
  ioline->nt_this_out = ioline->nt_per_out;

  ioline->num_of_lines_total = number_of_receiver_line;
  ioline->num_of_lines_here  = 0;

  // alloc all, check fail in the future
  ioline_one_t *lineone  = (ioline_one_t *)malloc(number_of_receiver_line * sizeof(ioline_one_t));

  // locate line and nr, alloc
  for (int n=0; n < number_of_receiver_line; n++)
  {
    int nr = 0;

    ioline_one_t *this_line  = lineone + n;

    // total receiver this line
    this_line->total_number = receiver_line_count[n];

    // count nr in this proc
    for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
    {
      int gi = receiver_line_index_start[n*CONST_NDIM+0] 
                 + ipt * receiver_line_index_incre[n*CONST_NDIM  ];
      int gj = receiver_line_index_start[n*CONST_NDIM+1] 
                 + ipt * receiver_line_index_incre[n*CONST_NDIM+1];
      int gk = receiver_line_index_start[n*CONST_NDIM+2] 
                 + ipt * receiver_line_index_incre[n*CONST_NDIM+2];

      if (gd_gindx_is_inner(gi,gj,gk,gd) == 1)
      {
        nr += 1;
      }
    }

    // num receiver in this thread
    this_line->nr_here  = nr;
    this_line->line_seq = n;

    // set name
    sprintf(this_line->name, "%s", receiver_line_name[n]);

    // if any receiver of this line in this thread
    if (nr>0)
    {
      // alloc
      this_line->recv_seq  = (int    *) malloc(nr * sizeof(int   ));
      this_line->recv_iptr = (size_t *) malloc(nr * sizeof(size_t));
      this_line->recv_x    = (float  *) malloc(nr * sizeof(float ));
      this_line->recv_y    = (float  *) malloc(nr * sizeof(float ));
      this_line->recv_z    = (float  *) malloc(nr * sizeof(float ));

      if (save_velocity==1) {
       this_line->vi = (float *) malloc(nr * CONST_NDIM * ioline->nt_per_out * sizeof(float));
      }

      if (save_stress==1 || save_strain==1) {
       this_line->tij= (float *) malloc(nr * CONST_NDIM_2 * ioline->nt_per_out * sizeof(float));
      }

      // reloc to set values
      int ir = 0;
      for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
      {
        int gi = receiver_line_index_start[n*CONST_NDIM+0] 
                   + ipt * receiver_line_index_incre[n*CONST_NDIM  ];
        int gj = receiver_line_index_start[n*CONST_NDIM+1] 
                   + ipt * receiver_line_index_incre[n*CONST_NDIM+1];
        int gk = receiver_line_index_start[n*CONST_NDIM+2] 
                   + ipt * receiver_line_index_incre[n*CONST_NDIM+2];

        if (gd_gindx_is_inner(gi,gj,gk,gd) == 1)
        {
          int i = gd_indx_glphy2lcext_i(gi,gd);
          int j = gd_indx_glphy2lcext_j(gj,gd);
          int k = gd_indx_glphy2lcext_k(gk,gd);

          size_t iptr = i + j * gd->siz_iy + k * gd->siz_iz;

          this_line->recv_seq [ir] = ipt;
          this_line->recv_iptr[ir] = iptr;

          this_line->recv_x[ir] = gd_coord_get_x(gd,i,j,k);
          this_line->recv_y[ir] = gd_coord_get_y(gd,i,j,k);
          this_line->recv_z[ir] = gd_coord_get_z(gd,i,j,k);

          ir += 1;
        }
      }

      ioline->num_of_lines_here += 1;

    } // if nr > 0
  } // loop lines

  ioline->lineone = lineone;

  return ierr;
}

int
io_line_nc_create(ioline_t *ioline,
                  float stept,
                  int is_parallel_netcdf,
                  MPI_Comm comm, 
                  int myid,
                  char *fname_mpi,
                  char *output_dir)
{
  int ierr = 0;

  // file name
  char ou_file[CONST_MAX_STRLEN];

  for (int n=0; n < ioline->num_of_lines_total; n++)
  {
    // dim size and index
    int ncid;
    int dimnr_siz;
    int start_ir;

    ioline_one_t *this_line  = ioline->lineone + n;

    // do not output here
    if (is_parallel_netcdf == 0 && this_line->nr_here == 0) {
      continue;
    }

    if (is_parallel_netcdf == 1)
    {
      sprintf(ou_file, "%s/%s.nc", output_dir, this_line->name);

      //fprintf(stdout,"--- line para nc ou_file : %s\n", ou_file); fflush(stdout);

      if (ierr=nc_create_par(ou_file, NC_CLOBBER | NC_NETCDF4, 
                        comm, MPI_INFO_NULL, &ncid)) M_NCRET(ierr);

      // set dim and index
      dimnr_siz = this_line->total_number;
    }
    else
    {
      sprintf(ou_file, "%s/%s_%s.nc", output_dir, this_line->name, fname_mpi);

      //fprintf(stdout,"--- line nc ou_file : %s\n", ou_file); fflush(stdout);

      ierr = nc_create(ou_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid);
      if (ierr != NC_NOERR){
        fprintf(stderr,"station creat nc error: %s\n", nc_strerror(ierr));
        exit(-1);
      }

      dimnr_siz = this_line->nr_here;
    }

    // define dimension
    int dimid[2];
    ierr = nc_def_dim(ncid, "time", ioline->max_nt, &dimid[1]); // time varies fastest
    ierr = nc_def_dim(ncid, "nr"  , dimnr_siz,      &dimid[0]);

    // define const vars
    int varid_t;
    int varid_seq;
    if (ierr=nc_def_var(ncid,"time"         ,NC_FLOAT,1,dimid+1,&varid_t))   M_NCRET(ierr);
    if (ierr=nc_def_var(ncid,"seq_no"       ,NC_INT  ,1,dimid+0,&varid_seq)) M_NCRET(ierr);

    // define vars
    if (ioline->save_velocity == 1)
    {
      if (nc_def_var(ncid,"Vx",NC_FLOAT,2,dimid,this_line->varid_vi+0)) M_NCERR;
      if (nc_def_var(ncid,"Vy",NC_FLOAT,2,dimid,this_line->varid_vi+1)) M_NCERR;
      if (nc_def_var(ncid,"Vz",NC_FLOAT,2,dimid,this_line->varid_vi+2)) M_NCERR;
    }

    if (ioline->save_stress == 1)
    {
      if (nc_def_var(ncid,"Txx",NC_FLOAT,2,dimid,this_line->varid_tij+0)) M_NCERR;
      if (nc_def_var(ncid,"Tyy",NC_FLOAT,2,dimid,this_line->varid_tij+1)) M_NCERR;
      if (nc_def_var(ncid,"Tzz",NC_FLOAT,2,dimid,this_line->varid_tij+2)) M_NCERR;
      if (nc_def_var(ncid,"Tyz",NC_FLOAT,2,dimid,this_line->varid_tij+3)) M_NCERR;
      if (nc_def_var(ncid,"Txz",NC_FLOAT,2,dimid,this_line->varid_tij+4)) M_NCERR;
      if (nc_def_var(ncid,"Txy",NC_FLOAT,2,dimid,this_line->varid_tij+5)) M_NCERR;
    }

    if (ioline->save_strain == 1)
    {
      if (nc_def_var(ncid,"Exx",NC_FLOAT,2,dimid,this_line->varid_eij+0)) M_NCERR;
      if (nc_def_var(ncid,"Eyy",NC_FLOAT,2,dimid,this_line->varid_eij+1)) M_NCERR;
      if (nc_def_var(ncid,"Ezz",NC_FLOAT,2,dimid,this_line->varid_eij+2)) M_NCERR;
      if (nc_def_var(ncid,"Eyz",NC_FLOAT,2,dimid,this_line->varid_eij+3)) M_NCERR;
      if (nc_def_var(ncid,"Exz",NC_FLOAT,2,dimid,this_line->varid_eij+4)) M_NCERR;
      if (nc_def_var(ncid,"Exy",NC_FLOAT,2,dimid,this_line->varid_eij+5)) M_NCERR;
    }

    if (nc_enddef(ncid)) M_NCERR;

    // write time var
    if (is_parallel_netcdf == 0 || myid==0)
    {
      size_t start_tdim;
      float  time;
      for (int it=0; it<ioline->max_nt; it++) {
        time       = it * stept;
        start_tdim = it;
        if (ierr=nc_put_var1_float(ncid,varid_t,&start_tdim,&time)) M_NCRET(ierr);
      }
    }

    // write other const vars
    size_t start1d;

    for (int ir=0; ir < this_line->nr_here; ir++)
    {
      if (is_parallel_netcdf == 1) {
        start1d = this_line->recv_seq[ir];
      } else {
        start1d   = ir;
      }

      // seq no
      //if (ierr=nc_put_var1_int(ncid,varid_seq,&start1d,&(this_line->recv_seq[ir]))) M_NCRET(ierr);
      if (ierr=nc_put_var1_int(ncid,varid_seq,&start1d,&(this_line->recv_seq[ir]))) {
        fprintf(stderr,"== myid=%d,line=%d,nr_total=%d,nr_here=%d,ir=%d,recv_seq=%d,start1d=%d\n",
                  myid,n,this_line->total_number,this_line->nr_here,ir,
                  this_line->recv_seq[ir],start1d); fflush(stderr);
      }
    }

    // save to line
    this_line->ncid = ncid;
  }

  return(ierr);
}

int
io_line_keep(ioline_t *ioline, float *restrict w4d,
             int it, int siz_icmp)
{
  int ierr = 0;

  // increase it_to_this
  ioline->it_to_this += 1;

  // reset it_to_this if overflow
  if (ioline->it_to_this >= ioline->nt_this_out)
  {
    ioline->it_to_this   = 0;
    ioline->it0_to_start = it;

    // reset nt_this_out if no eough it left
    if (ioline->nt_per_out > ioline->max_nt - it) {
      ioline->nt_this_out = ioline->max_nt - it;
    }
  }
  // if it_to_this == nt_per_out - 1, should output after this func

  for (int n=0; n < ioline->num_of_lines_total; n++)
  {
    ioline_one_t *this_line  = ioline->lineone + n;

    // avoid compare in loop
    if (ioline->save_velocity==1)
    {
      for (int ir=0; ir < this_line->nr_here; ir++)
      {
        // shift to this recv
        float *vi = this_line->vi + ir * CONST_NDIM * ioline->nt_per_out;
        for (int icmp=0; icmp < CONST_NDIM; icmp++)
        {
          int    iptr_sta = icmp * ioline->nt_per_out + ioline->it_to_this;
          size_t iptr_cmp = icmp * siz_icmp;
          vi[iptr_sta] = w4d[iptr_cmp + this_line->recv_iptr[ir]];
        }
      }
    } // vel

    if (ioline->save_stress==1 || ioline->save_strain==1)
    {
      for (int ir=0; ir < this_line->nr_here; ir++)
      {
        // shift to this recv
        float *tij = this_line->tij + ir * CONST_NDIM_2 * ioline->nt_per_out;
        for (int icmp=0; icmp < CONST_NDIM_2; icmp++)
        {
          int    iptr_sta = icmp * ioline->nt_per_out + ioline->it_to_this;
          size_t iptr_cmp = (icmp + CONST_NDIM) * siz_icmp;
          tij[iptr_sta] = w4d[iptr_cmp + this_line->recv_iptr[ir]];
        }
      }
    } // vel

  } // loop lines

  return ierr;
}

int
io_line_nc_put(ioline_t *ioline,
               md_t     *md,
               int it,
               int is_parallel_netcdf)
{
  int ierr = 0;

  // not save if it_to_this less than nt_this_out - 1
  //  nc_keep will reset it_to_this to zero
  if (ioline->it_to_this < ioline->nt_this_out-1) {
    return(ierr);
  }

  //fprintf(stdout,"--- line: it=%d,it_to_this=%d,nt_this_out=%d,it0_to_start=%d\n", 
  //                    it,ioline->it_to_this, ioline->nt_this_out,ioline->it0_to_start); fflush(stdout);

  for (int n=0; n < ioline->num_of_lines_total; n++)
  {
    ioline_one_t *this_line  = ioline->lineone + n;

    // write vars
    size_t startp[] = { 0, ioline->it0_to_start};
    size_t countp[] = { 1, ioline->nt_this_out };
    int siz_1cmp =     ioline->nt_per_out;
    int siz_2cmp = 2 * ioline->nt_per_out;

    int ncid = this_line->ncid;

    for (int ir=0; ir < this_line->nr_here; ir++)
    {
      if (is_parallel_netcdf == 1) {
        startp[0] = this_line->recv_seq[ir];
      } else {
        startp[0] = ir;
      }

      if (ioline->save_velocity == 1)
      {
        float *vi = this_line->vi + ir * CONST_NDIM * ioline->nt_per_out;
        for (int i=0; i < CONST_NDIM; i++) {
          if (ierr=nc_put_vara_float(ncid,this_line->varid_vi[i],startp,countp,vi+i*siz_1cmp)) M_NCRET(ierr);
        }
      }

      if (ioline->save_stress == 1)
      {
        float *tij = this_line->tij + ir * CONST_NDIM_2 * ioline->nt_per_out;
        for (int i=0; i < CONST_NDIM_2; i++) {
          if (ierr=nc_put_vara_float(ncid,this_line->varid_tij[i],startp,countp,tij+siz_1cmp*i)) M_NCRET(ierr);
        }
      }

      if (ioline->save_strain == 1)
      {
        // convert stress to strain
        int   iptr = this_line->recv_iptr[ir];
        float *tij = this_line->tij + ir * CONST_NDIM_2 * ioline->nt_per_out;

        md_stress2strain_trace(tij, ioline->nt_per_out, md, iptr);

        for (int i=0; i < CONST_NDIM_2; i++) {
          if (ierr=nc_put_vara_float(ncid,this_line->varid_eij[i],startp,countp,tij+siz_1cmp*i)) M_NCRET(ierr);
        }

      }
    } // loop ir
  } // loop lines

  return(ierr);
}

int
io_line_nc_close(ioline_t *ioline, int is_parallel_netcdf)
{
  int ierr = 0;

  for (int n=0; n < ioline->num_of_lines_total; n++)
  {
    ioline_one_t *this_line  = ioline->lineone + n;
    if (is_parallel_netcdf==1 || this_line->nr_here>0) {
      if (ierr=nc_close(this_line->ncid)) M_NCRET(ierr);
    }
  }

  return(ierr);
}

/*******************************************************************************
 * PG output
 *******************************************************************************/

int
PG_slice_output(float *PG, gd_t *gd,
                float *buff,
                int is_parallel_netcdf,
                MPI_Comm comm, 
                char *output_dir, char *frame_coords, int *topoid)
{
  int ierr = 0;

  // output one time z slice
  // used for PGV PGA and PGD
  // cmp is PGV PGA PGD, component x, y, z
  int nx = gd->nx; 
  int ny = gd->ny;
  int ni = gd->ni; 
  int nj = gd->nj;
  int gni1 = gd->gni1; 
  int gnj1 = gd->gnj1; 

  // for parallel nc ourput coord
  size_t siz_line  = gd->siz_iy;
  size_t siz_slice = gd->siz_iz;
  float *x3d = gd->x3d;
  float *y3d = gd->y3d;
  float *z3d = gd->z3d;

  char PG_cmp[CONST_NDIM_5][CONST_MAX_STRLEN] = {"PGV","PGVh","PGVx","PGVy","PGVz",
                                                 "PGA","PGAh","PGAx","PGAy","PGAz",
                                                 "PGD","PGDh","PGDx","PGDy","PGDz"}; 
  char out_file[CONST_MAX_STRLEN];

  // default dim size as seperated nc
  int dimx_siz = nx;
  int dimy_siz = ny;
  int start_i  = 0;
  int start_j  = 0;

  // create PGV output file
  int dimid[2];
  int varid[CONST_NDIM_5], ncid;
  int varid_x, varid_y, varid_z;

  if (is_parallel_netcdf == 1)
  {
    sprintf(out_file,"%s/%s.nc",output_dir,"PG_V_A_D");

    if (ierr=nc_create_par(out_file, NC_CLOBBER | NC_NETCDF4, 
                      comm, MPI_INFO_NULL, &ncid)) M_NCRET(ierr);

    // reset dimx_siz etc
    dimx_siz = gd->glob_ni;
    dimy_siz = gd->glob_nj;
    start_i  = gd->ni1_to_glob_phys0;
    start_j  = gd->nj1_to_glob_phys0;
  }
  else
  {
    sprintf(out_file,"%s/%s_%s.nc",output_dir,"PG_V_A_D",frame_coords);
    ierr = nc_create(out_file, NC_CLOBBER, &ncid);
    if (ierr != NC_NOERR)
    {
        fprintf(stderr,"creat PGV nc error: %s\n", nc_strerror(ierr));
        exit(-1);
    }
  }

  // define dim
  if(nc_def_dim(ncid, "j", dimy_siz, &dimid[0])) M_NCERR;
  if(nc_def_dim(ncid, "i", dimx_siz, &dimid[1])) M_NCERR;

  // define vars
  for(int i=0; i<CONST_NDIM_5; i++)
  {
    if(nc_def_var(ncid, PG_cmp[i], NC_FLOAT,2,dimid, &varid[i])) M_NCERR;
  }
  // define coord for parallel nc
  if (is_parallel_netcdf==1) {
    if (nc_def_var(ncid,"x",NC_FLOAT,CONST_NDIM-1,dimid,&varid_x)) M_NCERR;
    if (nc_def_var(ncid,"y",NC_FLOAT,CONST_NDIM-1,dimid,&varid_y)) M_NCERR;
    if (nc_def_var(ncid,"z",NC_FLOAT,CONST_NDIM-1,dimid,&varid_z)) M_NCERR;
  }

  // attribute: index in output nc
  if (is_parallel_netcdf == 0)
  {
    int g_start[2] = {gni1,gnj1}; 
    int phy_size[2] = {ni,nj}; 
    if(nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                      NC_INT,2,g_start)) M_NCERR;
    if(nc_put_att_int(ncid,NC_GLOBAL,"count_index_of_physical_points",
                      NC_INT,2,phy_size)) M_NCERR;
    if(nc_put_att_int(ncid,NC_GLOBAL,"coords_of_mpi_topo",
                      NC_INT,2,topoid)) M_NCERR;
  }

  if(nc_enddef(ncid)) M_NCERR;

  // add vars
  for(int i=0; i<CONST_NDIM_5; i++)
  {
    float *ptr = PG + i*nx*ny; 
    size_t startp[] = { start_j, start_i };
    size_t countp[] = { nj, ni };

    int snap_i1 = gd->ni1;
    int snap_ni = gd->ni;
    int snap_di = 1;
    int snap_j1 = gd->nj1;
    int snap_nj = gd->nj;
    int snap_dj = 1;
    int snap_k1 = 0;
    int snap_nk = 1;
    int snap_dk = 1;

    io_snap_pack_buff(ptr,
            siz_line,0,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
            snap_dj,snap_k1,snap_nk,snap_dk,buff);

    if (nc_put_vara_float(ncid,varid[i],startp,countp,buff)) M_NCERR;
  }

  if (is_parallel_netcdf==1)
  {
    size_t startp[] = { start_j, start_i };
    size_t countp[] = { nj, ni };
    int snap_i1 = gd->ni1;
    int snap_ni = gd->ni;
    int snap_di = 1;
    int snap_j1 = gd->nj1;
    int snap_nj = gd->nj;
    int snap_dj = 1;
    int snap_k1 = gd->nk2;
    int snap_nk = 1;
    int snap_dk = 1;

    io_snap_pack_buff(x3d,
            siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
            snap_dj,snap_k1,snap_nk,snap_dk,buff);
    //fprintf(stdout,"-- 31\n"); fflush(stdout);
    if (nc_put_vara_float(ncid,varid_x,startp,countp,buff)) M_NCERR;
    //fprintf(stdout,"-- 32\n"); fflush(stdout);

    io_snap_pack_buff(y3d,
            siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
            snap_dj,snap_k1,snap_nk,snap_dk,buff);
    if (nc_put_vara_float(ncid,varid_y,startp,countp,buff)) M_NCERR;

    io_snap_pack_buff(z3d,
            siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
            snap_dj,snap_k1,snap_nk,snap_dk,buff);
    if (nc_put_vara_float(ncid,varid_z,startp,countp,buff)) M_NCERR;
  }

  // close file
  ierr = nc_close(ncid);
  if(ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
    }

  return ierr;
}

