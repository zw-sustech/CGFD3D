/*
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "netcdf.h"

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "io_funcs.h"

/*
 * export to a single file
 */

void
io_build_fname(char *out_dir,
               char *prefix,
               char *sufix,
               int *myid2,
               char *ou_fname)
{
  sprintf(ou_fname,"%s/%s_%d_%d%s",out_dir,prefix,myid2[0],myid2[1],sufix);

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
                    int *myid2,
                    int  it,
                    char *ou_fname)
{
  sprintf(ou_fname,"%s/%s_%d_%d_it%d%s",out_dir,prefix,myid2[0],myid2[1],it,sufix);
}

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
  int dimid[FD_NDIM];

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
    ierr = nc_def_var(ncid, v3d_name[ivar], NC_FLOAT, FD_NDIM, dimid, &varid[ivar]);
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
