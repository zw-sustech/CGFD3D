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
//  sta_name = (char **) fdlib_mem_malloc_2l_char(sta_num,FD_MAX_STRLEN,"io_read_station_list");
//  sta_xyz  = (float *) fdlib_mem_malloc_1d(sta_num*FD_NDIM*sizeof(float),"io_read_station_list");
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
//                      sta_xzy+ir*FD_NDIM,
//                      sta_xzy+ir*FD_NDIM+1,
//                      sta_xzy+ir*FD_NDIM+2);
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

/*
 * read in station list file and locate station
 */

int
io_read_locate_station(char *in_filenm, 
                       int   glob_phys_ix1, // gloabl start index along x this thread
                       int   glob_phys_ix2, // gloabl end index along x
                       int   glob_phys_iy1,
                       int   glob_phys_iy2,
                       int   glob_phys_iz1,
                       int   glob_phys_iz2,
                       int   ni1,
                       int   nj1,
                       int   nk1,
                       size_t siz_line,
                       size_t siz_slice,
                       float *x3d, float *y3d, float *z3d,
                       int *num_of_sta,
                       char ***p_sta_name,
                       float **p_sta_coord,
                       int   **p_sta_point,
                       float **p_sta_shift)
{
  FILE *fp;
  char line[500];

  if (!(fp = fopen (in_filenm, "rt")))
	{
	    fprintf (stdout, "Cannot open input file %s\n", in_filenm);
	    fflush (stdout);
	    return(1);
	}

  // number of station
  int nr, nr_loc, nr_point;

  fgets(line,500,fp);
  sscanf(line, "%d %d",&nr_loc, &nr_point);

  nr = nr_loc + nr_point;

  //fprintf(stdout, "-- nr_loc=%d, nr_point=%d, nr=%d\n", nr_loc, nr_point, nr);

  if (nr_loc > 0) {
    fprintf(stderr, "ERROR: nr_loc=%d, coord loc not implement yet\n", nr_loc);
    fflush(stderr);
    exit(1);
  }

  // alloc by maximum number to reduce code
  char **sta_name;
  float *sta_shift;
  float *sta_coord;
  int   *sta_point;

  sta_name = (char **) fdlib_mem_malloc_2l_char(nr,FD_MAX_STRLEN,"io_read_locate_station");
  sta_coord= (float *) fdlib_mem_malloc_1d(nr*FD_NDIM*sizeof(float),"io_read_locate_station");
  sta_shift= (float *) fdlib_mem_malloc_1d(nr*FD_NDIM*sizeof(float),"io_read_locate_station");
  sta_point= (int   *) fdlib_mem_malloc_1d(nr*FD_NDIM*sizeof(int),"io_read_locate_station");

  // read coord and locate

  int ir=0;
  int nr_this = 0; // in this thread

  //for (ir=0; ir<nr_loc; ir++)
  //{
  //  float rx, ry, rz;
  //  // read one line
  //  fgets(line,500,fp);
  //  // get values
  //  sscanf(line, "%s %f %f %f", sta_name[ir], &rx, &ry, &rz);
  //  // locate
  //  if (is_coord_in_phys_region(rx,ry,rz,nx,ny,nz,ni1,ni2,nj1,nj2,nk1,nk2,x3d,y3d,z3d)==1)
  //  {
  //    int ptr_this = nr_this * FD_NDIM;
  //    sprintf(sta_name[nr_this], "%s", sta_name[ir]);
  //    sta_coord[ptr_this+0]=rx;
  //    sta_coord[ptr_this+1]=ry;
  //    sta_coord[ptr_this+2]=rz;
  //    // set point and shift

  //    nr_this += 1;
  //  }
  //}

  // read index and locate

  for (ir=0; ir<nr_point; ir++)
  {
    int ix, iy, iz;
    int i_local, j_local, k_local;
    // read one line
    fgets(line,500,fp);
    // get values
    sscanf(line, "%s %d %d %d", sta_name[ir], &ix, &iy, &iz);
    //fprintf(stdout,"== in: %s %d %d %d\n", sta_name[ir],ix,iy,iz); fflush(stdout);
    // locate
    int ptr_this = nr_this * FD_NDIM;
    if ( ix >= glob_phys_ix1 && ix <= glob_phys_ix2 &&
         iy >= glob_phys_iy1 && iy <= glob_phys_iy2 &&
         iz >= glob_phys_iz1 && iz <= glob_phys_iz2 )
    {
      // convert to local index w ghost
      int i_local = ix - glob_phys_ix1 + ni1;
      int j_local = iy - glob_phys_iy1 + nj1;
      int k_local = iz - glob_phys_iz1 + nk1;

      int ptr_this = nr_this * FD_NDIM;
      int ptr_point = i_local + j_local * siz_line + k_local * siz_slice;

      sprintf(sta_name[nr_this], "%s", sta_name[ir]);
      // get coord
      sta_coord[ptr_this+0]=x3d[ptr_point];
      sta_coord[ptr_this+1]=y3d[ptr_point];
      sta_coord[ptr_this+2]=z3d[ptr_point];
      // set point and shift
      sta_point[ptr_this+0]=i_local;
      sta_point[ptr_this+1]=j_local;
      sta_point[ptr_this+2]=k_local;
      sta_shift[ptr_this+0]=0.0;
      sta_shift[ptr_this+1]=0.0;
      sta_shift[ptr_this+2]=0.0;

      //fprintf(stdout,"== ir_this=%d,name=%s,i=%d,j=%d,k=%d\n",
      //      nr_this,sta_name[nr_this],i_local,j_local,k_local); fflush(stdout);

      nr_this += 1;
    }
  }

  fclose(fp);

  *num_of_sta  = nr_this;
  *p_sta_name = sta_name;
  *p_sta_coord  = sta_coord;
  *p_sta_point  = sta_point;
  *p_sta_shift  = sta_shift;

  return(0);
}
