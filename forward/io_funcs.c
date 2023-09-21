/*
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "netcdf.h"
#include "sacLib.h"

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "constants.h"
#include "fd_t.h"
#include "gd_t.h"
#include "io_funcs.h"

//#define M_NCERR(ierr) {fprintf(stderr,"sv_ nc error: %s\n", nc_strerror(ierr)); exit(1);}
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

/*
 * read in station list file and locate station
 */

int
io_recv_read_locate(
                    gd_t *gd,
                    iorecv_t  *iorecv,
                    int       nt_total,
                    int       num_of_vars,
		    int       medium_type,
                    char *in_filenm,
                    MPI_Comm  comm,
                    int       myid,
                    int       verbose)
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
  int nr;

  io_get_nextline(fp, line, 500);
  sscanf(line, "%d", &nr);

  //fprintf(stdout, "-- nr=%d\n", nr);

  // check fail in the future
  iorecv_one_t *recvone = (iorecv_one_t *)malloc(nr * sizeof(iorecv_one_t));

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

      int ptr_this = nr_this * CONST_NDIM;
      iorecv_one_t *this_recv = recvone + nr_this;

      sprintf(this_recv->name, "%s", recvone[ir].name);

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

  iorecv->total_number = nr_this;
  iorecv->recvone      = recvone;
  iorecv->max_nt       = nt_total;
  iorecv->ncmp         = num_of_vars;
  if(medium_type == CONST_MEDIUM_VISCOELASTIC_ISO)
  {
    iorecv->ncmp = 9;//not output J
  }

  // malloc seismo
  for (int ir=0; ir < iorecv->total_number; ir++)
  {
    recvone = iorecv->recvone + ir;
    recvone->seismo = (float *) malloc(num_of_vars * nt_total * sizeof(float));
  }

  return(0);
}

int
io_line_locate(
               gd_t *gd,
               ioline_t *ioline,
               int    num_of_vars,
	       int    medium_type,
               int    nt_total,
               int    number_of_receiver_line,
               int   *receiver_line_index_start,
               int   *receiver_line_index_incre,
               int   *receiver_line_count,
               char **receiver_line_name)
{
  int ierr = 0;

  // init
  ioline->num_of_lines  = 0;
  ioline->max_nt        = nt_total;
  ioline->ncmp          = num_of_vars;
  if(medium_type == CONST_MEDIUM_VISCOELASTIC_ISO)
  {
    ioline->ncmp = 9;//not output J
  }

  // alloc as max num to keep nr and seq values, easy for second round
  ioline->line_nr  = (int *) malloc(number_of_receiver_line * sizeof(int));
  ioline->line_seq = (int *) malloc(number_of_receiver_line * sizeof(int));

  // first run to count line and nr
  for (int n=0; n < number_of_receiver_line; n++)
  {
    int nr = 0;
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

    // if any receiver of this line in this thread
    if (nr>0)
    {
      ioline->line_nr [ ioline->num_of_lines ] = nr;
      ioline->line_seq[ ioline->num_of_lines ] = n;
      ioline->num_of_lines += 1;
    }
  }

  // alloc
  if (ioline->num_of_lines>0)
  {
    ioline->line_name   = (char **)fdlib_mem_malloc_2l_char(ioline->num_of_lines,
                                    CONST_MAX_STRLEN, "io_line_locate");

    ioline->recv_seq    = (int **) malloc(ioline->num_of_lines * sizeof(int*));
    //ioline->recv_indx   = (int **) malloc(ioline->num_of_lines * sizeof(int*));
    ioline->recv_iptr   = (int **) malloc(ioline->num_of_lines * sizeof(int*));
    ioline->recv_x  = (float **) malloc(ioline->num_of_lines * sizeof(float*));
    ioline->recv_y  = (float **) malloc(ioline->num_of_lines * sizeof(float*));
    ioline->recv_z  = (float **) malloc(ioline->num_of_lines * sizeof(float*));
    ioline->recv_seismo = (float **) malloc(ioline->num_of_lines * sizeof(float*));

    for (int n=0; n < ioline->num_of_lines; n++)
    {
      int nr = ioline->line_nr[n];
      //ioline->recv_indx[n] = (int *)malloc(nr * CONST_NDIM * sizeof(int)); 
      ioline->recv_seq [n]  = (int *)malloc( nr * sizeof(int) ); 
      ioline->recv_iptr[n]  = (int *)malloc( nr * sizeof(int) ); 
      ioline->recv_x[n] = (float *)malloc( nr * sizeof(float) );
      ioline->recv_y[n] = (float *)malloc( nr * sizeof(float) );
      ioline->recv_z[n] = (float *)malloc( nr * sizeof(float) );
      ioline->recv_seismo[n] = (float *)malloc(
                                nr * num_of_vars * nt_total * sizeof(float) );
    }
  }

  // second run for value
  //  only loop lines in this thread
  for (int m=0; m < ioline->num_of_lines; m++)
  {
    int n = ioline->line_seq[m];

    sprintf(ioline->line_name[m], "%s", receiver_line_name[n]);

    int ir = 0;
    for (int ipt=0; ipt<receiver_line_count[n]; ipt++)
    {
      int gi = receiver_line_index_start[n*CONST_NDIM+0] + ipt * receiver_line_index_incre[n*CONST_NDIM  ];
      int gj = receiver_line_index_start[n*CONST_NDIM+1] + ipt * receiver_line_index_incre[n*CONST_NDIM+1];
      int gk = receiver_line_index_start[n*CONST_NDIM+2] + ipt * receiver_line_index_incre[n*CONST_NDIM+2];

      if (gd_gindx_is_inner(gi,gj,gk,gd) == 1)
      {
        int i = gd_indx_glphy2lcext_i(gi,gd);
        int j = gd_indx_glphy2lcext_j(gj,gd);
        int k = gd_indx_glphy2lcext_k(gk,gd);

        int iptr = i + j * gd->siz_iy + k * gd->siz_iz;

        ioline->recv_seq [m][ir] = ipt;
        ioline->recv_iptr[m][ir] = iptr;

        ioline->recv_x[m][ir] = gd_coord_get_x(gd,i,j,k);
        ioline->recv_y[m][ir] = gd_coord_get_y(gd,i,j,k);
        ioline->recv_z[m][ir] = gd_coord_get_z(gd,i,j,k);

        ir += 1;
      }
    }
  }

  return ierr;
}

int
io_slice_locate(gd_t  *gdinfo,
                ioslice_t *ioslice,
                int  number_of_slice_x,
                int  number_of_slice_y,
                int  number_of_slice_z,
                int *slice_x_index,
                int *slice_y_index,
                int *slice_z_index,
                char *output_fname_part,
                char *output_dir)
{
  int ierr = 0;

  ioslice->siz_max_wrk = 0;

  if (number_of_slice_x>0) {
    ioslice->slice_x_fname = (char **) fdlib_mem_malloc_2l_char(number_of_slice_x,
                                                            CONST_MAX_STRLEN,
                                                            "slice_x_fname");
    ioslice->slice_x_indx = (int *) malloc(number_of_slice_x * sizeof(int));
  }
  if (number_of_slice_y>0) {
    ioslice->slice_y_fname = (char **) fdlib_mem_malloc_2l_char(number_of_slice_y,
                                                            CONST_MAX_STRLEN,
                                                            "slice_y_fname");
    ioslice->slice_y_indx = (int *) malloc(number_of_slice_y * sizeof(int));
  }
  if (number_of_slice_z>0) {
    ioslice->slice_z_fname = (char **) fdlib_mem_malloc_2l_char(number_of_slice_z,
                                                            CONST_MAX_STRLEN,
                                                            "slice_z_fname");
    ioslice->slice_z_indx = (int *) malloc(number_of_slice_z * sizeof(int));
  }

  // init
  ioslice->num_of_slice_x = 0;
  ioslice->num_of_slice_y = 0;
  ioslice->num_of_slice_z = 0;

  // x slice
  for (int n=0; n < number_of_slice_x; n++)
  {
    int gi = slice_x_index[n];
    if (gd_gindx_is_inner_i(gi, gdinfo)==1)
    {
      int islc = ioslice->num_of_slice_x;

      ioslice->slice_x_indx[islc]  = gd_indx_glphy2lcext_i(gi, gdinfo);
      sprintf(ioslice->slice_x_fname[islc],"%s/slicex_i%d_%s.nc",
                output_dir,gi,output_fname_part);

      ioslice->num_of_slice_x += 1;

      size_t slice_siz = gdinfo->nj * gdinfo->nk;
      ioslice->siz_max_wrk = slice_siz > ioslice->siz_max_wrk ? 
                             slice_siz : ioslice->siz_max_wrk;
    }
  }

  // y slice
  for (int n=0; n < number_of_slice_y; n++)
  {
    int gj = slice_y_index[n];
    if (gd_gindx_is_inner_j(gj, gdinfo)==1)
    {
      int islc = ioslice->num_of_slice_y;

      ioslice->slice_y_indx[islc]  = gd_indx_glphy2lcext_j(gj, gdinfo);
      sprintf(ioslice->slice_y_fname[islc],"%s/slicey_j%d_%s.nc",
                output_dir,gj,output_fname_part);

      ioslice->num_of_slice_y += 1;

      size_t slice_siz = gdinfo->ni * gdinfo->nk;
      ioslice->siz_max_wrk = slice_siz > ioslice->siz_max_wrk ? 
                             slice_siz : ioslice->siz_max_wrk;
    }
  }

  // z slice
  for (int n=0; n < number_of_slice_z; n++)
  {
    int gk = slice_z_index[n];
    if (gd_gindx_is_inner_k(gk, gdinfo)==1)
    {
      int islc = ioslice->num_of_slice_z;

      ioslice->slice_z_indx[islc]  = gd_indx_glphy2lcext_k(gk, gdinfo);
      sprintf(ioslice->slice_z_fname[islc],"%s/slicez_k%d_%s.nc",
                output_dir,gk,output_fname_part);

      ioslice->num_of_slice_z += 1;

      size_t slice_siz = gdinfo->ni * gdinfo->nj;
      ioslice->siz_max_wrk = slice_siz > ioslice->siz_max_wrk ? 
                             slice_siz : ioslice->siz_max_wrk;
    }
  }

  return ierr;
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
                    char *output_fname_part,
                    char *output_dir)
{
  // malloc to max, num of snap will not be large
  if (number_of_snapshot > 0)
  {
    iosnap->fname = (char **) fdlib_mem_malloc_2l_char(number_of_snapshot,
                                    CONST_MAX_STRLEN,"snap_fname");
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

    iosnap->i1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->j1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
    iosnap->k1_to_glob = (int *) malloc(number_of_snapshot * sizeof(int));
  }

  // init

  iosnap->siz_max_wrk = 0;

  int isnap = 0;

  for (int n=0; n < number_of_snapshot; n++)
  {
    int iptr0 = n*CONST_NDIM;

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
      iosnap->i1[isnap]  = gd_indx_glphy2lcext_i(gi1, gdinfo);
      iosnap->j1[isnap]  = gd_indx_glphy2lcext_j(gj1, gdinfo);
      iosnap->k1[isnap]  = gd_indx_glphy2lcext_k(gk1, gdinfo);
      iosnap->ni[isnap]  = ngi;
      iosnap->nj[isnap]  = ngj;
      iosnap->nk[isnap]  = ngk;
      iosnap->di[isnap]  = snapshot_index_incre[iptr0+0];
      iosnap->dj[isnap]  = snapshot_index_incre[iptr0+1];
      iosnap->dk[isnap]  = snapshot_index_incre[iptr0+2];

      iosnap->it1[isnap]  = snapshot_time_start[n];
      iosnap->dit[isnap]  = snapshot_time_incre[n];

      iosnap->out_vel   [isnap] = snapshot_save_velocity[n];
      iosnap->out_stress[isnap] = snapshot_save_stress[n];
      iosnap->out_strain[isnap] = snapshot_save_strain[n];

      iosnap->i1_to_glob[isnap] = i_in_nc;
      iosnap->j1_to_glob[isnap] = j_in_nc;
      iosnap->k1_to_glob[isnap] = k_in_nc;

      sprintf(iosnap->fname[isnap],"%s/%s_%s.nc",output_dir,
                                                 snapshot_name[n],
                                                 output_fname_part);

      // for max wrk
      size_t snap_siz =  ngi * ngj * ngk;
      iosnap->siz_max_wrk = snap_siz > iosnap->siz_max_wrk ? 
                            snap_siz : iosnap->siz_max_wrk;

      isnap += 1;
    } // if in this
  } // loop all snap

  iosnap->num_of_snap = isnap;
}

/*
 * combine init and creat to reduce funcs call
 */

int
io_slice_nc_create(ioslice_t *ioslice, 
                  int num_of_vars, int visco_type, char **w3d_name,
                  int ni, int nj, int nk,
                  int *topoid, ioslice_nc_t *ioslice_nc)
{
  int ierr = 0;

  int num_of_slice_x = ioslice->num_of_slice_x;
  int num_of_slice_y = ioslice->num_of_slice_y;
  int num_of_slice_z = ioslice->num_of_slice_z;

  ioslice_nc->num_of_slice_x = num_of_slice_x;
  ioslice_nc->num_of_slice_y = num_of_slice_y;
  ioslice_nc->num_of_slice_z = num_of_slice_z;
  ioslice_nc->num_of_vars    = num_of_vars   ;
  if (visco_type == CONST_VISCO_GMB){
    num_of_vars = 9;
    ioslice_nc->num_of_vars = 9;
  }

  // malloc vars
  ioslice_nc->ncid_slx = (int *)malloc(num_of_slice_x*sizeof(int));
  ioslice_nc->ncid_sly = (int *)malloc(num_of_slice_y*sizeof(int));
  ioslice_nc->ncid_slz = (int *)malloc(num_of_slice_z*sizeof(int));

  ioslice_nc->timeid_slx = (int *)malloc(num_of_slice_x*sizeof(int));
  ioslice_nc->timeid_sly = (int *)malloc(num_of_slice_y*sizeof(int));
  ioslice_nc->timeid_slz = (int *)malloc(num_of_slice_z*sizeof(int));

  ioslice_nc->varid_slx = (int *)malloc(num_of_vars*num_of_slice_x*sizeof(int));
  ioslice_nc->varid_sly = (int *)malloc(num_of_vars*num_of_slice_y*sizeof(int));
  ioslice_nc->varid_slz = (int *)malloc(num_of_vars*num_of_slice_z*sizeof(int));

  // slice x
  for (int n=0; n<num_of_slice_x; n++)
  {
    int dimid[3];
    if (nc_create(ioslice->slice_x_fname[n], NC_CLOBBER,
                  &(ioslice_nc->ncid_slx[n]))) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_slx[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_slx[n], "k"   , nk          , &dimid[1])) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_slx[n], "j"   , nj          , &dimid[2])) M_NCERR;
    // time var
    if (nc_def_var(ioslice_nc->ncid_slx[n], "time", NC_FLOAT, 1, dimid+0,
                   &(ioslice_nc->timeid_slx[n]))) M_NCERR;
    // other vars
    for (int ivar=0; ivar<num_of_vars; ivar++) {
      if (nc_def_var(ioslice_nc->ncid_slx[n], w3d_name[ivar], NC_FLOAT, 3, dimid,
                     &(ioslice_nc->varid_slx[ivar+n*num_of_vars]))) M_NCERR;
    }
    // attribute: index info for plot
    nc_put_att_int(ioslice_nc->ncid_slx[n],NC_GLOBAL,"i_index_with_ghosts_in_this_thread",
                   NC_INT,1,ioslice->slice_x_indx+n);
    nc_put_att_int(ioslice_nc->ncid_slx[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,topoid);
    // end def
    if (nc_enddef(ioslice_nc->ncid_slx[n])) M_NCERR;
  }

  // slice y
  for (int n=0; n<num_of_slice_y; n++)
  {
    int dimid[3];
    if (nc_create(ioslice->slice_y_fname[n], NC_CLOBBER,
                  &(ioslice_nc->ncid_sly[n]))) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_sly[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_sly[n], "k"   , nk          , &dimid[1])) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_sly[n], "i"   , ni          , &dimid[2])) M_NCERR;
    // time var
    if (nc_def_var(ioslice_nc->ncid_sly[n], "time", NC_FLOAT, 1, dimid+0,
                   &(ioslice_nc->timeid_sly[n]))) M_NCERR;
    // other vars
    for (int ivar=0; ivar<num_of_vars; ivar++) {
      if (nc_def_var(ioslice_nc->ncid_sly[n], w3d_name[ivar], NC_FLOAT, 3, dimid,
                     &(ioslice_nc->varid_sly[ivar+n*num_of_vars]))) M_NCERR;
    }
    // attribute: index info for plot
    nc_put_att_int(ioslice_nc->ncid_sly[n],NC_GLOBAL,"j_index_with_ghosts_in_this_thread",
                   NC_INT,1,ioslice->slice_y_indx+n);
    nc_put_att_int(ioslice_nc->ncid_sly[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,topoid);
    // end def
    if (nc_enddef(ioslice_nc->ncid_sly[n])) M_NCERR;
  }

  // slice z
  for (int n=0; n<num_of_slice_z; n++)
  {
    int dimid[3];
    if (nc_create(ioslice->slice_z_fname[n], NC_CLOBBER,
                  &(ioslice_nc->ncid_slz[n]))) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_slz[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_slz[n], "j"   , nj          , &dimid[1])) M_NCERR;
    if (nc_def_dim(ioslice_nc->ncid_slz[n], "i"   , ni          , &dimid[2])) M_NCERR;
    // time var
    if (nc_def_var(ioslice_nc->ncid_slz[n], "time", NC_FLOAT, 1, dimid+0,
                   &(ioslice_nc->timeid_slz[n]))) M_NCERR;
    // other vars
    for (int ivar=0; ivar<num_of_vars; ivar++) {
      if (nc_def_var(ioslice_nc->ncid_slz[n], w3d_name[ivar], NC_FLOAT, 3, dimid,
                     &(ioslice_nc->varid_slz[ivar+n*num_of_vars]))) M_NCERR;
    }
    // attribute: index info for plot
    nc_put_att_int(ioslice_nc->ncid_slz[n],NC_GLOBAL,"k_index_with_ghosts_in_this_thread",
                   NC_INT,1,ioslice->slice_z_indx+n);
    nc_put_att_int(ioslice_nc->ncid_slz[n],NC_GLOBAL,"coords_of_mpi_topo",
                   NC_INT,2,topoid);
    // end def
    if (nc_enddef(ioslice_nc->ncid_slz[n])) M_NCERR;
  }

  return ierr;
}

int
io_snap_nc_create(iosnap_t *iosnap, iosnap_nc_t *iosnap_nc, int *topoid)
{
  int ierr = 0;

  int num_of_snap = iosnap->num_of_snap;
  char **snap_fname = iosnap->fname;

  iosnap_nc->num_of_snap = num_of_snap;
  iosnap_nc->ncid = (int *)malloc(num_of_snap*sizeof(int));
  iosnap_nc->timeid = (int *)malloc(num_of_snap*sizeof(int));

  iosnap_nc->varid_V = (int *)malloc(num_of_snap*CONST_NDIM*sizeof(int));
  iosnap_nc->varid_T = (int *)malloc(num_of_snap*CONST_NDIM_2*sizeof(int));
  iosnap_nc->varid_E = (int *)malloc(num_of_snap*CONST_NDIM_2*sizeof(int));

  // will be used in put step
  iosnap_nc->cur_it = (int *)malloc(num_of_snap*sizeof(int));
  for (int n=0; n<num_of_snap; n++) {
    iosnap_nc->cur_it[n] = 0;
  }

  int *ncid   = iosnap_nc->ncid;
  int *timeid = iosnap_nc->timeid;
  int *varid_V = iosnap_nc->varid_V;
  int *varid_T = iosnap_nc->varid_T;
  int *varid_E = iosnap_nc->varid_E;

  for (int n=0; n<num_of_snap; n++)
  {
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

    if (nc_create(snap_fname[n], NC_CLOBBER, &ncid[n])) M_NCERR;
    if (nc_def_dim(ncid[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid[n], "k", snap_nk     , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid[n], "j", snap_nj     , &dimid[2])) M_NCERR;
    if (nc_def_dim(ncid[n], "i", snap_ni     , &dimid[3])) M_NCERR;
    // time var
    if (nc_def_var(ncid[n], "time", NC_FLOAT, 1, dimid+0, &timeid[n])) M_NCERR;
    // other vars
    if (snap_out_V==1) {
       if (nc_def_var(ncid[n],"Vx",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Vy",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Vz",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+2])) M_NCERR;
    }
    if (snap_out_T==1) {
       if (nc_def_var(ncid[n],"Txx",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Tyy",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Tzz",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+2])) M_NCERR;
       if (nc_def_var(ncid[n],"Txz",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+3])) M_NCERR;
       if (nc_def_var(ncid[n],"Tyz",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+4])) M_NCERR;
       if (nc_def_var(ncid[n],"Txy",NC_FLOAT,4,dimid,&varid_T[n*CONST_NDIM_2+5])) M_NCERR;
    }
    if (snap_out_E==1) {
       if (nc_def_var(ncid[n],"Exx",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Eyy",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Ezz",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+2])) M_NCERR;
       if (nc_def_var(ncid[n],"Exz",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+3])) M_NCERR;
       if (nc_def_var(ncid[n],"Eyz",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+4])) M_NCERR;
       if (nc_def_var(ncid[n],"Exy",NC_FLOAT,4,dimid,&varid_E[n*CONST_NDIM_2+5])) M_NCERR;
    }
    // attribute: index in output snapshot, index w ghost in thread
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

    if (nc_enddef(ncid[n])) M_NCERR;
  } // loop snap

  return ierr;
}

int
io_slice_nc_put(ioslice_t    *ioslice,
                ioslice_nc_t *ioslice_nc,
                gd_t     *gdinfo,
                float *restrict w4d,
                float *restrict buff,
                int   it,
                float time,
                int   i1_cmp,
                int   i2_cmp)
{
  int ierr = 0;

  int   ni1 = gdinfo->ni1;
  int   ni2 = gdinfo->ni2;
  int   nj1 = gdinfo->nj1;
  int   nj2 = gdinfo->nj2;
  int   nk1 = gdinfo->nk1;
  int   nk2 = gdinfo->nk2;
  int   ni  = gdinfo->ni ;
  int   nj  = gdinfo->nj ;
  int   nk  = gdinfo->nk ;
  size_t siz_line = gdinfo->siz_iy;
  size_t siz_slice= gdinfo->siz_iz;
  size_t siz_icmp = gdinfo->siz_icmp;

  int  num_of_vars = ioslice_nc->num_of_vars;

  //-- slice x, 
  for (int n=0; n < ioslice_nc->num_of_slice_x; n++)
  {
    size_t startp[] = { it, 0, 0 };
    size_t countp[] = { 1, nk, nj};
    size_t start_tdim = it;

    nc_put_var1_float(ioslice_nc->ncid_slx[n], ioslice_nc->timeid_slx[n],
                        &start_tdim, &time);

    for (int ivar=i1_cmp; ivar <= i2_cmp; ivar++)
    {
      float *var = w4d + ivar * siz_icmp;
      int iptr_slice = 0;
      for (int k=nk1; k<=nk2; k++) {
        for (int j=nj1; j<=nj2; j++) {
          int i = ioslice->slice_x_indx[n];
          int iptr = i + j * siz_line + k * siz_slice;
          buff[iptr_slice] = var[iptr];
          iptr_slice++;
        }
      }
      nc_put_vara_float(ioslice_nc->ncid_slx[n], 
                        ioslice_nc->varid_slx[n*num_of_vars + ivar],
                        startp, countp, buff);
    }
  }

  // slice y
  for (int n=0; n < ioslice_nc->num_of_slice_y; n++)
  {
    size_t startp[] = { it, 0, 0 };
    size_t countp[] = { 1, nk, ni};
    size_t start_tdim = it;

    nc_put_var1_float(ioslice_nc->ncid_sly[n], ioslice_nc->timeid_sly[n],
                        &start_tdim, &time);

    for (int ivar=i1_cmp; ivar <= i2_cmp; ivar++)
    {
      float *var = w4d + ivar * siz_icmp;
      int iptr_slice = 0;
      for (int k=nk1; k<=nk2; k++) {
        int j = ioslice->slice_y_indx[n];
        for (int i=ni1; i<=ni2; i++) {
          int iptr = i + j * siz_line + k * siz_slice;
          buff[iptr_slice] = var[iptr];
          iptr_slice++;
        }
      }
      nc_put_vara_float(ioslice_nc->ncid_sly[n], 
                        ioslice_nc->varid_sly[n*num_of_vars + ivar],
                        startp, countp, buff);
    }
  }

  // slice z
  for (int n=0; n < ioslice_nc->num_of_slice_z; n++)
  {
    size_t startp[] = { it, 0, 0 };
    size_t countp[] = { 1, nj, ni};
    size_t start_tdim = it;

    nc_put_var1_float(ioslice_nc->ncid_slz[n], ioslice_nc->timeid_slz[n],
                        &start_tdim, &time);

    for (int ivar=i1_cmp; ivar <= i2_cmp; ivar++)
    {
      float *var = w4d + ivar * siz_icmp;
      int iptr_slice = 0;
      int k = ioslice->slice_z_indx[n];
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni1; i<=ni2; i++) {
          int iptr = i + j * siz_line + k * siz_slice;
          buff[iptr_slice] = var[iptr];
          iptr_slice++;
        }
      }
      nc_put_vara_float(ioslice_nc->ncid_slz[n], 
                        ioslice_nc->varid_slz[n*num_of_vars + ivar],
                        startp, countp, buff);
    }
  }

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

  int num_of_snap = iosnap->num_of_snap;
  size_t siz_line  = gdinfo->siz_iy;
  size_t siz_slice = gdinfo->siz_iz;

  for (int n=0; n<num_of_snap; n++)
  {
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
      size_t startp[] = { iosnap_nc->cur_it[n], 0, 0, 0 };
      size_t countp[] = { 1, snap_nk, snap_nj, snap_ni };
      size_t start_tdim = iosnap_nc->cur_it[n];

      // put time var
      nc_put_var1_float(iosnap_nc->ncid[n],iosnap_nc->timeid[n],&start_tdim,&time);

      // vel
      if (is_run_out_vel == 1 && snap_out_V==1)
      {
        io_snap_pack_buff(w4d + wav->Vx_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+0],
              startp,countp,buff);

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
        io_snap_stress_to_strain_eliso(md->lambda,md->mu,
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
io_snap_nc_create_ac(iosnap_t *iosnap, iosnap_nc_t *iosnap_nc, int *topoid)
{
  int ierr = 0;

  int num_of_snap = iosnap->num_of_snap;
  char **snap_fname = iosnap->fname;

  iosnap_nc->num_of_snap = num_of_snap;
  iosnap_nc->ncid = (int *)malloc(num_of_snap*sizeof(int));
  iosnap_nc->timeid = (int *)malloc(num_of_snap*sizeof(int));

  iosnap_nc->varid_V = (int *)malloc(num_of_snap*CONST_NDIM*sizeof(int));
  iosnap_nc->varid_T = (int *)malloc(num_of_snap*1*sizeof(int));
  iosnap_nc->varid_E = (int *)malloc(num_of_snap*1*sizeof(int));

  // will be used in put step
  iosnap_nc->cur_it = (int *)malloc(num_of_snap*sizeof(int));
  for (int n=0; n<num_of_snap; n++) {
    iosnap_nc->cur_it[n] = 0;
  }

  int *ncid   = iosnap_nc->ncid;
  int *timeid = iosnap_nc->timeid;
  int *varid_V = iosnap_nc->varid_V;
  int *varid_T = iosnap_nc->varid_T;
  int *varid_E = iosnap_nc->varid_E;

  for (int n=0; n<num_of_snap; n++)
  {
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

    if (nc_create(snap_fname[n], NC_CLOBBER, &ncid[n])) M_NCERR;
    if (nc_def_dim(ncid[n], "time", NC_UNLIMITED, &dimid[0])) M_NCERR;
    if (nc_def_dim(ncid[n], "k", snap_nk     , &dimid[1])) M_NCERR;
    if (nc_def_dim(ncid[n], "j", snap_nj     , &dimid[2])) M_NCERR;
    if (nc_def_dim(ncid[n], "i", snap_ni     , &dimid[3])) M_NCERR;
    // time var
    if (nc_def_var(ncid[n], "time", NC_FLOAT, 1, dimid+0, &timeid[n])) M_NCERR;
    // other vars
    if (snap_out_V==1) {
       if (nc_def_var(ncid[n],"Vx",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+0])) M_NCERR;
       if (nc_def_var(ncid[n],"Vy",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+1])) M_NCERR;
       if (nc_def_var(ncid[n],"Vz",NC_FLOAT,4,dimid,&varid_V[n*CONST_NDIM+2])) M_NCERR;
    }
    if (snap_out_T==1) {
       if (nc_def_var(ncid[n],"Tii",NC_FLOAT,4,dimid,&varid_T[n])) M_NCERR;
    }
    if (snap_out_E==1) {
       if (nc_def_var(ncid[n],"Eii",NC_FLOAT,4,dimid,&varid_E[n])) M_NCERR;
    }
    // attribute: index in output snapshot, index w ghost in thread
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

    if (nc_enddef(ncid[n])) M_NCERR;
  } // loop snap

  return ierr;
}

int
io_snap_nc_put_ac(iosnap_t *iosnap,
               iosnap_nc_t *iosnap_nc,
               gd_t    *gdinfo,
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

  int num_of_snap = iosnap->num_of_snap;
  size_t siz_line  = gdinfo->siz_iy;
  size_t siz_slice = gdinfo->siz_iz;

  for (int n=0; n<num_of_snap; n++)
  {
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
      size_t startp[] = { iosnap_nc->cur_it[n], 0, 0, 0 };
      size_t countp[] = { 1, snap_nk, snap_nj, snap_ni };
      size_t start_tdim = iosnap_nc->cur_it[n];

      // put time var
      nc_put_var1_float(iosnap_nc->ncid[n],iosnap_nc->timeid[n],&start_tdim,&time);

      // vel
      if (is_run_out_vel == 1 && snap_out_V==1)
      {
        io_snap_pack_buff(w4d + wav->Vx_pos,
              siz_line,siz_slice,snap_i1,snap_ni,snap_di,snap_j1,snap_nj,
              snap_dj,snap_k1,snap_nk,snap_dk,buff);
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_V[n*CONST_NDIM+0],
              startp,countp,buff);

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
        nc_put_vara_float(iosnap_nc->ncid[n],iosnap_nc->varid_T[n],
              startp,countp,buff);
      }
      if (is_run_out_stress==1 && snap_out_E==1)
      {
        // need to implement
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
  int iptr_snap=0;
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
io_slice_nc_close(ioslice_nc_t *ioslice_nc)
{
  for (int n=0; n < ioslice_nc->num_of_slice_x; n++) {
    nc_close(ioslice_nc->ncid_slx[n]);
  }
  for (int n=0; n < ioslice_nc->num_of_slice_y; n++) {
    nc_close(ioslice_nc->ncid_sly[n]);
  }
  for (int n=0; n < ioslice_nc->num_of_slice_z; n++) {
    nc_close(ioslice_nc->ncid_slz[n]);
  }
  return(0);
}

int
io_snap_nc_close(iosnap_nc_t *iosnap_nc)
{
  for (int n=0; n < iosnap_nc->num_of_snap; n++)
  {
    nc_close(iosnap_nc->ncid[n]);
  }
  return(0);
}

int
io_recv_keep(iorecv_t *iorecv, float *restrict w4d,
             int it, int ncmp, int siz_icmp)
{
  float Lx1, Lx2, Ly1, Ly2, Lz1, Lz2;

  for (int n=0; n < iorecv->total_number; n++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + n;
    int *indx1d = this_recv->indx1d;

    // get coef of linear interp
    Lx2 = this_recv->di; Lx1 = 1.0 - Lx2;
    Ly2 = this_recv->dj; Ly1 = 1.0 - Ly2;
    Lz2 = this_recv->dk; Lz1 = 1.0 - Lz2;
    for (int icmp=0; icmp < ncmp; icmp++)
    {
      int iptr_sta = icmp * iorecv->max_nt + it;
      size_t iptr_cmp = icmp * siz_icmp;
      this_recv->seismo[iptr_sta] = 
          w4d[iptr_cmp + indx1d[0]] * Lx1 * Ly1 * Lz1
        + w4d[iptr_cmp + indx1d[1]] * Lx2 * Ly1 * Lz1
        + w4d[iptr_cmp + indx1d[2]] * Lx1 * Ly2 * Lz1
        + w4d[iptr_cmp + indx1d[3]] * Lx2 * Ly2 * Lz1
        + w4d[iptr_cmp + indx1d[5]] * Lx1 * Ly1 * Lz2
        + w4d[iptr_cmp + indx1d[5]] * Lx2 * Ly1 * Lz2
        + w4d[iptr_cmp + indx1d[6]] * Lx1 * Ly2 * Lz2
        + w4d[iptr_cmp + indx1d[7]] * Lx2 * Ly2 * Lz2;
    }

    // need to implement interp, flollowing just take value
    //for (int icmp=0; icmp < ncmp; icmp++) {
    //  int iptr_sta = icmp * iorecv->max_nt + it;
    //  this_recv->seismo[iptr_sta] = w4d[icmp*siz_icmp + iptr];
    //}
    
  }

  return(0);
}

int
io_line_keep(ioline_t *ioline, float *restrict w4d,
             int it, int ncmp, int siz_icmp)
{
  for (int n=0; n < ioline->num_of_lines; n++)
  {
    int   *this_line_iptr   = ioline->recv_iptr[n];
    float *this_line_seismo = ioline->recv_seismo[n];
  
    for (int ir=0; ir < ioline->line_nr[n]; ir++)
    {
      int iptr = this_line_iptr[ir];
      float *this_seismo = this_line_seismo + ir * ioline->max_nt * ncmp;
      for (int icmp=0; icmp < ncmp; icmp++)
      {
        int iptr_seismo = icmp * ioline->max_nt + it;
        this_seismo[iptr_seismo] = w4d[icmp*siz_icmp + iptr];
      }
    }
  }
}

int
io_recv_output_sac(iorecv_t *iorecv,
                   float dt,
                   int num_of_vars,
                   int visco_type,
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

  for (int ir=0; ir < iorecv->total_number; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;
    
    if (visco_type == CONST_VISCO_GMB)
        num_of_vars = 9; // not output J for attenuation

    //fprintf(stdout,"=== Debug: num_of_vars=%d\n",num_of_vars);fflush(stdout);
    for (int icmp=0; icmp < num_of_vars; icmp++)
    {
      //fprintf(stdout,"=== Debug: icmp=%d\n",icmp);fflush(stdout);

      float *this_trace = this_recv->seismo + icmp * iorecv->max_nt;

      sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm,
                      this_recv->name, cmp_name[icmp]);

      //fprintf(stdout,"=== Debug: icmp=%d,ou_file=%s\n",icmp,ou_file);fflush(stdout);

      sacExport1C1R(ou_file,
            this_trace,
            evt_x, evt_y, evt_z, evt_d,
            this_recv->x, this_recv->y, this_recv->z,
            dt, dt, iorecv->max_nt, err_message);
    }
  }

  return 0;
}

// calculate and output strain cmp for elastic medium
//   do not find a better file to hold this func
//   temporarily put here

int
io_recv_output_sac_el_iso_strain(iorecv_t *iorecv,
                     float *restrict lam3d,
                     float *restrict mu3d,
                     float dt,
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

  for (int ir=0; ir < iorecv->total_number; ir++)
  {
    iorecv_one_t *this_recv = iorecv->recvone + ir;
    int iptr = this_recv->indx1d[0];

    float lam = lam3d[iptr];
    float mu  =  mu3d[iptr];

    // cmp seq hard-coded, need to revise in the future
    float *Txx = this_recv->seismo + 3 * iorecv->max_nt;
    float *Tyy = this_recv->seismo + 4 * iorecv->max_nt;
    float *Tzz = this_recv->seismo + 5 * iorecv->max_nt;
    float *Tyz = this_recv->seismo + 6 * iorecv->max_nt;
    float *Txz = this_recv->seismo + 7 * iorecv->max_nt;
    float *Txy = this_recv->seismo + 8 * iorecv->max_nt;

    float E1 = (lam + mu) / (mu * ( 3.0 * lam + 2.0 * mu));
    float E2 = - lam / ( 2.0 * mu * (3.0 * lam + 2.0 * mu));
    float E3 = 1.0 / mu;

    // conver to strain per time step
    for (int it = 0; it < iorecv->max_nt; it++)
    {
      float E0 = E2 * (Txx[it] + Tyy[it] + Tzz[it]);

      Txx[it] = E0 - (E2 - E1) * Txx[it];
      Tyy[it] = E0 - (E2 - E1) * Tyy[it];
      Tzz[it] = E0 - (E2 - E1) * Tzz[it];
      Tyz[it] = 0.5 * E3 * Tyz[it];
      Txz[it] = 0.5 * E3 * Txz[it];
      Txy[it] = 0.5 * E3 * Txy[it];
    }

    // output to sca file
    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exx");
    sacExport1C1R(ou_file,Txx,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, this_recv->y, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Eyy");
    sacExport1C1R(ou_file,Tyy,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, this_recv->y, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Ezz");
    sacExport1C1R(ou_file,Tzz,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, this_recv->y, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Eyz");
    sacExport1C1R(ou_file,Tyz,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, this_recv->y, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exz");
    sacExport1C1R(ou_file,Txz,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, this_recv->y, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);

    sprintf(ou_file,"%s/%s.%s.%s.sac", output_dir, evtnm, this_recv->name, "Exy");
    sacExport1C1R(ou_file,Txy,evt_x, evt_y, evt_z, evt_d,
          this_recv->x, this_recv->y, this_recv->z,
          dt, dt, iorecv->max_nt, err_message);
  } // loop ir

  return 0;
}

int
io_recv_output_sac_el_vti_strain(iorecv_t *iorecv,
                        float *restrict c11, float *restrict c13,
                        float *restrict c33, float *restrict c55,
                        float *restrict c66,
                        float dt,
                        char *evtnm,
                        char *output_dir,
                        char *err_message)
{
  //not implement



  return 0;
}
int
io_recv_output_sac_el_aniso_strain(iorecv_t *iorecv,
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
                        float *restrict c66d,
                        float dt,
                        char *evtnm,
                        char *output_dir,
                        char *err_message)
{
  //not implement



  return 0;
}

int
io_line_output_sac(ioline_t *ioline,
      float dt, char **cmp_name, char *evtnm, char *output_dir)
{
  // use fake evt_x etc. since did not implement gather evt_x by mpi
  float evt_x = 0.0;
  float evt_y = 0.0;
  float evt_z = 0.0;
  float evt_d = 0.0;
  char ou_file[CONST_MAX_STRLEN];
  char err_message[CONST_MAX_STRLEN];
  
  for (int n=0; n < ioline->num_of_lines; n++)
  {
    int   *this_line_iptr   = ioline->recv_iptr[n];
    float *this_line_seismo = ioline->recv_seismo[n];

    for (int ir=0; ir < ioline->line_nr[n]; ir++)
    {
      float *this_seismo = this_line_seismo + ir * ioline->max_nt * ioline->ncmp;

      for (int icmp=0; icmp < ioline->ncmp; icmp++)
      {
        float *this_trace = this_seismo + icmp * ioline->max_nt;

        sprintf(ou_file,"%s/%s.%s.no%d.%s.sac", output_dir,evtnm,
                  ioline->line_name[n],ioline->recv_seq[n][ir],
                  cmp_name[icmp]);

        sacExport1C1R(ou_file,
              this_trace,
              evt_x, evt_y, evt_z, evt_d,
              ioline->recv_x[n][ir],
              ioline->recv_y[n][ir],
              ioline->recv_z[n][ir],
              dt, dt, ioline->max_nt, err_message);
      } // icmp
    } // ir
  } // line
}

int
ioslice_print(ioslice_t *ioslice)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> slice output information:\n");
  fprintf(stdout, "-------------------------------------------------------\n");

  fprintf(stdout, "--> num_of_slice_x = %d\n", ioslice->num_of_slice_x);
  for (int n=0; n<ioslice->num_of_slice_x; n++)
  {
    fprintf(stdout, "  #%d, i=%d, fname=%s\n", n, ioslice->slice_x_indx[n],ioslice->slice_x_fname[n]);
  }
  fprintf(stdout, "--> num_of_slice_y = %d\n", ioslice->num_of_slice_y);
  for (int n=0; n<ioslice->num_of_slice_y; n++)
  {
    fprintf(stdout, "  #%d, j=%d, fname=%s\n", n, ioslice->slice_y_indx[n],ioslice->slice_y_fname[n]);
  }
  fprintf(stdout, "--> num_of_slice_z = %d\n", ioslice->num_of_slice_z);
  for (int n=0; n<ioslice->num_of_slice_z; n++)
  {
    fprintf(stdout, "  #%d, k=%d, fname=%s\n", n, ioslice->slice_z_indx[n],ioslice->slice_z_fname[n]);
  }

  return(0);
}

int
iosnap_print(iosnap_t *iosnap)
{    
  fprintf(stdout, "--> num_of_snap = %d\n", iosnap->num_of_snap);
  fprintf(stdout, "#   i0 j0 k0 ni nj nk di dj dk it0 dit vel stress strain gi1 gj1 gk1\n");
  for (int n=0; n < iosnap->num_of_snap; n++)
  {
    fprintf(stdout, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
              n,
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

int
PG_slice_output(float *PG, gd_t *gdinfo, char *output_dir, char *frame_coords, int *topoid)
{
  // output one time z slice
  // used for PGV PGA and PGD
  // cmp is PGV PGA PGD, component x, y, z
  int nx = gdinfo->nx; 
  int ny = gdinfo->ny;
  int ni = gdinfo->ni; 
  int nj = gdinfo->nj;
  int gni1 = gdinfo->gni1; 
  int gnj1 = gdinfo->gnj1; 
  char PG_cmp[CONST_NDIM_5][CONST_MAX_STRLEN] = {"PGV","PGVh","PGVx","PGVy","PGVz",
                                                 "PGA","PGAh","PGAx","PGAy","PGAz",
                                                 "PGD","PGDh","PGDx","PGDy","PGDz"}; 
  char out_file[CONST_MAX_STRLEN];
  sprintf(out_file,"%s/%s_%s.nc",output_dir,"PG_V_A_D",frame_coords);

  // create PGV output file
  int dimid[2];
  int varid[CONST_NDIM_5], ncid;
  int i,ierr;
  ierr = nc_create(out_file, NC_CLOBBER, &ncid);
  if (ierr != NC_NOERR){
      fprintf(stderr,"creat PGV nc error: %s\n", nc_strerror(ierr));
      exit(-1);
     }

  if(nc_def_dim(ncid, "j", ny, &dimid[0])) M_NCERR;
  if(nc_def_dim(ncid, "i", nx, &dimid[1])) M_NCERR;

  // define vars
  for(int i=0; i<CONST_NDIM_5; i++)
  {
    if(nc_def_var(ncid, PG_cmp[i], NC_FLOAT,2,dimid, &varid[i])) M_NCERR;
  }
  int g_start[2] = {gni1,gnj1}; 
  int phy_size[2] = {ni,nj}; 
  if(nc_put_att_int(ncid,NC_GLOBAL,"global_index_of_first_physical_points",
                    NC_INT,2,g_start)) M_NCERR;
  if(nc_put_att_int(ncid,NC_GLOBAL,"count_index_of_physical_points",
                    NC_INT,2,phy_size)) M_NCERR;
  if(nc_put_att_int(ncid,NC_GLOBAL,"coords_of_mpi_topo",
                    NC_INT,2,topoid)) M_NCERR;

  if(nc_enddef(ncid)) M_NCERR;

  // add vars
  for(int i=0; i<CONST_NDIM_5; i++)
  {
  float *ptr = PG + i*nx*ny; 
  ierr = nc_put_var_float(ncid,varid[i],ptr);
  }
  // close file
  ierr = nc_close(ncid);
  if(ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    exit(-1);
    }
  return 0;
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
