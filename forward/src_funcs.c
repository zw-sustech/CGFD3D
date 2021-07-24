/*
 * source term related processing
 */

// todo:

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "fdlib_math.h"
#include "interp.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "src_funcs.h"
#include "isPointInHexahedron.h"
 
/*
 * for single force or moment source term, with Gaussian spatial smoothing
 */
void
src_gen_single_point_gauss(size_t siz_line,
                           size_t siz_slice,
                           float t0,
                           float dt,
                           int   num_of_stages,
                           float *rk_stage_time,
                           int   glob_phys_ix1, // gloabl start index along x this thread
                           int   glob_phys_ix2, // gloabl end index along x
                           int   glob_phys_iy1,
                           int   glob_phys_iy2,
                           int   glob_phys_iz1,
                           int   glob_phys_iz2,
                           int   ni1,
                           int   ni2,
                           int   nj1,
                           int   nj2,
                           int   nk1,
                           int   nk2,
                           int   npoint_half_ext,
                           int   npoint_ghosts,
                           float *x3d,
                           float *y3d,
                           float *z3d,
                           int   *source_gridindex,
                           float *source_coords,
                           float *force_vector,
                           float *moment_tensor,
                           char  *wavelet_name,
                           float *wavelet_coefs,
                           float wavelet_tstart,
                           float wavelet_tend,
                           MPI_Comm comm, 
                           int myid,
                           // following output
                           int  *num_of_force, // inout: if force source, if in this thread
                           int **restrict p_force_info,
                           float  **restrict p_force_vec_stf,
                           int    **restrict p_force_ext_indx,
                           float  **restrict p_force_ext_coef,
                           int  *num_of_moment, // inout: if moment source, if in this thread
                           int    **restrict p_moment_info,
                           float  **restrict p_moment_ten_rate,
                           int    **restrict p_moment_ext_indx,
                           float  **restrict p_moment_ext_coef,
                           int verbose)
{
  // get total elem of exted src region for a single point
  //int siz_ext = 7 * 7 * 7;
  int siz_ext = (2*npoint_half_ext+1)*(2*npoint_half_ext+1)*(2*npoint_half_ext+1);

  // get initial number
  int nforce  = *num_of_force;
  int nmoment = *num_of_moment;

  // use local var to represent index
  int sgpi = source_gridindex[0]; // g: global, p: phys
  int sgpj = source_gridindex[1];
  int sgpk = source_gridindex[2];

  // convert time to index
  int  it_begin = (int) (wavelet_tstart / dt);
  int  it_end   = (int) ((wavelet_tend / dt + 0.5));
  int  nt_total_wavelet = it_end - it_begin + 1;
  float *wavelet_values = NULL;

  // default local index and relative shift to -1 and 0
  int si=-1; int sj=-1; int sk=-1;
  float sx_inc = 0.0; float sy_inc = 0.0; float sz_inc = 0.0;

  // workspace 3d var for distance calculation
  float *wrk3d=NULL;

  // locate into this thread
  if (sgpi < 0 || sgpj < 0 || sgpk < 0) {
    // use coord to locate local index
    fprintf(stdout,"locate source by coord (%f,%f,%f) ...\n",
        source_coords[0],source_coords[1],source_coords[2]); 
    //fprintf(stdout,"   not implemented yet\n"); 
    fflush(stdout);

    int nx = (ni2-ni1+1)+2*npoint_ghosts;
    int ny = (nj2-nj1+1)+2*npoint_ghosts;
    int nz = (nk2-nk1+1)+2*npoint_ghosts;
    wrk3d = (float *) fdlib_mem_calloc_1d_float(nx*ny*nz,0.0,
        "src_gen_single_point_gauss");

    // if located in this thread
    int is_here = src_coord2index(source_coords[0],source_coords[1],source_coords[2],
                                  nx, ny, nz, 
                                  ni1,ni2,nj1,nj2,nk1,nk2,
                                  x3d, y3d, z3d, wrk3d,
                                  &si, &sj, &sk,
                                  &sx_inc, &sy_inc, &sz_inc);
    if ( is_here == 1)
    {
      // conver to global index
      sgpi = si - npoint_ghosts + glob_phys_ix1;
      sgpj = sj - npoint_ghosts + glob_phys_iy1;
      sgpk = sk - npoint_ghosts + glob_phys_iz1;
      fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
      fprintf(stdout," -- located to global index = %d %d %d\n", sgpi,sgpj,sgpk);
      fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
    } else {
      fprintf(stdout," -- not in this thread %d\n", myid);
    }

    // free allocated vars
    free(wrk3d);
  }

  // reduce global index and shift values
  int sendbufi = sgpi;
  MPI_Allreduce(&sendbufi, &sgpi, 1, MPI_INT, MPI_MAX, comm);
  sendbufi = sgpj;
  MPI_Allreduce(&sendbufi, &sgpj, 1, MPI_INT, MPI_MAX, comm);
  sendbufi = sgpk;
  MPI_Allreduce(&sendbufi, &sgpk, 1, MPI_INT, MPI_MAX, comm);

  float sendbuf = sx_inc;
  MPI_Allreduce(&sendbuf, &sx_inc, 1, MPI_INT, MPI_SUM, comm);
  sendbuf = sy_inc;
  MPI_Allreduce(&sendbuf, &sy_inc, 1, MPI_INT, MPI_SUM, comm);
  sendbuf = sz_inc;
  MPI_Allreduce(&sendbuf, &sz_inc, 1, MPI_INT, MPI_SUM, comm);

  fprintf(stdout," --myid=%d,index=%d %d %d,shift = %f %f %f\n",
      myid,sgpi,sgpj,sgpk, sx_inc,sy_inc,sz_inc);

  // use grid index to check if in this thread after extend
  if (sgpi-npoint_half_ext <= glob_phys_ix2 && // exted left point is less than right bdry
      sgpi+npoint_half_ext >= glob_phys_ix1 && // exted right point is larger than left bdry
      sgpj-npoint_half_ext <= glob_phys_iy2 && 
      sgpj+npoint_half_ext >= glob_phys_iy1 &&
      sgpk-npoint_half_ext <= glob_phys_iz2 && 
      sgpk+npoint_half_ext >= glob_phys_iz1)
  {
    // at least one extend point in this thread
    // convert to local index
    si = sgpi - glob_phys_ix1 + npoint_ghosts;
    sj = sgpj - glob_phys_iy1 + npoint_ghosts;
    sk = sgpk - glob_phys_iz1 + npoint_ghosts;
  }
  else
  {  // no source in this thread
    nforce  = 0;
    nmoment = 0;
  }
  // if source here, cal wavelet series
  if (nt_total_wavelet > 0) {
    wavelet_values = (float *)fdlib_mem_calloc_1d_float(
        nt_total_wavelet*num_of_stages,0.0,"src_gen_single_point_gauss");
  }

  if (nforce > 0 || nmoment > 0)
  {
    for (int it=it_begin; it<=it_end; it++)
    {
      int it_to_itbegin = (it - it_begin);
      for (int istage=0; istage<num_of_stages; istage++)
      {
        float t = it * dt + t0 + rk_stage_time[istage] * dt;
        float stf_val;
        if (strcmp(wavelet_name, "ricker")==0) {
          stf_val = fun_ricker(t, wavelet_coefs[0], wavelet_coefs[1]);
        } else if (strcmp(wavelet_name, "gaussian")==0) {
          stf_val = fun_gauss(t, wavelet_coefs[0], wavelet_coefs[1]);
        } else if (strcmp(wavelet_name, "ricker_deriv")==0) {
          stf_val = fun_ricker_deriv(t, wavelet_coefs[0], wavelet_coefs[1]);
        } else if (strcmp(wavelet_name, "gaussian_deriv")==0) {
          stf_val = fun_gauss_deriv(t, wavelet_coefs[0], wavelet_coefs[1]);
        } else{
          fprintf(stderr,"wavelet_name=%s\n", wavelet_name); 
          fprintf(stderr,"   not implemented yet\n"); 
          fflush(stderr);
        }
        // save to vector
        wavelet_values[it_to_itbegin*num_of_stages+istage] = stf_val;
      }
    }
  }
  // set moment
  int *moment_info = NULL;
  if (nmoment > 0) // source type is moment
  {
    // allocate info and return to main function
    moment_info = (int *)fdlib_mem_calloc_1d_int(nmoment*M_SRC_INFO_NVAL, 1,
        "src_gen_single_point_gauss");
    // alloc ten_rate
    int nt_moment = nt_total_wavelet;
    float *moment_ten_rate = (float *)fdlib_mem_calloc_1d_float(
        nmoment*nt_moment*6*num_of_stages,0.0,"src_gen_single_point_gauss");
    // ext
    float *moment_ext_coef = (float *)malloc(nmoment*siz_ext * sizeof(float));
    int   *moment_ext_indx = (int   *)malloc(nmoment*siz_ext * sizeof(int  ));

    // save values to inner var
    moment_info[M_SRC_INFO_SEQ_SI   ] = si;
    moment_info[M_SRC_INFO_SEQ_SJ   ] = sj;
    moment_info[M_SRC_INFO_SEQ_SK   ] = sk;
    moment_info[M_SRC_INFO_SEQ_POS  ] = 0;
    moment_info[M_SRC_INFO_SEQ_ITBEG] = it_begin;
    moment_info[M_SRC_INFO_SEQ_ITEND] = it_end;

    // for this point
    int ipos = moment_info[M_SRC_INFO_SEQ_POS];
    int it1  = moment_info[M_SRC_INFO_SEQ_ITBEG];
    int it2  = moment_info[M_SRC_INFO_SEQ_ITEND];
    float *this_ten_rate = moment_ten_rate + ipos;

    for (int icmp=0; icmp<6; icmp++)
    {
      for (int it=it1; it<=it2; it++)
      {
        int it_to_it1 = (it - it1);
        for (int istage=0; istage<num_of_stages; istage++)
        {
          int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_moment,num_of_stages);
          float stf_val = wavelet_values[it_to_it1*num_of_stages+istage];
          this_ten_rate[iptr] = stf_val * moment_tensor[icmp];
        }
      }
    }

    cal_norm_delt3d(moment_ext_coef, sx_inc, sy_inc, sz_inc, 1.5, 1.5, 1.5, 3);

    size_t iptr_s = 0;
    for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
    {
      if (k<nk1 || k>nk2) continue;

      for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
      {
        if (j<nj1 || j>nj2) continue;

        for (int i=si-npoint_half_ext; i<=si+npoint_half_ext; i++)
        {
          if (i<ni1 || i>ni2) continue;

          // Note index need match coef
          int iptr = i + j * siz_line + k * siz_slice;
          int iptr1 = (i-(si-npoint_half_ext)) + 7 * (j-(sj-npoint_half_ext)) + 7 * 7 *(k-(sk-npoint_half_ext));
          moment_ext_indx[iptr_s] = iptr;
          moment_ext_coef[iptr_s] = moment_ext_coef[iptr1];
          iptr_s++;
        }
      }
    }
    // only count index inside phys region for this thread
    moment_info[M_SRC_INFO_SEQ_NEXT_MAX ] = siz_ext;
    moment_info[M_SRC_INFO_SEQ_NEXT_THIS] = iptr_s;

    *p_moment_ten_rate = moment_ten_rate;
    *p_moment_ext_indx = moment_ext_indx;
    *p_moment_ext_coef = moment_ext_coef;
  }
  // set force
  int *force_info = NULL;
  if (nforce > 0) // source type is force
  {
    // allocate info and return to main function
    force_info = (int *)fdlib_mem_calloc_1d_int(nforce*M_SRC_INFO_NVAL, 1,
        "src_gen_single_point_gauss");
    // stf
    int nt_force = nt_total_wavelet;
    float *force_vec_stf = (float *)fdlib_mem_calloc_1d_float(
        nforce*nt_force*FD_NDIM*num_of_stages,0.0,"src_gen_single_point_gauss");
    // ext
    float *force_ext_coef = (float *)malloc(nforce*siz_ext * sizeof(float));
    int   *force_ext_indx = (int   *)malloc(nforce*siz_ext * sizeof(int  ));

    // save values to inner var
    force_info[M_SRC_INFO_SEQ_SI   ] = si;
    force_info[M_SRC_INFO_SEQ_SJ   ] = sj;
    force_info[M_SRC_INFO_SEQ_SK   ] = sk;
    force_info[M_SRC_INFO_SEQ_POS  ] = 0;
    force_info[M_SRC_INFO_SEQ_ITBEG] = it_begin;
    force_info[M_SRC_INFO_SEQ_ITEND] = it_end;

    int ipos = force_info[M_SRC_INFO_SEQ_POS];
    int it1  = force_info[M_SRC_INFO_SEQ_ITBEG];
    int it2  = force_info[M_SRC_INFO_SEQ_ITEND];
    float *this_vec_stf = force_vec_stf + ipos;

    for (int icmp=0; icmp<FD_NDIM; icmp++)
    {
      for (int it=it1; it<=it2; it++)
      {
        int it_to_it1 = (it - it1);
        for (int istage=0; istage<num_of_stages; istage++)
        {
          int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_force,num_of_stages);
          float stf_val = wavelet_values[it_to_it1*num_of_stages+istage];
          this_vec_stf[iptr] = stf_val * force_vector[icmp];
        }
      }
    }
    cal_norm_delt3d(force_ext_coef, sx_inc, sy_inc, sz_inc, 1.5, 1.5, 1.5, 3);

    size_t iptr_s = 0;
    for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
    {
      if (k<nk1 || k>nk2) continue;

      for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
      {
        if (j<nj1 || j>nj2) continue;

        for (int i=si-npoint_half_ext; i<=si+npoint_half_ext; i++)
        {
          if (i<ni1 || i>ni2) continue;

          // Note index need match coef
          int iptr = i + j * siz_line + k * siz_slice;
          int iptr1 = (i-(si-npoint_half_ext)) + 7 * (j-(sj-npoint_half_ext)) + 7 * 7 *(k-(sk-npoint_half_ext));
          force_ext_indx[iptr_s] = iptr;
          force_ext_coef[iptr_s] = force_ext_coef[iptr1];
          iptr_s++;
        }
      }
    }

    force_info[M_SRC_INFO_SEQ_NEXT_MAX ] = siz_ext;
    force_info[M_SRC_INFO_SEQ_NEXT_THIS] = iptr_s;

    *p_force_vec_stf = force_vec_stf;
    *p_force_ext_indx = force_ext_indx;
    *p_force_ext_coef = force_ext_coef;
  }

  *num_of_force = nforce;
  *p_force_info = force_info;
  *num_of_moment = nmoment;
  *p_moment_info = moment_info;

  if (wavelet_values!=NULL) free(wavelet_values);

}


int
src_read_locate_valsrc(char *pfilepath,
                       size_t siz_line,
                       size_t siz_slice,
                       float t0,
                       float dt,
                       int   num_of_stages,
                       float *rk_stage_time,
                       int   glob_phys_ix1, // gloabl start index along x this thread
                       int   glob_phys_ix2, // gloabl end index along x
                       int   glob_phys_iy1,
                       int   glob_phys_iy2,
                       int   glob_phys_iz1,
                       int   glob_phys_iz2,
                       int   ni1,
                       int   ni2,
                       int   nj1,
                       int   nj2,
                       int   nk1,
                       int   nk2,
                       int   npoint_half_ext,
                       int   npoint_ghosts,
                       float *x3d,
                       float *y3d,
                       float *z3d,
                       MPI_Comm comm,
                       int myid,
                       // following output
                       char **p_event_name,
                       int  *num_of_force, // inout: if force source, if in this thread
                       int **restrict p_force_info,
                       float  **restrict p_force_vec_stf,
                       int    **restrict p_force_ext_indx,
                       float  **restrict p_force_ext_coef,
                       int  *num_of_moment, // inout: if moment source, if in this thread
                       int    **restrict p_moment_info,
                       float  **restrict p_moment_ten_rate,
                       int    **restrict p_moment_ext_indx,
                       float  **restrict p_moment_ext_coef,
                       int verbose)
{
  FILE *fp =NULL;
  char str[500];
  char *event_name = (char *)malloc(500*sizeof(char));
  int num_force;
  int num_moment;
  int nt_in;       // numbers_time_steps from inputfile
  float dt_in;     // time_step from inputfile
  // read sample source value from inputfile
  if ((fp = fopen(pfilepath, "r"))==NULL) fprintf(stdout,"fail to open");

  fgets(str,500,fp);
  sscanf(str,"%s",event_name);
  fgets(str,500,fp);
  sscanf(str,"%d %d",&num_force,&num_moment);
  fgets(str,500,fp);
  sscanf(str,"%f %d",&dt_in,&nt_in);

  *p_event_name = event_name;
  float *force_coords = NULL;
  int   *force_sur = NULL;
  float *moment_coords = NULL;
  int   *moment_sur = NULL;

  if(num_force>0)
  {
    force_coords = (float *) malloc(3 * num_force * sizeof(float));
    force_sur = (int *) malloc(num_force * sizeof(int)); // relocated_to_surface
  }
  if(num_moment>0)
  {
    moment_coords = (float *) malloc(3 * num_moment * sizeof(float));
    moment_sur = (int *) malloc(num_moment * sizeof(int));                                         
  }

  // first read coords to determine src whether in this thread
  for (int i=0; i<num_force; i++)
  {
    fgets(str,500,fp);
    sscanf(str,"%f %f %f %d",&force_coords[3*i+0],&force_coords[3*i+1],&force_coords[3*i+2],&force_sur[i]);
  }
  for (int i=0; i<num_moment; i++)
  {
    fgets(str,500,fp);
    sscanf(str,"%f %f %f %d",&moment_coords[3*i+0],&moment_coords[3*i+1],&moment_coords[3*i+2],&moment_sur[i]);
  }

  int nx = (ni2-ni1+1)+2*npoint_ghosts;
  int ny = (nj2-nj1+1)+2*npoint_ghosts;
  int nz = (nk2-nk1+1)+2*npoint_ghosts;  
  // workspace 3d var for distance calculation
  float *wrk3d = (float *) fdlib_mem_calloc_1d_float(nx*ny*nz,0.0, "src_read_locat_valsrc");

  int *force_global_index = NULL;
  int nforce = 0;
  int *index_force = NULL;
  int *force_local_index = NULL;
  float *force_index_inc = NULL;
  if(num_force>0)
  {
    // judge force in this thread 
    force_global_index = (int *)malloc(3 * num_force * sizeof(int));
    force_index_inc = (float *)malloc(3 * num_force*sizeof(float));
    for (int i=0; i<num_force; i++)
    {
      fprintf(stdout,"locate force by coord  ...\n"); 
      // default global index and relative shift to -1 and 0
      int si, sj, sk;
      float sx_inc = 0.0; float sy_inc = 0.0; float sz_inc = 0.0;
      int sgpi=-1; int sgpj=-1; int sgpk=-1;
      // if located in this thread
      int is_here = src_coord2index(force_coords[3*i+0],force_coords[3*i+1],force_coords[3*i+2],
                                    nx, ny, nz, 
                                    ni1,ni2,nj1,nj2,nk1,nk2,
                                    x3d, y3d, z3d, wrk3d,
                                    &si, &sj, &sk,
                                    &sx_inc, &sy_inc, &sz_inc);
      if ( is_here == 1)
      {
        // conver to global index
        sgpi = si - npoint_ghosts + glob_phys_ix1;
        sgpj = sj - npoint_ghosts + glob_phys_iy1;
        sgpk = sk - npoint_ghosts + glob_phys_iz1;
        fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
        fprintf(stdout," -- located to global index = %d %d %d\n", sgpi,sgpj,sgpk);
        fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
      } else {
        fprintf(stdout," this force not in this thread %d\n", myid);
      }
      force_index_inc[3*i+0] = sx_inc;
      force_index_inc[3*i+1] = sy_inc;
      force_index_inc[3*i+2] = sz_inc;
      force_global_index[3*i+0] = sgpi;
      force_global_index[3*i+1] = sgpj;
      force_global_index[3*i+2] = sgpk;
    }

    // reduce all source points global index and shift values
    // Note select MPI_MAX for global index and MPI_SUM for shift value.
    // 0 -> x, 1->y, 2->z.
    int *reduce_force_global_index = (int*)malloc(3*num_force*sizeof(int));
    MPI_Allreduce(force_global_index, reduce_force_global_index, 3*num_force, MPI_INT, MPI_MAX, comm);

    float *reduce_force_index_inc = (float*)malloc(3*num_force*sizeof(float));
    MPI_Allreduce(force_index_inc,reduce_force_index_inc, 3*num_force, MPI_INT, MPI_SUM, comm);

    for (int i=0; i<num_force; i++)
    {
      fprintf(stdout,"force %d --myid=%d,index=%d %d %d,shift = %f %f %f\n",
              i, myid,reduce_force_global_index[3*i+0],reduce_force_global_index[3*i+1],reduce_force_global_index[3*i+2],
              reduce_force_index_inc[3*i+0],reduce_force_index_inc[3*i+1],reduce_force_index_inc[3*i+2]);
    }

    index_force = (int *)malloc(num_force*sizeof(int));
    force_local_index = (int*)malloc(3*num_force*sizeof(int));
    for (int i=0; i<num_force; i++)
    {
      // use grid index to check if in this thread after extend
      if( reduce_force_global_index[3*i+0]-npoint_half_ext <= glob_phys_ix2 && // exted left point is less than right bdry
          reduce_force_global_index[3*i+0]+npoint_half_ext >= glob_phys_ix1 && // exted right point is larger than left bdry
          reduce_force_global_index[3*i+1]-npoint_half_ext <= glob_phys_iy2 && 
          reduce_force_global_index[3*i+1]+npoint_half_ext >= glob_phys_iy1 &&
          reduce_force_global_index[3*i+2]-npoint_half_ext <= glob_phys_iz2 && 
          reduce_force_global_index[3*i+2]+npoint_half_ext >= glob_phys_iz1)
      {
        // at least one extend point in this thread
        // convert to local index
        force_local_index[3*i+0] = reduce_force_global_index[3*i+0] - glob_phys_ix1 + npoint_ghosts;
        force_local_index[3*i+1] = reduce_force_global_index[3*i+1] - glob_phys_iy1 + npoint_ghosts;
        force_local_index[3*i+2] = reduce_force_global_index[3*i+2] - glob_phys_iz1 + npoint_ghosts;
        index_force[nforce] = i;
        nforce++;
        //fprintf(stdout,"myid is %d,si,sj,sk is %d,%d,%d\n",myid,force_local_index[3*i+0],force_local_index[3*i+1],force_local_index[3*i+2]);
      }
      else
      {
        //default local index is -1
        force_local_index[3*i+0] = -1;
        force_local_index[3*i+1] = -1; 
        force_local_index[3*i+2] = -1;
      }
    }
  }

  int *moment_global_index = NULL;
  int nmoment = 0;
  int *index_moment = NULL;
  int *moment_local_index = NULL;
  float *moment_index_inc = NULL;
  if(num_moment>0)
  {
    // judge moment in this thread
    moment_global_index = (int *)malloc(3 * num_moment*sizeof(int));
    moment_index_inc = (float *)malloc(3 * num_moment*sizeof(float));
    for (int i =0; i<num_moment; i++)
    {
      fprintf(stdout,"locate moment by coord  ...\n"); 
      // default global index and relative shift to -1 and 0
      int si, sj, sk;
      float sx_inc = 0.0; float sy_inc = 0.0; float sz_inc = 0.0;
      int sgpi=-1; int sgpj=-1; int sgpk=-1;
      // if located in this thread
      int is_here = src_coord2index(moment_coords[3*i+0],moment_coords[3*i+1],moment_coords[3*i+2],
                                    nx, ny, nz, 
                                    ni1,ni2,nj1,nj2,nk1,nk2,
                                    x3d, y3d, z3d, wrk3d,
                                    &si, &sj, &sk,
                                    &sx_inc, &sy_inc, &sz_inc);
      if ( is_here == 1)
      {
        // conver to global index
        sgpi = si - npoint_ghosts + glob_phys_ix1;
        sgpj = sj - npoint_ghosts + glob_phys_iy1;
        sgpk = sk - npoint_ghosts + glob_phys_iz1;
        fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
        fprintf(stdout," -- located to global index = %d %d %d\n", sgpi,sgpj,sgpk);
        fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
      } else {
        fprintf(stdout," this moment not in this thread %d\n", myid);
      }
      moment_index_inc[3*i+0] = sx_inc;
      moment_index_inc[3*i+1] = sy_inc;
      moment_index_inc[3*i+2] = sz_inc;
      moment_global_index[3*i+0] = sgpi;
      moment_global_index[3*i+1] = sgpj;
      moment_global_index[3*i+2] = sgpk;
    }

    // reduce all source points global index and shift values
    // Note select MPI_MAX for global index and MPI_SUM for shift value.
    // 0 -> x, 1->y, 2->z.
    int *reduce_moment_global_index = (int*)malloc(3*num_moment*sizeof(int));
    MPI_Allreduce(moment_global_index, reduce_moment_global_index, 3*num_moment, MPI_INT, MPI_MAX, comm);

    float *reduce_moment_index_inc = (float*)malloc(3*num_moment*sizeof(float));
    MPI_Allreduce(moment_index_inc,reduce_moment_index_inc, 3*num_moment, MPI_INT, MPI_SUM, comm);

    for (int i=0; i<num_moment; i++)
    {
      fprintf(stdout,"moment %d --myid=%d,index=%d %d %d,shift = %f %f %f\n",
          i, myid,reduce_moment_global_index[3*i+0],reduce_moment_global_index[3*i+1],reduce_moment_global_index[3*i+2],
          reduce_moment_index_inc[3*i+0],reduce_moment_index_inc[3*i+1],reduce_moment_index_inc[3*i+2]);
    }

    index_moment = (int *)malloc(num_moment*sizeof(int));
    moment_local_index = (int *)malloc(3 * num_moment*sizeof(int));
    // use grid index to check if in this thread after extend
    for(int i=0; i<num_moment; i++)
    {
      if (reduce_moment_global_index[3*i+0]-npoint_half_ext <= glob_phys_ix2 && // exted left point is less than right bdry
          reduce_moment_global_index[3*i+0]+npoint_half_ext >= glob_phys_ix1 && // exted right point is larger than left bdry
          reduce_moment_global_index[3*i+1]-npoint_half_ext <= glob_phys_iy2 && 
          reduce_moment_global_index[3*i+1]+npoint_half_ext >= glob_phys_iy1 &&
          reduce_moment_global_index[3*i+2]-npoint_half_ext <= glob_phys_iz2 && 
          reduce_moment_global_index[3*i+2]+npoint_half_ext >= glob_phys_iz1)
      {
        // at least one extend point in this thread
        // convert to local index
        moment_local_index[3*i+0]= reduce_moment_global_index[3*i+0] - glob_phys_ix1 + npoint_ghosts;
        moment_local_index[3*i+1]= reduce_moment_global_index[3*i+1] - glob_phys_iy1 + npoint_ghosts;
        moment_local_index[3*i+2]= reduce_moment_global_index[3*i+2] - glob_phys_iz1 + npoint_ghosts;
        index_moment[nmoment] = i;
        nmoment++;
        //fprintf(stdout,"myid is %d,si,sj,sk is %d,%d,%d\n",myid,moment_local_index[3*i+0],moment_local_index[3*i+1],moment_local_index[3*i+2]);
      }
      else
      {
        //default local index is -1
        moment_local_index[3*i+0] = -1;
        moment_local_index[3*i+1] = -1;
        moment_local_index[3*i+2] = -1;
      }
    }
  }

  if(num_force > 0)
  {
    free(force_coords);
    free(force_sur);
    free(force_global_index);
  }
  if(num_moment > 0)
  {
    free(moment_coords);
    free(moment_sur);
    free(moment_global_index);
  }

  float **force_value = NULL;
  float *force_wavelet_tstart = NULL;
  float **moment_value = NULL;
  float *moment_wavelet_tstart = NULL;
  char  **moment_wavelet_mechism = NULL;
  if(nforce>0)
  {
    force_value = (float **) fdlib_mem_malloc_2l_float(FD_NDIM, nforce*nt_in, "force_value");
    force_wavelet_tstart = (float *) malloc(nforce * sizeof(float));
  }
  if(nmoment>0)
  { 
    moment_value = (float **) fdlib_mem_malloc_2l_float(6, nmoment*nt_in, "moment_value");
    moment_wavelet_tstart = (float *) malloc(nmoment * sizeof(float));
    moment_wavelet_mechism = (char **) fdlib_mem_malloc_2l_char(nmoment,100,"mechanism type:tensor or angel"); 
  }
  // read force
  if (num_force>0)
  {
    for (int i=0; i<num_force; i++)
    {	
      int indx = 0;
      if (force_local_index[3*i+0] >= 0)
      {  
        fgets(str,500,fp);
        sscanf(str,"%f",&force_wavelet_tstart[indx]);
        for(int k=0; k<nt_in; k++)
        {
          fgets(str,500,fp);
          sscanf(str,"%f %f %f",&force_value[0][indx*nt_in+k],&force_value[1][indx*nt_in+k],&force_value[2][indx*nt_in+k]);
        }
        indx++;
      }

      if(force_local_index[3*i+0] == -1)
      {
        fgets(str,500,fp);
        for(int k=0; k<nt_in; k++)
        {
          fgets(str,500,fp);
        }
      }  
    }
  }

  //read moment
  if (num_moment>0)
  {
    for(int i=0; i<num_moment; i++)
    {
      int indx = 0;
      if (moment_local_index[3*i+0] >= 0)   //force is in this thread
      {
        fgets(str,500,fp);
        sscanf(str,"%f",&moment_wavelet_tstart[indx]);

        fgets(str,500,fp);
        sscanf(str,"%s",moment_wavelet_mechism[indx]);

        if (strcmp("mechanism_angle",moment_wavelet_mechism[indx]) == 0)
        {
          for (int k=0; k<nt_in; k++)
          {
            float angel[6];
            float temp_moment[6];
            fgets(str,500,fp);
            sscanf(str,"%f %f %f %f %f %f",&angel[0],&angel[1],&angel[2],&angel[3],&angel[4],&angel[5]);
            // M0 = u*D*A; 
            float M0 = angel[3]*angel[4]*angel[5];  
            angle2moment(angel[0],angel[1],angel[2],temp_moment);
            moment_value[0][indx*nt_in + k] = M0*temp_moment[0];
            moment_value[1][indx*nt_in + k] = M0*temp_moment[1];
            moment_value[2][indx*nt_in + k] = M0*temp_moment[2];
            moment_value[3][indx*nt_in + k] = M0*temp_moment[3];
            moment_value[4][indx*nt_in + k] = M0*temp_moment[4];
            moment_value[5][indx*nt_in + k] = M0*temp_moment[5];
          }
        }
        if (strcmp("moment_tensor",moment_wavelet_mechism[indx]) == 0)
        {
          for (int k=0;k<nt_in;k++)
          {
            fgets(str,500,fp);
            sscanf(str,"%f %f %f %f %f %f",&moment_value[0][indx*nt_in+k], &moment_value[1][indx*nt_in+k], 
                   &moment_value[2][indx*nt_in+k], &moment_value[3][indx*nt_in+k],
                   &moment_value[4][indx*nt_in+k], &moment_value[5][indx*nt_in+k]);
          }
        }
        indx++;
      }

      if (moment_local_index[3*i+0] == -1)
      {
        fgets(str,500,fp);
        fgets(str,500,fp);
        for (int k=0;k<nt_in;k++)
        {
          fgets(str,500,fp);
        }
      }  
    }  
  }
  fclose(fp); 

  // get total elem of exted src region for a single point
  //int siz_ext = 7 * 7 * 7;
  int siz_ext = (2*npoint_half_ext+1)*(2*npoint_half_ext+1)*(2*npoint_half_ext+1);
  float *t_in = (float *)malloc(nt_in*sizeof(float));
  // set force
  int   *force_info = NULL;
  float *force_vec_stf = NULL;
  float *force_ext_coef = NULL;
  int   *force_ext_indx = NULL;
  int   nt_force = 0;
  if(nforce > 0)
  {
    // allocate info and return to main function
    force_info = (int *)fdlib_mem_calloc_1d_int(nforce*M_SRC_INFO_NVAL, 1,"src_read_locat_valsrc");
    // stf
    nt_force = (int) (((nt_in-1)*dt_in / dt) + 0.5) + 1;
    force_vec_stf = (float *)fdlib_mem_calloc_1d_float(nforce*nt_force*FD_NDIM*num_of_stages,0.0,"src_read_locat_valsrc");
    // ext
    force_ext_coef = (float *)malloc(nforce*siz_ext * sizeof(float));
    force_ext_indx = (int   *)malloc(nforce*siz_ext * sizeof(int  ));
  }
  for (int i = 0 ; i < nforce; i++)
  {
    // convert time to index
    int  it_begin = (int) (force_wavelet_tstart[i] / dt);
    int  it_end   = it_begin + nt_force - 1;
    int  indx     = index_force[i]; 
    int  si = force_local_index[3 * indx + 0];
    int  sj = force_local_index[3 * indx + 1];
    int  sk = force_local_index[3 * indx + 2];
    // save values to inner var
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SI   ] = si;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SJ   ] = sj;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SK   ] = sk;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS  ] = 0 + i*nt_force*FD_NDIM*num_of_stages;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG] = it_begin;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND] = it_end;

    int ipos = force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS];
    int it1  = force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG];
    int it2  = force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND];
    float *this_vec_stf = force_vec_stf + ipos;
    float *this_force_value;
    int order = 3;
    for(int it_in=0; it_in<nt_in; it_in++)
    {
      t_in[it_in] = force_wavelet_tstart[i] + it_in*dt_in;
    } 

    for (int icmp=0; icmp<FD_NDIM; icmp++)
    {
      this_force_value = *(force_value+icmp) + i * nt_in;
      for (int it=it1; it<=it2; it++)
      {
        int it_to_it1 = (it - it1);
        for (int istage=0; istage<num_of_stages; istage++)
        {
          int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_force,num_of_stages);
          float t = it * dt + t0 + rk_stage_time[istage] * dt;
          // interp1d give t to get this_vec_stf
          this_vec_stf[iptr] = LagInterp_Piecewise_1d(t_in, this_force_value, nt_in, order, force_wavelet_tstart[i], dt_in, t);
          fprintf(stdout,"value is %f\n",this_vec_stf[iptr]);
        }
      }     
    }
    float *this_force_ext_coef = force_ext_coef + i * siz_ext;
    cal_norm_delt3d(this_force_ext_coef, force_index_inc[3*indx+0], force_index_inc[3*indx+1], force_index_inc[3*indx+2], 1.5, 1.5, 1.5, 3);

    size_t iptr_s = 0;
    // force_local_index si,sj,sk
    for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
    {
      if (k<nk1 || k>nk2) continue;

      for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
      {
        if (j<nj1 || j>nj2) continue;

        for (int ii=si-npoint_half_ext; ii<=si+npoint_half_ext; ii++)
        {
          if (ii<ni1 || ii>ni2) continue;

          int iptr = ii + j * siz_line + k * siz_slice;
          int iptr1 = (ii-(si-npoint_half_ext)) + 7 * (j-(sj-npoint_half_ext)) + 7 * 7 *(k-(sk-npoint_half_ext));
          force_ext_indx[iptr_s + i * siz_ext] = iptr;
          force_ext_coef[iptr_s + i * siz_ext] = this_force_ext_coef[iptr1];
          iptr_s++;
        }
      }
    }

    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_MAX ] = siz_ext;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_THIS] = iptr_s;
  }
  *num_of_force = nforce;
  *p_force_info = force_info;
  *p_force_vec_stf = force_vec_stf;
  *p_force_ext_indx = force_ext_indx;
  *p_force_ext_coef = force_ext_coef;

  // set moment
  int   *moment_info = NULL;
  float *moment_ten_rate = NULL;
  float *moment_ext_coef = NULL;
  int   *moment_ext_indx = NULL;
  int   nt_moment = 0;
  if(nmoment > 0)
  {
    // allocate info and return to main function
    moment_info = (int *)fdlib_mem_calloc_1d_int(nmoment*M_SRC_INFO_NVAL, 1 , "src_gen_multiple_moment_point_gauss");
    // stf
    nt_moment = (int) (((nt_in-1)*dt_in / dt) + 0.5) + 1;
    moment_ten_rate = (float *)fdlib_mem_calloc_1d_float(nmoment*nt_moment*6*num_of_stages,0.0,
        "src_read_locat_valsrc");
    // ext
    moment_ext_coef = (float *)malloc(nmoment*siz_ext * sizeof(float));
    moment_ext_indx = (int   *)malloc(nmoment*siz_ext * sizeof(int  ));
  }
  for (int i = 0 ; i < nmoment; i++)
  {
    // convert time to index
    int  it_begin = (int) (moment_wavelet_tstart[i] / dt);
    int  it_end   = it_begin + nt_moment - 1;
    int  indx     = index_moment[i]; 
    int  si = moment_local_index[3 * indx + 0];
    int  sj = moment_local_index[3 * indx + 1];
    int  sk = moment_local_index[3 * indx + 2];
    // save values to inner var
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SI   ] = moment_local_index[3 * indx + 0];
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SJ   ] = moment_local_index[3 * indx + 1];
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SK   ] = moment_local_index[3 * indx + 2];
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS  ] = 0 + i*nt_moment*6*num_of_stages;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG] = it_begin;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND] = it_end;
    int ipos = moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS];
    int it1  = moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG];
    int it2  = moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND];
    float *this_ten_rate = moment_ten_rate + ipos;
    float *this_moment_value;
    int order = 3;
    for(int it_in=0; it_in<nt_in; it_in++)
    {
      t_in[it_in] = moment_wavelet_tstart[i] + it_in*dt_in;
    } 

    for (int icmp=0; icmp<6; icmp++)
    {
      this_moment_value = *(moment_value+icmp) + i * nt_in;
      for (int it=it1; it<=it2; it++)
      {
        int it_to_it1 = (it - it1);
        for (int istage=0; istage<num_of_stages; istage++)
        {
          int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_moment,num_of_stages);
          float t = it * dt + t0 + rk_stage_time[istage] * dt;
          // interp give t to get this_ten_rate
          // save to vector
          this_ten_rate[iptr] = LagInterp_Piecewise_1d(t_in, this_moment_value, nt_in, order, moment_wavelet_tstart[i], dt_in, t);
        }
      }
    }
    float *this_moment_ext_coef = moment_ext_coef + i * siz_ext;
    cal_norm_delt3d(this_moment_ext_coef, moment_index_inc[3*indx+0], moment_index_inc[3*indx+1], moment_index_inc[3*indx+2],1.5, 1.5, 1.5, 3);

    size_t iptr_s = 0;
    for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
    {
      if (k<nk1 || k>nk2) continue;

      for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
      {
        if (j<nj1 || j>nj2) continue;

        for (int ii=si-npoint_half_ext; ii<=si+npoint_half_ext; ii++)
        {
          if (ii<ni1 || ii>ni2) continue;

          int iptr = ii + j * siz_line + k * siz_slice;
          int iptr1 = (ii-(si-npoint_half_ext)) + 7 * (j-(sj-npoint_half_ext)) + 7 * 7 *(k-(sk-npoint_half_ext));
          moment_ext_indx[iptr_s + i * siz_ext] = iptr;
          moment_ext_coef[iptr_s + i * siz_ext] = this_moment_ext_coef[iptr1];
          iptr_s++;
        }
      }
    }
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_MAX ] = siz_ext;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_THIS] = iptr_s;
  }
  *num_of_moment = nmoment;
  *p_moment_info = moment_info;
  *p_moment_ten_rate = moment_ten_rate;
  *p_moment_ext_indx = moment_ext_indx;
  *p_moment_ext_coef = moment_ext_coef;

  if(num_force>0)
  {
    free(force_value);
    free(force_wavelet_tstart);
    free(index_force);
    free(force_local_index);
    free(force_index_inc);
  }
  if(num_moment>0)
  {
    free(moment_value);
    free(moment_wavelet_tstart);
    free(moment_wavelet_mechism);
    free(index_moment);
    free(moment_local_index);
    free(moment_index_inc);
  }
  return 0;
}


int
src_read_locate_anasrc(char *pfilepath,
                       size_t siz_line,
                       size_t siz_slice,
                       float t0,
                       float dt,
                       int   num_of_stages,
                       float *rk_stage_time,
                       int   glob_phys_ix1, // gloabl start index along x this thread
                       int   glob_phys_ix2, // gloabl end index along x
                       int   glob_phys_iy1,
                       int   glob_phys_iy2,
                       int   glob_phys_iz1,
                       int   glob_phys_iz2,
                       int   ni1,
                       int   ni2,
                       int   nj1,
                       int   nj2,
                       int   nk1,
                       int   nk2,
                       int   npoint_half_ext,
                       int   npoint_ghosts,
                       float *x3d,
                       float *y3d,
                       float *z3d,
                       MPI_Comm comm,
                       int myid,
                       // following output
                       char **p_event_name,
                       int   *num_of_force, // force in this thread
                       int **restrict p_force_info,
                       float  **restrict p_force_vec_stf,
                       int    **restrict p_force_ext_indx,
                       float  **restrict p_force_ext_coef,
                       int   *num_of_moment, // moment in this thread
                       int    **restrict p_moment_info,
                       float  **restrict p_moment_ten_rate,
                       int    **restrict p_moment_ext_indx,
                       float  **restrict p_moment_ext_coef,
                       int verbose)
{
  FILE* fp =NULL;
  char str[500];
  int num_force;
  int num_moment;
  float wavelet_tlen;
  char *event_name = (char *)malloc(500 * sizeof(char));
  // read analysis source from inputfile
  if ((fp = fopen(pfilepath, "r"))==NULL) fprintf(stdout,"fail to open");
  fgets(str,500,fp);
  sscanf(str,"%s",event_name);
  fgets(str,500,fp);
  sscanf(str,"%d %d",&num_force,&num_moment);
  fgets(str,500,fp);
  sscanf(str,"%f",&wavelet_tlen);
  *p_event_name = event_name;
  float *force_coords = NULL;
  int   *force_sur = NULL;
  float *force_vector = NULL;
  float *force_wavelet_coefs = NULL;
  float *force_eight_paras = NULL;
  float *force_wavelet_tstart = NULL;
  char  **force_wavelet_name = NULL;
  float *moment_coords = NULL;
  int   *moment_sur = NULL;
  float *moment_tensor = NULL;
  float *moment_wavelet_coefs = NULL;
  float *moment_eight_paras = NULL; 
  float *moment_wavelet_tstart = NULL;
  char  **moment_wavelet_name = NULL;
  char  **moment_wavelet_mechism = NULL;
  if(num_force>0)
  {
    force_coords = (float *) malloc(3 * num_force * sizeof(float));
    force_sur = (int *) malloc(num_force * sizeof(int)); // relocated_to_surface
    force_vector = (float *) malloc(3 * num_force * sizeof(float)); // fx fy fz
    force_wavelet_coefs = (float *) malloc(2 * num_force * sizeof(float)); // 2 paras of src function
    force_eight_paras = (float *) malloc(8 * num_force * sizeof(float)); // rest eight paras for function
    force_wavelet_tstart = (float *) malloc(num_force * sizeof(float));
    force_wavelet_name = (char **) fdlib_mem_malloc_2l_char(num_force,100,"source_time_function type");
  }
  if(num_moment>0)
  {
    moment_coords = (float *) malloc(3 * num_moment * sizeof(float));
    moment_sur = (int *) malloc(num_moment * sizeof(int)); 
    moment_tensor = (float *) malloc(6 * num_moment * sizeof(float)); // Mxx Myy Mzz Mxz Myz Mxy
    moment_wavelet_coefs = (float *) malloc(2 * num_moment * sizeof(float)); 
    moment_eight_paras = (float *) malloc(8 * num_moment * sizeof(float)); 
    moment_wavelet_tstart = (float *) malloc(num_moment * sizeof(float));
    moment_wavelet_name = (char **) fdlib_mem_malloc_2l_char(num_moment,100,"source_time_function type");
    moment_wavelet_mechism = (char **) fdlib_mem_malloc_2l_char(num_moment,100,"mechanism type:tensor or angle");
  }
  // read force
  if (num_force > 0)
  {
    for (int i=0;i<num_force;i++)
    {
      fgets(str,500,fp);
      sscanf(str,"%f %f %f %d",&force_coords[3*i+0],&force_coords[3*i+1],&force_coords[3*i+2],&force_sur[i]);

      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&force_vector[3*i+0], &force_vector[3*i+1], &force_vector[3*i+2]);

      fgets(str,500,fp);
      sscanf(str,"%s",force_wavelet_name[i]);

      fgets(str,500,fp);
      sscanf(str,"%f %f",&force_wavelet_coefs[2*i+0], &force_wavelet_coefs[2*i+1]);

      fgets(str,500,fp);
      sscanf(str,"%f %f %f %f %f %f %f %f",&force_eight_paras[8*i+0],&force_eight_paras[8*i+1],
          &force_eight_paras[8*i+2],&force_eight_paras[8*i+3],&force_eight_paras[8*i+4],
          &force_eight_paras[8*i+5],&force_eight_paras[8*i+6],&force_eight_paras[8*i+7]);

      fgets(str,500,fp);
      sscanf(str,"%f",&force_wavelet_tstart[i]);
    }
  }
  // read moment
  if (num_moment > 0)
  {
    for(int i=0;i<num_moment;i++)
    {
      fgets(str,500,fp);
      sscanf(str,"%f %f %f %d",&moment_coords[3*i+0],&moment_coords[3*i+1],
          &moment_coords[3*i+2],&moment_sur[i]);

      fgets(str,500,fp);
      sscanf(str, "%s", moment_wavelet_mechism[i]);

      if (strcmp("mechanism_angle",moment_wavelet_mechism[i])==0)
      {
        float angle[6];
        float temp_moment[6];
        fgets(str,500,fp);
        sscanf(str,"%f %f %f %f %f %f",&angle[0],&angle[1],&angle[2],&angle[3],&angle[4],&angle[5]);
        // M0 = mu * D * A
        float M0 = angle[3]*angle[4]*angle[5];    
        angle2moment(angle[0],angle[1],angle[2],temp_moment);
        for (int j=0;j<6;j++)
        {
          moment_tensor[6*i +j]=M0*temp_moment[j];
          fprintf(stdout,"moment_tensor is %f \n",moment_tensor[6*i +j]);
        }
      } 
      if (strcmp("moment_tensor",moment_wavelet_mechism[i])==0)
      {  
        fgets(str,500,fp);
        sscanf(str,"%f %f %f %f %f %f",&moment_tensor[6*i +0],&moment_tensor[6*i +1],
            &moment_tensor[6*i +2],&moment_tensor[6*i +3],&moment_tensor[6*i +4],&moment_tensor[6*i +5]);
      } 
      fgets(str,500,fp);
      sscanf(str,"%s",moment_wavelet_name[i]);

      fgets(str,500,fp);
      sscanf(str,"%f %f",&moment_wavelet_coefs[2*i+0], &moment_wavelet_coefs[2*i+1]);

      fgets(str,500,fp);
      sscanf(str,"%f %f %f %f %f %f %f %f",&moment_eight_paras[8*i+0],&moment_eight_paras[8*i+1],&moment_eight_paras[8*i+2],
          &moment_eight_paras[8*i+3],&moment_eight_paras[8*i+4],&moment_eight_paras[8*i+5],
          &moment_eight_paras[8*i+6],&moment_eight_paras[8*i+7]);

      fgets(str,500,fp);
      sscanf(str,"%f",&moment_wavelet_tstart[i]);
    }
  }
  fclose(fp);

  // get total elem of exted src region for a single point
  // int siz_ext = 7 * 7 * 7;
  int siz_ext = (2*npoint_half_ext+1)*(2*npoint_half_ext+1)*(2*npoint_half_ext+1);

  int nx = (ni2-ni1+1)+2*npoint_ghosts;
  int ny = (nj2-nj1+1)+2*npoint_ghosts;
  int nz = (nk2-nk1+1)+2*npoint_ghosts;  
  // workspace 3d var for distance calculation
  float *wrk3d = (float *) fdlib_mem_calloc_1d_float(nx*ny*nz,0.0, "src_read_locat_valsrc");
  int   nforce = 0;  //default number of force in this thread
  int   *force_global_index = NULL;
  int   *force_local_index = NULL;
  int   *index_force =NULL;
  float *force_index_inc =NULL;
  if(num_force>0)
  {
    // judge force in this thread
    force_global_index = (int*)malloc(3*num_force*sizeof(int));
    force_index_inc = (float*)malloc(3*num_force*sizeof(float));
    for (int i =0; i<num_force; i++)
    {
      fprintf(stdout,"locate force by coord  ...\n"); 
      // default global index and relative shift to -1 and 0.0

      int si,sj,sk;
      float sx_inc = 0.0; float sy_inc = 0.0; float sz_inc = 0.0;
      int sgpi=-1; int sgpj=-1; int sgpk=-1;
      // if located in this thread
      int is_here = src_coord2index(force_coords[3*i+0],force_coords[3*i+1],force_coords[3*i+2],
                                    nx, ny, nz, 
                                    ni1,ni2,nj1,nj2,nk1,nk2,
                                    x3d, y3d, z3d, wrk3d,
                                    &si, &sj, &sk,
                                    &sx_inc, &sy_inc, &sz_inc);
      if ( is_here == 1)
      {
        // conver to global index
        sgpi = si - npoint_ghosts + glob_phys_ix1;
        sgpj = sj - npoint_ghosts + glob_phys_iy1;
        sgpk = sk - npoint_ghosts + glob_phys_iz1;
        fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
        fprintf(stdout," -- located to global index = %d %d %d\n", sgpi,sgpj,sgpk);
        fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
      } else {
        fprintf(stdout," this force not in this thread %d\n", myid);
      }
      force_index_inc[3*i+0] = sx_inc;
      force_index_inc[3*i+1] = sy_inc;
      force_index_inc[3*i+2] = sz_inc;
      force_global_index[3*i+0] = sgpi;
      force_global_index[3*i+1] = sgpj;
      force_global_index[3*i+2] = sgpk;
    }

    // reduce all source points global index and shift values
    // Note select MPI_MAX for global index and MPI_SUM for shift value.
    // 0 -> x, 1->y, 2->z.
    int *reduce_force_global_index = (int*)malloc(3*num_force*sizeof(int));
    MPI_Allreduce(force_global_index, reduce_force_global_index, 3*num_force, MPI_INT, MPI_MAX, comm);

    float *reduce_force_index_inc = (float*)malloc(3*num_force*sizeof(float));
    MPI_Allreduce(force_index_inc,reduce_force_index_inc, 3*num_force, MPI_INT, MPI_SUM, comm);

    for (int i=0; i<num_force; i++)
    {
      fprintf(stdout,"force %d --myid=%d,index=%d %d %d,shift = %f %f %f\n",
          i, myid,reduce_force_global_index[3*i+0],reduce_force_global_index[3*i+1],reduce_force_global_index[3*i+2],
          reduce_force_index_inc[3*i+0],reduce_force_index_inc[3*i+1],reduce_force_index_inc[3*i+2]);
    }

    index_force = (int *)malloc(num_force*sizeof(int));
    force_local_index = (int*)malloc(3*num_force*sizeof(int));
    for (int i=0; i<num_force; i++)
    {
      // use grid index to check if in this thread after extend
      if( reduce_force_global_index[3*i+0]-npoint_half_ext <= glob_phys_ix2 && // exted left point is less than right bdry
          reduce_force_global_index[3*i+0]+npoint_half_ext >= glob_phys_ix1 && // exted right point is larger than left bdry
          reduce_force_global_index[3*i+1]-npoint_half_ext <= glob_phys_iy2 && 
          reduce_force_global_index[3*i+1]+npoint_half_ext >= glob_phys_iy1 &&
          reduce_force_global_index[3*i+2]-npoint_half_ext <= glob_phys_iz2 && 
          reduce_force_global_index[3*i+2]+npoint_half_ext >= glob_phys_iz1)
      {
        // at least one extend point in this thread
        // convert to local index
        force_local_index[3*i+0] = reduce_force_global_index[3*i+0] - glob_phys_ix1 + npoint_ghosts;
        force_local_index[3*i+1] = reduce_force_global_index[3*i+1] - glob_phys_iy1 + npoint_ghosts;
        force_local_index[3*i+2] = reduce_force_global_index[3*i+2] - glob_phys_iz1 + npoint_ghosts;
        index_force[nforce] = i;
        nforce++;
        //fprintf(stdout,"myid is %d,si,sj,sk is %d,%d,%d\n",myid,force_local_index[3*i+0],force_local_index[3*i+1],force_local_index[3*i+2]);
      }
      else
      {
        //default local index is -1
        force_local_index[3*i+0] = -1;
        force_local_index[3*i+1] = -1; 
        force_local_index[3*i+2] = -1;
      }
    }
  }
  int *force_info = NULL;
  float *force_vec_stf = NULL;
  float *force_ext_coef = NULL;
  int   *force_ext_indx = NULL;
  int nt_force = 0;
  // set force
  if(nforce > 0)
  {
    // allocate info and return to main function
    force_info = (int *)fdlib_mem_calloc_1d_int(nforce*M_SRC_INFO_NVAL, 1, "src_read_locat_anasrc");
    // stf
    nt_force = (int) (wavelet_tlen / dt +0.5) + 1;
    force_vec_stf = (float *)fdlib_mem_calloc_1d_float(nforce*nt_force*FD_NDIM*num_of_stages,0.0,
        "src_read_loact_anasrc");
    // ext
    force_ext_coef = (float *)malloc(nforce*siz_ext * sizeof(float));
    force_ext_indx = (int   *)malloc(nforce*siz_ext * sizeof(int  ));
  }
  for (int i = 0 ; i < nforce; i++)
  {
    int  indx     = index_force[i]; 
    // convert time to index
    int  it_begin = (int) (force_wavelet_tstart[indx] / dt);
    int  it_end   = it_begin + nt_force - 1;
    int  si = force_local_index[3 * indx + 0];
    int  sj = force_local_index[3 * indx + 1];
    int  sk = force_local_index[3 * indx + 2];
    // save values to inner var
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SI   ] = si;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SJ   ] = sj;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SK   ] = sk;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS  ] = 0 + i * nt_force*FD_NDIM*num_of_stages;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG] = it_begin;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND] = it_end;

    int ipos = force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS];
    int it1  = force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG];
    int it2  = force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND];
    float *this_vec_stf = force_vec_stf + ipos;
    float *this_force_vector = force_vector + indx*3; 
    for (int icmp=0; icmp<FD_NDIM; icmp++)
    {
      for (int it=it1; it<=it2; it++)
      {
        int it_to_it1 = (it - it1);
        for (int istage=0; istage<num_of_stages; istage++)
        {
          int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_force,num_of_stages);
          float t = it * dt + t0 + rk_stage_time[istage] * dt;
          float stf_val;
          if (strcmp(force_wavelet_name[indx], "ricker")==0) {
            stf_val = fun_ricker(t, force_wavelet_coefs[2*indx+0], force_wavelet_coefs[2*indx+1]);
          } else if (strcmp(force_wavelet_name[indx], "gaussian")==0) {
            stf_val = fun_gauss(t, force_wavelet_coefs[2*indx+0], force_wavelet_coefs[2*indx+1]);
          } else {
            fprintf(stderr,"wavelet_name=%s\n", force_wavelet_name[indx]); 
            fprintf(stderr,"   not implemented yet\n"); 
            fflush(stderr);
          }  
          // save to vector   
          this_vec_stf[iptr] = stf_val * this_force_vector[icmp];
        }
      }     
    }

    float *this_force_ext_coef = force_ext_coef + i * siz_ext;
    cal_norm_delt3d(this_force_ext_coef, force_index_inc[3*indx+0], force_index_inc[3*indx+1], force_index_inc[3*indx+2], 1.5, 1.5, 1.5, 3);

    size_t iptr_s = 0;
    for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
    {
      if (k<nk1 || k>nk2) continue;

      for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
      {
        if (j<nj1 || j>nj2) continue;

        for (int ii=si-npoint_half_ext; ii<=si+npoint_half_ext; ii++)
        {
          if (ii<ni1 || ii>ni2) continue;

          // Note index need match coef
          int iptr = ii + j * siz_line + k * siz_slice;
          int iptr1 = (ii-(si-npoint_half_ext)) + 7 * (j-(sj-npoint_half_ext)) + 7 * 7 *(k-(sk-npoint_half_ext));
          force_ext_indx[iptr_s + i * siz_ext] = iptr;
          force_ext_coef[iptr_s + i * siz_ext] = this_force_ext_coef[iptr1];
          iptr_s++;
        }
      }
    }

    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_MAX ] = siz_ext;
    force_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_THIS] = iptr_s;
  }
  *num_of_force = nforce;
  *p_force_info = force_info;
  *p_force_vec_stf = force_vec_stf;
  *p_force_ext_indx = force_ext_indx;
  *p_force_ext_coef = force_ext_coef;

  if(num_force>0)
  {
    free(force_coords);
    free(force_sur);
    free(force_vector);
    free(force_wavelet_coefs);
    free(force_eight_paras);
    free(force_wavelet_tstart);
    free(force_wavelet_name);
    free(index_force);
    free(force_global_index);
    free(force_local_index);
    free(force_index_inc);
  }

  int *moment_global_index = NULL; 
  int nmoment = 0;
  int *index_moment =NULL;
  int *moment_local_index = NULL;
  float *moment_index_inc = NULL;

  if(num_moment>0)
  {
    // judge moment in this thread
    moment_global_index = (int*)malloc(3*num_moment*sizeof(int));
    moment_index_inc = (float*)malloc(3*num_moment*sizeof(float));
    for (int i =0; i<num_moment; i++)
    {
      fprintf(stdout,"locate moment by coord  ...\n"); 
      // default global index and relative shift to -1 and 0
      int si, sj, sk;
      float sx_inc = 0.0; float sy_inc = 0.0; float sz_inc = 0.0;
      int sgpi=-1; int sgpj=-1; int sgpk=-1;
      // if located in this thread
      int is_here = src_coord2index(moment_coords[3*i+0],moment_coords[3*i+1],moment_coords[3*i+2],
                                    nx, ny, nz, 
                                    ni1,ni2,nj1,nj2,nk1,nk2,
                                    x3d, y3d, z3d, wrk3d,
                                    &si, &sj, &sk,
                                    &sx_inc, &sy_inc, &sz_inc);
      if ( is_here == 1)
      {
        // conver to global index
        sgpi = si - npoint_ghosts + glob_phys_ix1;
        sgpj = sj - npoint_ghosts + glob_phys_iy1;
        sgpk = sk - npoint_ghosts + glob_phys_iz1;
        fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
        fprintf(stdout," -- located to global index = %d %d %d\n", sgpi,sgpj,sgpk);
        fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
      } else {
        fprintf(stdout," this moment not in this thread %d\n", myid);
      }
      moment_index_inc[3*i+0] = sx_inc;
      moment_index_inc[3*i+1] = sy_inc;
      moment_index_inc[3*i+2] = sz_inc;
      moment_global_index[3*i+0] = sgpi;
      moment_global_index[3*i+1] = sgpj;
      moment_global_index[3*i+2] = sgpk;
    }

    // reduce all source points global index and shift values
    // Note select MPI_MAX for global index and MPI_SUM for shift value.
    // 0 -> x, 1->y, 2->z.
    int *reduce_moment_global_index = (int*)malloc(3*num_moment*sizeof(int));
    MPI_Allreduce(moment_global_index, reduce_moment_global_index, 3*num_moment, MPI_INT, MPI_MAX, comm);

    float *reduce_moment_index_inc = (float*)malloc(3*num_moment*sizeof(float));
    MPI_Allreduce(moment_index_inc,reduce_moment_index_inc, 3*num_moment, MPI_INT, MPI_SUM, comm);

    for (int i=0; i<num_moment; i++)
    {
      fprintf(stdout,"moment %d --myid=%d,index=%d %d %d,shift = %f %f %f\n",
          i, myid,reduce_moment_global_index[3*i+0],reduce_moment_global_index[3*i+1],reduce_moment_global_index[3*i+2],
          reduce_moment_index_inc[3*i+0],reduce_moment_index_inc[3*i+1],reduce_moment_index_inc[3*i+2]);
    }

    index_moment = (int *)malloc(num_moment*sizeof(int));
    moment_local_index = (int*)malloc(3*num_moment*sizeof(int));
    // use grid index to check if in this thread after extend
    for(int i=0; i<num_moment; i++)
    {
      if (reduce_moment_global_index[3*i+0]-npoint_half_ext <= glob_phys_ix2 && // exted left point is less than right bdry
          reduce_moment_global_index[3*i+0]+npoint_half_ext >= glob_phys_ix1 && // exted right point is larger than left bdry
          reduce_moment_global_index[3*i+1]-npoint_half_ext <= glob_phys_iy2 && 
          reduce_moment_global_index[3*i+1]+npoint_half_ext >= glob_phys_iy1 &&
          reduce_moment_global_index[3*i+2]-npoint_half_ext <= glob_phys_iz2 && 
          reduce_moment_global_index[3*i+2]+npoint_half_ext >= glob_phys_iz1)
      {
        // at least one extend point in this thread
        // convert to local index
        moment_local_index[3*i+0]= reduce_moment_global_index[3*i+0] - glob_phys_ix1 + npoint_ghosts;
        moment_local_index[3*i+1]= reduce_moment_global_index[3*i+1] - glob_phys_iy1 + npoint_ghosts;
        moment_local_index[3*i+2]= reduce_moment_global_index[3*i+2] - glob_phys_iz1 + npoint_ghosts;
        index_moment[nmoment] = i;
        nmoment++;
        fprintf(stdout,"myid is %d,si,sj,sk is %d,%d,%d\n",myid,moment_local_index[3*i+0],moment_local_index[3*i+1],moment_local_index[3*i+2]);
      }
      else
      {
        //default local index is -1
        moment_local_index[3*i+0] = -1;
        moment_local_index[3*i+1] = -1;
        moment_local_index[3*i+2] = -1;
      }
    }
  }
  // set moment
  int   *moment_info = NULL;
  float *moment_ten_rate = NULL;
  float *moment_ext_coef = NULL;
  int   *moment_ext_indx = NULL;
  int   nt_moment = 0;
  if(nmoment > 0)
  {
    // allocate info and return to main function
    moment_info = (int *)fdlib_mem_calloc_1d_int(nmoment*M_SRC_INFO_NVAL, 1 , "src_read_loact_anasrc");
    // stf
    nt_moment = (int) (wavelet_tlen / dt + 0.5) + 1;
    moment_ten_rate = (float *)fdlib_mem_calloc_1d_float(nmoment*nt_moment*6*num_of_stages,0.0,
        "src_read_locat_anasrc");
    // ext
    moment_ext_coef = (float *)malloc(nmoment*siz_ext * sizeof(float));
    moment_ext_indx = (int   *)malloc(nmoment*siz_ext * sizeof(int  ));
  }
  for (int i = 0 ; i < nmoment; i++)
  {
    int  indx     = index_moment[i];
    // convert time to index
    int  it_begin = (int) (moment_wavelet_tstart[indx] / dt);
    int  it_end   = it_begin + nt_moment - 1;
    int  si = moment_local_index[3 * indx + 0];
    int  sj = moment_local_index[3 * indx + 1];
    int  sk = moment_local_index[3 * indx + 2];
    // save values to inner var
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SI   ] = si;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SJ   ] = sj;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_SK   ] = sk;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS  ] = 0 + i*nt_moment*6*num_of_stages;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG] = it_begin;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND] = it_end;

    int ipos = moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_POS];
    int it1  = moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITBEG];
    int it2  = moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_ITEND];
    float *this_ten_rate = moment_ten_rate + ipos;
    float *this_moment_tensor = moment_tensor + indx*6; 
    for (int icmp=0; icmp<6; icmp++)
    {
      for (int it=it1; it<=it2; it++)
      {
        int it_to_it1 = (it - it1);
        for (int istage=0; istage<num_of_stages; istage++)
        {
          int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_moment,num_of_stages);
          float t = it * dt + t0 + rk_stage_time[istage] * dt;
          float stf_val;
          if (strcmp(moment_wavelet_name[indx], "ricker_deriv")==0) {
            stf_val = fun_ricker_deriv(t, moment_wavelet_coefs[2*indx+0], moment_wavelet_coefs[2*indx+1]);
          } else if (strcmp(moment_wavelet_name[indx], "gaussian_deriv")==0) {
            stf_val = fun_gauss_deriv(t, moment_wavelet_coefs[2*indx+0], moment_wavelet_coefs[2*indx+1]);
          } else {
            fprintf(stderr,"wavelet_name=%s\n", moment_wavelet_name[indx]); 
            fprintf(stderr,"   not implemented yet\n"); 
            fflush(stderr);
          }
          // save to vector
          this_ten_rate[iptr] = stf_val * this_moment_tensor[icmp];
        }
      }
    }
    float *this_moment_ext_coef = moment_ext_coef + i * siz_ext;
    cal_norm_delt3d(this_moment_ext_coef, moment_index_inc[3*indx+0], moment_index_inc[3*indx+1], moment_index_inc[3*indx+2],1.5, 1.5, 1.5, 3);

    size_t iptr_s = 0;
    for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
    {
      if (k<nk1 || k>nk2) continue;

      for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
      {
        if (j<nj1 || j>nj2) continue;

        for (int ii=si-npoint_half_ext; ii<=si+npoint_half_ext; ii++)
        {
          if (ii<ni1 || ii>ni2) continue;

          // Note index need match coef
          int iptr = ii + j * siz_line + k * siz_slice;
          int iptr1 = (ii-(si-npoint_half_ext)) + 7 * (j-(sj-npoint_half_ext)) + 7 * 7 *(k-(sk-npoint_half_ext));
          moment_ext_indx[iptr_s + i * siz_ext] = iptr;
          moment_ext_coef[iptr_s + i * siz_ext] = this_moment_ext_coef[iptr1];
          iptr_s++;
        }
      }
    }
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_MAX ] = siz_ext;
    moment_info[M_SRC_INFO_NVAL * i + M_SRC_INFO_SEQ_NEXT_THIS] = iptr_s;
  }
  *num_of_moment = nmoment;
  *p_moment_info = moment_info;
  *p_moment_ten_rate = moment_ten_rate;
  *p_moment_ext_indx = moment_ext_indx;
  *p_moment_ext_coef = moment_ext_coef;
  if(num_moment>0)
  {
    free(moment_coords);
    free(moment_sur);
    free(moment_tensor);
    free(moment_wavelet_coefs);
    free(moment_eight_paras);
    free(moment_wavelet_tstart);
    free(moment_wavelet_name);
    free(moment_wavelet_mechism);
    free(index_moment);
    free(moment_local_index);
    free(moment_index_inc);
    free(moment_global_index);
  }
  return 0;
}


/*
 * 3d spatial smoothing
 */

  void
cal_norm_delt3d(float *delt, float x0, float y0, float z0, float rx0, float ry0, float rz0, int LenDelt)
{
  float SUM = 0.0 ;

  int iptr = 0;
  for(int k=-LenDelt; k<=LenDelt; k++) {
    for(int j=-LenDelt; j<=LenDelt; j++) {
      for(int i=-LenDelt; i<=LenDelt; i++) {
        float D1 = fun_gauss(i-x0, rx0 ,0.0);           
        float D2 = fun_gauss(j-y0, ry0 ,0.0);          
        float D3 = fun_gauss(k-z0, rz0 ,0.0);          
        delt[iptr] = D1 * D2 * D3;
        SUM += delt[iptr];
        iptr++;
      }
    }               
  }

  if( SUM < 1e-20 )
  {
    fprintf(stderr, "cal_norm_delt is zero\n");
    exit(1);
  }

  int siz_1d = 2 * LenDelt + 1;
  for (int iptr=0; iptr< siz_1d*siz_1d*siz_1d; iptr++) {
    delt[iptr] /= SUM;
  }
} 

/*
 * wavelet functions
 */

// ricker and it deriv.
  float 
fun_ricker(float t, float fc, float t0)
{
  //float pi = acos(-1.0);
  float f0 = sqrtf(PI)/2.0;
  float u = (t-t0)*2.0*PI*fc;
  float v = (u*u/4-0.5)*exp(-u*u/4)*f0;

  return v;
}

  float 
fun_ricker_deriv(float t, float fc, float t0)
{
  //float pi = acos(-1.0);
  float f0 = sqrtf(PI)/2.0;
  float u = (t-t0)*2.0*PI*fc;
  float v = u*(1.5-u*u/4)*exp(-u*u/4)*f0*PI*fc;

  return v;
}
//gauss and it deriv
  float
fun_gauss(float t, float a, float t0)
{
  float f;
  f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrtf(PI)*a);
  return f;
}

  float
fun_gauss_deriv(float t, float a, float t0)
{
  float f;
  f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrtf(PI)*a)*(-2*(t-t0)/(a*a));
  return f;
}


/*
 * get the stf and moment rate for one stage
 */

  void
src_get_stage_stf(
    int num_of_force,
    int *restrict force_info, // num_of_force * 6 : si,sj,sk,start_pos_in_stf,start_it, end_it
    float *restrict force_vec_stf,
    int num_of_moment,
    int *restrict moment_info, // num_of_force * 6 : si,sj,sk,start_pos_in_rate,start_it, end_it
    float *restrict moment_ten_rate,
    int it, int istage, int num_of_stages,
    float *restrict force_vec_value,
    float *restrict moment_ten_value,
    const int myid, const int verbose)
{
  for (int n=0; n<num_of_force; n++)
  {
    int ipos = force_info[M_SRC_INFO_NVAL*n + M_SRC_INFO_SEQ_POS];
    int it1  = force_info[M_SRC_INFO_NVAL*n + M_SRC_INFO_SEQ_ITBEG];
    int it2  = force_info[M_SRC_INFO_NVAL*n + M_SRC_INFO_SEQ_ITEND];
    int nt_force = it2 - it1 + 1;

    // point tho this force in vec_stf
    float *ptr_force = force_vec_stf + ipos;

    for (int icmp=0; icmp<FD_NDIM; icmp++)
    {
      int iptr_value = n * FD_NDIM + icmp;
      if (it < it1 || it > it2)
      {
        force_vec_value[iptr_value] = 0.0;
      }
      else
      {
        int it_to_it1 = it - it1;
        int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_force,num_of_stages);

        force_vec_value[iptr_value] = ptr_force[iptr];
      }
    }
  }

  for (int n=0; n<num_of_moment; n++)
  {
    int ipos = moment_info[M_SRC_INFO_NVAL*n + M_SRC_INFO_SEQ_POS];
    int it1  = moment_info[M_SRC_INFO_NVAL*n + M_SRC_INFO_SEQ_ITBEG];
    int it2  = moment_info[M_SRC_INFO_NVAL*n + M_SRC_INFO_SEQ_ITEND];
    int nt_moment = it2 - it1 + 1;
    // point tho this moment in ten_rate
    float *ptr_moment = moment_ten_rate + ipos;

    for (int icmp=0; icmp<6; icmp++)
    {
      int iptr_value = n * 6 + icmp;

      if (it < it1 || it > it2)
      {
        moment_ten_value[iptr_value] = 0.0;
      }
      else
      {
        int it_to_it1 = it - it1;
        int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_moment,num_of_stages);
        moment_ten_value[iptr_value] = ptr_moment[iptr];
      }
    }
  }

  // for test, reset only at it=0
  //int n = 0;
  //for (int icmp=0; icmp<6; icmp++)
  //{
  //  int iptr_value = n * 6 + icmp;
  //  if (it==0 && icmp<3) {
  //    moment_ten_value[iptr_value] = 1.0e9;
  //  }
  //  else
  //  {
  //    moment_ten_value[iptr_value] = 0.0;
  //  }
  //}

}

/*
 * convert angles (defined as Aki and Richards) to moment tensor 
 *  in the cartesian coordinate: x-east, y-north, z-upward
 */

  void 
angle2moment(float strike, float dip, float rake, float* source_moment_tensor)
{
  float strike_pi,dip_pi,rake_pi; 
  float M11,M22,M33,M12,M13,M23;

  dip_pi    = dip    / 180.0 * PI; 
  strike_pi = strike / 180.0 * PI;
  rake_pi   = rake   / 180.0 * PI;

  // Angles are defined same as in Aki and Richard's book
  M11 = - (  sin(dip_pi) * cos(rake_pi) * sin(2.0*strike_pi) 
      + sin(2.0*dip_pi) * sin(rake_pi) * sin(strike_pi) * sin(strike_pi) );

  M22 =  sin(dip_pi) * cos(rake_pi) * sin(2.0 * strike_pi)     
    -sin(2.0*dip_pi) * sin(rake_pi) * cos(strike_pi) * cos(strike_pi) ;

  M33 = - ( M11 + M22 );

  M12 =   sin(dip_pi) * cos(rake_pi) * cos(2.0 * strike_pi)     
    + 0.5 * sin(2.0 * dip_pi) * sin(rake_pi) * sin(2.0 * strike_pi) ;

  M13 = - (  cos(dip_pi) * cos(rake_pi) * cos(strike_pi)  
      + cos(2.0 * dip_pi) * sin(rake_pi) * sin(strike_pi) ) ;

  M23 = - (  cos(dip_pi) * cos(rake_pi) * sin(strike_pi) 
      - cos(2.0*dip_pi) * sin(rake_pi) * cos(strike_pi) );

  // attention: the order may be different with outside
  // Mxz=-Mxz;Mxy=-Mxy !for upward positive z axis
  // x->2 y->1 z->3
  source_moment_tensor[0] =  M22 ;  // Mxx
  source_moment_tensor[1] =  M11 ;  // Myy 
  source_moment_tensor[2] =  M33 ;  // Mzz
  source_moment_tensor[3] = -M23 ;  // Mxz 
  source_moment_tensor[4] = -M13 ;  // Myz
  source_moment_tensor[5] =  M12 ;  // Mxy 
}

/* 
 * if the nearest point in this thread then search its grid index
 *   return value:
 *      1 - in this thread
 *      0 - not in this thread
 */

  int
src_coord2index(float sx, float sy, float sz,
    int nx, int ny, int nz,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    float *restrict x3d,
    float *restrict y3d,
    float *restrict z3d,
    float *restrict wrk3d,
    int *si, int *sj, int *sk,
    float *sx_inc, float *sy_inc, float *sz_inc)
{
  int is_here = 0; // default outside

  // range of coord of this thread w ghosts
  int x_min = x3d[0];
  int x_max = x3d[0];
  int y_min = y3d[0];
  int y_max = y3d[0];
  int z_min = z3d[0];
  int z_max = z3d[0];
  size_t siz_line  = nx;
  size_t siz_slice = nx * ny;

  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;
        x_min = (x_min < x3d[iptr]) ? x_min : x3d[iptr]; 
        x_max = (x_max > x3d[iptr]) ? x_max : x3d[iptr]; 
        y_min = (y_min < y3d[iptr]) ? y_min : y3d[iptr]; 
        y_max = (y_max > y3d[iptr]) ? y_max : y3d[iptr]; 
        z_min = (z_min < z3d[iptr]) ? z_min : z3d[iptr]; 
        z_max = (z_max > z3d[iptr]) ? z_max : z3d[iptr]; 
      }
    }
  }

  // outside coord range
  if ( sx < x_min || sx > x_max ||
      sy < y_min || sy > y_max ||
      sz < z_min || sz > z_max)
  {
    is_here = 0;
    return is_here;
  }

  // init closest point
  float min_dist = sqrtf(  (sx - x3d[0]) * (sx - x3d[0])
      + (sy - y3d[0]) * (sy - y3d[0])
      + (sz - z3d[0]) * (sz - z3d[0]) );
  int min_dist_i = 0 ;
  int min_dist_j = 0 ;
  int min_dist_k = 0 ;

  // compute distance to each point
  for (int k=0; k<nz; k++) {
    for (int j=0; j<ny; j++) {
      for (int i=0; i<nx; i++)
      {
        size_t iptr = i + j * siz_line + k * siz_slice;

        float x = x3d[iptr];
        float y = y3d[iptr];
        float z = z3d[iptr];

        float DistInt = sqrtf(  (sx - x) * (sx - x)
            + (sy - y) * (sy - y)
            + (sz - z) * (sz - z) );
        wrk3d[iptr] =  DistInt;

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = i;
          min_dist_j = j;
          min_dist_k = k;
        }
      }
    }
  }

  // if nearest index is outside phys region, not here
  if ( min_dist_i < ni1 || min_dist_i > ni2 ||
      min_dist_j < nj1 || min_dist_j > nj2 ||
      min_dist_k < nk1 || min_dist_k > nk2 )
  {
    is_here = 0;
    return is_here;
  }

  // in this thread
  is_here = 1;

  // inverse distance interp curv coord
  //float points_x[3*3*3];
  //float points_y[3*3*3];
  //float points_z[3*3*3];
  //float points_i[3*3*3];
  //float points_j[3*3*3];
  //float points_k[3*3*3];

  //for (int kk=0; kk<3; kk++)
  //{
  //  for (int jj=0; jj<3; jj++)
  //  {
  //    for (int ii=0; ii<3; ii++)
  //    {
  //      int iptr_cube = ii + jj * 3 + kk * 9;
  //      int iptr = (min_dist_i-1+ii) + (min_dist_j-1+jj) * siz_line +
  //                 (min_dist_k-1+kk) * siz_slice;

  //      points_x[iptr_cube] = x3d[iptr];
  //      points_y[iptr_cube] = y3d[iptr];
  //      points_z[iptr_cube] = z3d[iptr];

  //      points_i[iptr_cube] = min_dist_i-1+ii;
  //      points_j[iptr_cube] = min_dist_j-1+jj;
  //      points_k[iptr_cube] = min_dist_k-1+kk;
  //    }
  //  }
  //}

  float points_x[8];
  float points_y[8];
  float points_z[8];
  float points_i[8];
  float points_j[8];
  float points_k[8];

  for (int kk=0; kk<2; kk++)
  {
    for (int jj=0; jj<2; jj++)
    {
      for (int ii=0; ii<2; ii++)
      {
        int cur_i = min_dist_i-1+ii;
        int cur_j = min_dist_j-1+jj;
        int cur_k = min_dist_k-1+kk;

        for (int n3=0; n3<2; n3++) {
          for (int n2=0; n2<2; n2++) {
            for (int n1=0; n1<2; n1++) {
              int iptr_cube = n1 + n2 * 2 + n3 * 4;
              int iptr = (cur_i+n1) + (cur_j+n2) * siz_line +
                (cur_k+n3) * siz_slice;
              points_x[iptr_cube] = x3d[iptr];
              points_y[iptr_cube] = y3d[iptr];
              points_z[iptr_cube] = z3d[iptr];
              points_i[iptr_cube] = cur_i+n1;
              points_j[iptr_cube] = cur_j+n2;
              points_k[iptr_cube] = cur_k+n3;
            }
          }
        }

        if (isPointInHexahedron(sx,sy,sz,points_x,points_y,points_z) == true)
        {
          float si_curv, sj_curv, sk_curv;
          //src_cart2curv_rdinterp(sx,sy,sz,
          //              8,
          //              points_x,points_y,points_z,
          //              points_i,points_j,points_k,
          //              &si_curv, &sj_curv, &sk_curv);

          src_cart2curv_sample(sx,sy,sz,
              8,
              points_x,points_y,points_z,
              points_i,points_j,points_k,
              100,100,100,
              &si_curv, &sj_curv, &sk_curv);

          // convert to return values
          *si = min_dist_i;
          *sj = min_dist_j;
          *sk = min_dist_k;
          *sx_inc = si_curv - min_dist_i;
          *sy_inc = sj_curv - min_dist_j;
          *sz_inc = sk_curv - min_dist_k;

          return is_here;
        }
      }
    }
  }

  // if not in any cube due to bug, set default value
  //    if everything is right, should be return 10 line before
  *si = min_dist_i;
  *sj = min_dist_j;
  *sk = min_dist_k;
  *sx_inc = 0.0;
  *sy_inc = 0.0;
  *sz_inc = 0.0;

  return is_here;
}

/* 
 * interp curv coord using inverse distance interp
 */

  int
src_cart2curv_rdinterp(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_j,
    float *points_k,
    float *si_curv, // interped curv coord
    float *sj_curv,
    float *sk_curv)
{
  float weight[num_points];
  float total_weight = 0.0 ;

  // cal weight
  int at_point_indx = -1;
  for (int i=0; i<num_points; i++)
  {
    float dist = sqrtf ((sx - points_x[i]) * (sx - points_x[i])
        + (sy - points_y[i]) * (sy - points_y[i])
        + (sz - points_z[i]) * (sz - points_z[i])
        );
    if (dist < 1e-9) {
      at_point_indx = i;
    } else {
      weight[i]   = 1.0 / dist;
      total_weight += weight[i];
    }
  }
  // if at a point
  if (at_point_indx > 0) {
    total_weight = 1.0;
    // other weight 0
    for (int i=0; i<num_points; i++) {
      weight[i] = 0.0;
    }
    // point weight 1
    weight[at_point_indx] = 1.0;
  }

  // interp

  *si_curv = 0.0;
  *sj_curv = 0.0;
  *sk_curv = 0.0;

  for (int i=0; i<num_points; i++)
  {
    weight[i] *= 1.0 / total_weight ;

    (*si_curv) += weight[i] * points_i[i];
    (*sj_curv) += weight[i] * points_j[i]; 
    (*sk_curv) += weight[i] * points_k[i];  

    fprintf(stdout,"---- i=%d,weight=%f,points_i=%f,points_j=%f,points_k=%f\n",
        i,weight[i],points_i[i],points_j[i],points_k[i]);
  }

  return 0;
}

/* 
 * find curv coord using sampling
 */

  int
src_cart2curv_sample(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_j,
    float *points_k,
    int    nx_sample,
    int    ny_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sj_curv,
    float *sk_curv)
{
  float Lx[2], Ly[2], Lz[2];

  // init closest point
  float min_dist = sqrtf(  (sx - points_x[0]) * (sx - points_x[0])
      + (sy - points_y[0]) * (sy - points_y[0])
      + (sz - points_z[0]) * (sz - points_z[0]) );
  int min_dist_i = 0 ;
  int min_dist_j = 0 ;
  int min_dist_k = 0 ;

  // linear interp for all sample
  for (int n3=0; n3<nz_sample+1; n3++)
  {
    Lz[1] = (float)(n3) / (float)(nz_sample);
    Lz[0] = 1.0 - Lz[1];
    for (int n2=0; n2<ny_sample+1; n2++)
    {
      Ly[1] = (float)(n2) / (float)(ny_sample);
      Ly[0] = 1.0 - Ly[1];
      for (int n1=0; n1<nx_sample+1; n1++)
      {
        Lx[1] = (float)(n1) / (float)(nx_sample);
        Lx[0] = 1.0 - Lx[1];

        // interp
        float x_pt=0;
        float y_pt=0;
        float z_pt=0;
        for (int kk=0; kk<2; kk++) {
          for (int jj=0; jj<2; jj++) {
            for (int ii=0; ii<2; ii++)
            {
              int iptr_cube = ii + jj * 2 + kk * 4;
              x_pt += Lx[ii]*Ly[jj]*Lz[kk] * points_x[iptr_cube];
              y_pt += Lx[ii]*Ly[jj]*Lz[kk] * points_y[iptr_cube];
              z_pt += Lx[ii]*Ly[jj]*Lz[kk] * points_z[iptr_cube];
            }
          }
        }

        // find min dist
        float DistInt = sqrtf(  (sx - x_pt) * (sx - x_pt)
            + (sy - y_pt) * (sy - y_pt)
            + (sz - z_pt) * (sz - z_pt) );

        // replace closest point
        if (min_dist > DistInt)
        {
          min_dist = DistInt;
          min_dist_i = n1;
          min_dist_j = n2;
          min_dist_k = n3;
        }
      } // n1
    } // n2
  } // n3

  *si_curv = points_i[0] + (float)min_dist_i / (float)nx_sample;
  *sj_curv = points_j[0] + (float)min_dist_j / (float)ny_sample;
  *sk_curv = points_k[0] + (float)min_dist_k / (float)nz_sample;

  return 0;
}

