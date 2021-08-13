/*
 * source term related processing
 */

// todo:

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>

#include "fdlib_math.h"
#include "interp.h"
#include "fdlib_mem.h"
#include "src_t.h"
#include "isPointInHexahedron.h"

/*
 * src_t alloc
 */

int
src_init(src_t *src, int force_actived, int moment_actived,
         int num_of_src, int max_nt, int max_stage, int max_ext)
{
  // set default value
  src->total_number = num_of_src;
  src->max_nt    = max_nt;
  src->max_stage = max_stage;
  src->max_ext   = max_ext;

  src->force_actived   = force_actived;
  src->moment_actived   = moment_actived;

  // allocate var
  src->si = (int *)malloc(num_of_src*sizeof(int));
  src->sj = (int *)malloc(num_of_src*sizeof(int));
  src->sk = (int *)malloc(num_of_src*sizeof(int));
  src->it_begin = (int *)malloc(num_of_src*sizeof(int));
  src->it_end   = (int *)malloc(num_of_src*sizeof(int));
  src->ext_num  = (int *)malloc(num_of_src*sizeof(int));
  src->ext_indx = (int *)malloc(num_of_src*max_ext * sizeof(int  ));
  src->ext_coef = (float *)malloc(num_of_src*max_ext * sizeof(float));

  src->Fx = NULL;
  src->Fy = NULL;
  src->Fz = NULL;

  if (force_actived == 1) {
    src->Fx = (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Fy = (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Fz = (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    for (int iptr=0; iptr < max_stage * max_nt * num_of_src; iptr++) {
      src->Fx[iptr] = 0.0;
      src->Fy[iptr] = 0.0;
      src->Fz[iptr] = 0.0;
    }
  }

  src->Mxx = NULL;
  src->Myy = NULL;
  src->Mzz = NULL;
  src->Mxz = NULL;
  src->Myz = NULL;
  src->Mxy = NULL;

  if (moment_actived == 1) {
    src->Mxx= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Myy= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mzz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mxz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Myz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mxy= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    for (int iptr=0; iptr < max_stage * max_nt * num_of_src; iptr++) {
      src->Mxx[iptr] = 0.0;
      src->Myy[iptr] = 0.0;
      src->Mzz[iptr] = 0.0;
      src->Mxz[iptr] = 0.0;
      src->Myz[iptr] = 0.0;
      src->Mxy[iptr] = 0.0;
    }
  }

  return 0;
}

int
src_set_time(src_t *src, int it, int istage)
{
  src->it     = it;
  src->istage = istage;

  return 0;
}

/*
 * convert coord to global index using MPI
 */

int
src_coord_to_glob_indx(gdinfo_t *gdinfo,
                       gdcurv_t *gdcurv,
                       float sx,
                       float sy,
                       float sz,
                       MPI_Comm comm,
                       int myid,
                       int   *ou_si, int *ou_sj, int *ou_sk,
                       float *ou_sx_inc, float *ou_sy_inc, float *ou_sz_inc,
                       float *restrict wrk3d)
{
  int is_here = 0;

  int si_glob = 0;
  int sj_glob = 0;
  int sk_glob = 0;
  float sx_inc = 0.0;
  float sy_inc = 0.0;
  float sz_inc = 0.0;
  int si = 0;
  int sj = 0;
  int sk = 0;

  // if located in this thread
  is_here = src_coord_to_local_indx(gdinfo,gdcurv, sx,sy,sz,
                                    &si, &sj, &sk, &sx_inc, &sy_inc, &sz_inc,
                                    wrk3d);

  // if in this thread
  if ( is_here == 1)
  {
    // conver to global index
    si_glob = gd_info_ind_lcext2glphy_i(si, gdinfo);
    sj_glob = gd_info_ind_lcext2glphy_j(sj, gdinfo);
    sk_glob = gd_info_ind_lcext2glphy_k(sk, gdinfo);
    fprintf(stdout," -- located to local index = %d %d %d\n", si,sj,sk);
    fprintf(stdout," -- located to global index = %d %d %d\n", 
                          si_glob, sj_glob, sk_glob);
    fprintf(stdout," --  with shift = %f %f %f\n", sx_inc,sy_inc,sz_inc);
  } else {
    //fprintf(stdout," -- not in this thread %d\n", myid);
  }

  // reduce global index and shift values
  int sendbufi = si_glob;
  MPI_Allreduce(&sendbufi, &si_glob, 1, MPI_INT, MPI_MAX, comm);

  sendbufi = sj_glob;
  MPI_Allreduce(&sendbufi, &sj_glob, 1, MPI_INT, MPI_MAX, comm);

  sendbufi = sk_glob;
  MPI_Allreduce(&sendbufi, &sk_glob, 1, MPI_INT, MPI_MAX, comm);

  float sendbuf = sx_inc;
  MPI_Allreduce(&sendbuf, &sx_inc, 1, MPI_INT, MPI_SUM, comm);

  sendbuf = sy_inc;
  MPI_Allreduce(&sendbuf, &sy_inc, 1, MPI_INT, MPI_SUM, comm);

  sendbuf = sz_inc;
  MPI_Allreduce(&sendbuf, &sz_inc, 1, MPI_INT, MPI_SUM, comm);

  fprintf(stdout," --myid=%d,index=%d %d %d,shift = %f %f %f\n",
      myid,si_glob,sj_glob,sk_glob, sx_inc,sy_inc,sz_inc);

  *ou_si = si_glob;
  *ou_sj = sj_glob;
  *ou_sk = sk_glob;
  *ou_sx_inc = sx_inc;
  *ou_sy_inc = sy_inc;
  *ou_sz_inc = sz_inc;

  return is_here; 
}

/*
 * if extend global index in this thread
 */

int
src_glob_ext_ishere(int si, int sj, int sk, int half_ext, gdinfo_t *gdinfo)
{
  int is_here = 0;

  if (si-half_ext <= gdinfo->ni2_to_glob_phys0 && // exted left point is less than right bdry
      si+half_ext >= gdinfo->ni1_to_glob_phys0 && // exted right point is larger than left bdry
      sj-half_ext <= gdinfo->nj2_to_glob_phys0 && 
      sj+half_ext >= gdinfo->nj1_to_glob_phys0 &&
      sk-half_ext <= gdinfo->nk2_to_glob_phys0 && 
      sk+half_ext >= gdinfo->nk1_to_glob_phys0)
  {
    is_here = 1;
  }

  return is_here;
}

/*
 * for source info read from .json file
 */

int
src_set_by_par(gdinfo_t *gdinfo,
               gdcurv_t *gdcurv,
               src_t    *src,
               float t0,
               float dt,
               int   max_stage,
               float *rk_stage_time,
               int   npoint_half_ext,
               char  *in_source_name,
               int   in_num_of_src,
               int   **source_index,
               float **source_inc,
               float **source_coords,
               float **force_vector, 
               int   *source_force_actived,
               float **moment_tensor,
               int   *source_moment_actived,
               char  **wavelet_name,
               float **wavelet_coefs,
               float *wavelet_tstart,
               float *wavelet_tend,
               MPI_Comm comm, 
               int myid,
               int verbose)
{
  int ierr = 0;

  // get grid info from gdinfo
  int   ni1 = gdinfo->ni1;
  int   ni2 = gdinfo->ni2;
  int   nj1 = gdinfo->nj1;
  int   nj2 = gdinfo->nj2;
  int   nk1 = gdinfo->nk1;
  int   nk2 = gdinfo->nk2;
  int   nx  = gdinfo->nx ;
  int   ny  = gdinfo->ny ;
  int   nz  = gdinfo->nz ;
  int   npoint_ghosts = gdinfo->npoint_ghosts;
  size_t siz_line = gdinfo->siz_iy;
  size_t siz_slice= gdinfo->siz_iz;

  // get total elem of exted src region for a single point
  //    int max_ext = 7 * 7 * 7;
  int len_ext = 2*npoint_half_ext+1;
  int max_ext = len_ext * len_ext * len_ext;

  // local
  int si,sj,sk;
  int si_glob,sj_glob,sk_glob;
  float sx_inc, sy_inc, sz_inc;

  // workspace 3d var for distance calculation, only used for coord input
  float *wrk3d=NULL;
  wrk3d = (float *) fdlib_mem_calloc_1d_float(nx*ny*nz,0.0,"src_set_by_par");

  // set evtnm
  sprintf(src->evtnm,"%s",in_source_name);

  //
  // first run: loop all src to get info
  //
  int max_nt = 0;
  int num_of_src_here = 0;
  int force_actived = 0;
  int moment_actived = 0;

  for (int is=0; is < in_num_of_src; is++)
  {
    // get max_nt
    int  it_begin = (int) (wavelet_tstart[is] / dt);
    int  it_end   = (int) ((wavelet_tend[is] / dt + 0.5));
    int  nt_total_wavelet = it_end - it_begin + 1;
    max_nt = max_nt > nt_total_wavelet ? max_nt : nt_total_wavelet; 

    // check if force and moment used
    if (source_force_actived[is] == 1) force_actived = 1;
    if (source_moment_actived[is] == 1) moment_actived = 1;

    // count num of src in this thread

    // convert coord to glob index
    if (source_index[is][0] < 0)
    {
      float sx = source_coords[is][0];
      float sy = source_coords[is][1];
      float sz = source_coords[is][2];

      fprintf(stdout,"locate source by coord (%f,%f,%f) ...\n",sx,sy,sz);
      fflush(stdout);
      src_coord_to_glob_indx(gdinfo,gdcurv,sx,sy,sz,comm,myid,
                             &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc,
                             wrk3d);
      // keep index to avoid duplicat run
      source_index[is][0] = si_glob;
      source_index[is][1] = sj_glob;
      source_index[is][2] = sk_glob;
      source_inc[is][0] = sx_inc;
      source_inc[is][1] = sy_inc;
      source_inc[is][2] = sz_inc;
    } else {
      si_glob = source_index[is][0];
      sj_glob = source_index[is][1];
      sk_glob = source_index[is][2];
      source_inc[is][0] = 0.0;
      source_inc[is][1] = 0.0;
      source_inc[is][2] = 0.0;
    }
    // check if in this thread using index
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
    {
      num_of_src_here += 1;
    }
  }

  // alloc src_t
  src_init(src,force_actived,moment_actived,num_of_src_here,max_nt,max_stage,max_ext);

  //
  // second run to set each src in this thread
  //
  int is_local = 0;
  for (int is=0; is < in_num_of_src; is++)
  {
    si_glob = source_index[is][0];
    sj_glob = source_index[is][1];
    sk_glob = source_index[is][2];

    // check if in this thread
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
    {
      // convert to local index
      si = gd_info_ind_glphy2lcext_i(si_glob, gdinfo);
      sj = gd_info_ind_glphy2lcext_j(sj_glob, gdinfo);
      sk = gd_info_ind_glphy2lcext_k(sk_glob, gdinfo);
      // keep
      src->si[is_local] = si;
      src->sj[is_local] = sj;
      src->sk[is_local] = sk;

      // time step, considering t0
      int  it_begin = (int) ( (wavelet_tstart[is] - t0) / dt);
      int  it_end   = (int) ( ((wavelet_tend[is] - t0) / dt) + 0.5);

      src->it_begin[is_local] = it_begin;
      src->it_end  [is_local] = it_end  ;

      // for wavelet
      for (int it=it_begin; it<=it_end; it++)
      {
        int it_to_it1 = (it - it_begin);
        float t_shift = wavelet_tstart[is] - (it_begin * dt + t0);

        // order as: istage, is, it for better localization
        //    not work because it_begin may diff for diff source
        // int iptr_it = is_local * max_stage + it_to_it1 * max_stage * num_of_src_here;
        // use istage, it, is order
        int iptr_it = is_local * max_nt * max_stage + it_to_it1 * max_stage;

        for (int istage=0; istage<max_stage; istage++)
        {
          // time relative to start time of this source, considering diff from int conversion
          float t = it_to_it1 * dt + rk_stage_time[istage] * dt - t_shift;

          float stf_val;
          if (strcmp(wavelet_name[is], "ricker")==0) {
            stf_val = fun_ricker(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else if (strcmp(wavelet_name[is], "gaussian")==0) {
            stf_val = fun_gauss(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else if (strcmp(wavelet_name[is], "ricker_deriv")==0) {
            stf_val = fun_ricker_deriv(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else if (strcmp(wavelet_name[is], "gaussian_deriv")==0) {
            stf_val = fun_gauss_deriv(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else{
            fprintf(stderr,"wavelet_name=%s\n", wavelet_name[is]); 
            fprintf(stderr,"   not implemented yet\n"); 
            fflush(stderr);
          }

          int iptr = iptr_it + istage;

          if (source_force_actived[is]==1) {
            src->Fx[iptr]  = stf_val * force_vector[is][0];
            src->Fy[iptr]  = stf_val * force_vector[is][1];
            src->Fz[iptr]  = stf_val * force_vector[is][2];
          }

          if (source_moment_actived[is]==1) {
            src->Mxx[iptr] = stf_val * moment_tensor[is][0];
            src->Myy[iptr] = stf_val * moment_tensor[is][1];
            src->Mzz[iptr] = stf_val * moment_tensor[is][2];
            src->Mxz[iptr] = stf_val * moment_tensor[is][3];
            src->Myz[iptr] = stf_val * moment_tensor[is][4];
            src->Mxy[iptr] = stf_val * moment_tensor[is][5];
          }
        } // istage
      } // it

      // for extended points and coefs
      sx_inc = source_inc[is][0];
      sy_inc = source_inc[is][1];
      sz_inc = source_inc[is][2];
      float wid_gauss = npoint_half_ext / 2.0;
      float *this_ext_coef = src->ext_coef + is_local * max_ext;
      src_cal_norm_delt3d(this_ext_coef, sx_inc, sy_inc, sz_inc,
                          wid_gauss, wid_gauss, wid_gauss, npoint_half_ext);

      size_t iptr_ext = 0;
      for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
      {
        for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
        {
          for (int i=si-npoint_half_ext; i<=si+npoint_half_ext; i++)
          {
            if (gd_info_lindx_is_inner(i,j,k,gdinfo)==1)
            {
              // Note index need match coef
              int iptr_grid = i + j * siz_line + k * siz_slice;
              int iptr_coef =  (i-(si-npoint_half_ext))
                              + len_ext * (j-(sj-npoint_half_ext)) 
                              + len_ext * len_ext *(k-(sk-npoint_half_ext));
              src->ext_indx[iptr_ext + is_local * max_ext] = iptr_grid;
              src->ext_coef[iptr_ext + is_local * max_ext] = this_ext_coef[iptr_coef];
              iptr_ext++;
            }
          }
        }
      }
      // only count index inside phys region for this thread
      src->ext_num[is_local] = iptr_ext;

      is_local += 1;
    }
  }

  // free working space
  free(wrk3d);

  return ierr;
}


int
src_read_locate_valsrc(gdinfo_t *gdinfo,
                       gdcurv_t *gdcurv,
                       src_t    *src,
                       char *pfilepath,
                       float t0,
                       float dt,
                       int   max_stage,
                       float *rk_stage_time,
                       int   npoint_half_ext,
                       MPI_Comm comm,
                       int myid,
                       int verbose)
{
  int ierr = 0;

  FILE *fp =NULL;
  char str[500];
  char *in_source_name = (char *)malloc(500*sizeof(char));
  int in_num_force;
  int in_num_moment;
  int in_num_source;
  int nt_in;       // numbers_time_steps from inputfile
  float dt_in;     // time_step from inputfile
  // read sample source value from inputfile
  if ((fp = fopen(pfilepath, "r"))==NULL) fprintf(stdout,"fail to open");
  fgets(str,500,fp);
  sscanf(str,"%s",in_source_name);
  fgets(str,500,fp);
  sscanf(str,"%d %d",&in_num_force,&in_num_moment);
  fgets(str,500,fp);
  sscanf(str,"%f %d",&dt_in,&nt_in);

  // set evtnm
  sprintf(src->evtnm,"%s",in_source_name);
  in_num_source = in_num_force + in_num_moment;
  float **source_coords = NULL;

  if(in_num_source > 0)
  {
    source_coords = (float **) fdlib_mem_malloc_2l_float(in_num_source,CONST_NDIM,"source_coords");
  }

  // first read coords to determine src whether in this thread
  for (int i=0; i<in_num_source; i++)
  {
    fgets(str,500,fp);
    sscanf(str,"%f %f %f",&source_coords[i][0],&source_coords[i][1],&source_coords[i][2]);
  }

  // get grid info from gdinfo
  int   ni1 = gdinfo->ni1;
  int   ni2 = gdinfo->ni2;
  int   nj1 = gdinfo->nj1;
  int   nj2 = gdinfo->nj2;
  int   nk1 = gdinfo->nk1;
  int   nk2 = gdinfo->nk2;
  int   nx  = gdinfo->nx ;
  int   ny  = gdinfo->ny ;
  int   nz  = gdinfo->nz ;
  int   npoint_ghosts = gdinfo->npoint_ghosts;
  size_t siz_line = gdinfo->siz_iy;
  size_t siz_slice= gdinfo->siz_iz;

  // get total elem of exted src region for a single point
  //    int max_ext = 7 * 7 * 7;
  int len_ext = 2*npoint_half_ext+1;
  int max_ext = len_ext * len_ext * len_ext;
  // local
  int si,sj,sk;
  int si_glob,sj_glob,sk_glob;
  float sx_inc, sy_inc, sz_inc;

  // workspace 3d var for distance calculation, only used for coord input
  float *wrk3d=NULL;
  wrk3d = (float *) fdlib_mem_calloc_1d_float(nx*ny*nz,0.0,"src_set_by_par");

  //
  // first run: loop all src to get info
  //
  int num_of_src_here = 0;
  int force_actived = 0;
  int moment_actived = 0;
  int   max_nt = 0;
  int *source_in_thread = fdlib_mem_calloc_1d_int(in_num_source,-1, "source_in_this_thread");
  int **source_index = (int **) fdlib_mem_malloc_2l_int(in_num_source,CONST_NDIM,"source_index");
  float **source_inc = (float **) fdlib_mem_malloc_2l_float(in_num_source,CONST_NDIM,"source_inc");
  max_nt = (int) (((nt_in-1)*dt_in / dt)+ 0.5) + 1; 


  for (int is=0; is < in_num_source; is++)
  {
    // count num of src in this thread
    
    // convert coord to glob index
      float sx = source_coords[is][0];
      float sy = source_coords[is][1];
      float sz = source_coords[is][2];
      fprintf(stdout,"locate source by coord (%f,%f,%f) ...\n",sx,sy,sz);
      fflush(stdout);
      src_coord_to_glob_indx(gdinfo,gdcurv,sx,sy,sz,comm,myid,
                             &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc,
                             wrk3d);
      // keep index to avoid duplicat run
      source_index[is][0] = si_glob;
      source_index[is][1] = sj_glob;
      source_index[is][2] = sk_glob;
      source_inc[is][0] = sx_inc;
      source_inc[is][1] = sy_inc;
      source_inc[is][2] = sz_inc;
    // check if in this thread using index
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
    {
      num_of_src_here += 1;
      source_in_thread[is] = 1;
    }
  }

  // check if force and moment used
  if (in_num_force >= 1 && num_of_src_here > 0) force_actived = 1;
  if (in_num_moment >= 1 && num_of_src_here > 0) moment_actived = 1;

  float *wavelet_tstart = NULL;
  float ***force_vector = NULL;
  float ***moment_tensor = NULL;
  char  **moment_wavelet_mechism = NULL;
  if(in_num_source>0)
  {
    wavelet_tstart = (float *) fdlib_mem_calloc_1d_float(in_num_source, 0.0, "wavelet_tstart");
  }
  if(in_num_force>0)
  {
    force_vector = (float ***) fdlib_mem_calloc_3l_float(in_num_force, CONST_NDIM, nt_in, 0.0, "force_vector");
  }
  if(in_num_moment > 0)
  {
    moment_wavelet_mechism = fdlib_mem_malloc_2l_char(in_num_moment, 500, "moment_wavelet_mechism");
    moment_tensor = (float ***) fdlib_mem_calloc_3l_float(in_num_moment, CONST_NDIM_2, nt_in, 0.0, "moment_tensor");
  }

  //if (num_of_src_here>0)
  if (in_num_source>0)
  {
    // read force
    for (int is=0; is<in_num_force; is++)
    {	
      fgets(str,500,fp);
      sscanf(str,"%f",&wavelet_tstart[is]);
      if (source_in_thread[is] == 1)
      {  
        for(int k=0; k<nt_in; k++)
        {
          fgets(str,500,fp);
          sscanf(str,"%f %f %f",&force_vector[is][0][k],&force_vector[is][1][k],&force_vector[is][2][k]);
        }
      }

      if(source_in_thread[is] == -1)
      {
        for(int k=0; k<nt_in; k++)
        {
          fgets(str,500,fp);
        }
      }  
    }

    //read moment
    for(int is=0 + in_num_force; is<in_num_source; is++)
    {
      int index = is - in_num_force;
      fgets(str,500,fp);
      sscanf(str,"%f",&wavelet_tstart[is]);
      fgets(str,500,fp);
      sscanf(str,"%s",moment_wavelet_mechism[index]);
      if (source_in_thread[is] == 1)   //moment source is in this thread
      {

        if (strcmp("mechanism_angle",moment_wavelet_mechism[index]) == 0)
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
            moment_tensor[is-in_num_force][0][k] = M0*temp_moment[0];
            moment_tensor[is-in_num_force][1][k] = M0*temp_moment[1];
            moment_tensor[is-in_num_force][2][k] = M0*temp_moment[2];
            moment_tensor[is-in_num_force][3][k] = M0*temp_moment[3];
            moment_tensor[is-in_num_force][4][k] = M0*temp_moment[4];
            moment_tensor[is-in_num_force][5][k] = M0*temp_moment[5];
          }
        }
        if (strcmp("moment_tensor",moment_wavelet_mechism[index]) == 0)
        {
          for (int k=0;k<nt_in;k++)
          {
            fgets(str,500,fp);
            sscanf(str,"%f %f %f %f %f %f",&moment_tensor[is-in_num_force][0][k], &moment_tensor[is-in_num_force][1][k], 
                                           &moment_tensor[is-in_num_force][2][k], &moment_tensor[is-in_num_force][3][k],
                                           &moment_tensor[is-in_num_force][4][k], &moment_tensor[is-in_num_force][5][k]);
          }
        }
      }

      if (source_in_thread[is] == -1)
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
  //
  // second run to set each src in this thread
  //

  // alloc src_t
  src_init(src,force_actived,moment_actived,num_of_src_here,max_nt,max_stage,max_ext);

  float *t_in = (float *)malloc(nt_in*sizeof(float));
  int is_local = 0;
  for (int is=0; is < in_num_source; is++)
  {
    si_glob = source_index[is][0];
    sj_glob = source_index[is][1];
    sk_glob = source_index[is][2];

    // check if in this thread
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
    {
      // convert to local index
      si = gd_info_ind_glphy2lcext_i(si_glob, gdinfo);
      sj = gd_info_ind_glphy2lcext_j(sj_glob, gdinfo);
      sk = gd_info_ind_glphy2lcext_k(sk_glob, gdinfo);
      // keep
      src->si[is_local] = si;
      src->sj[is_local] = sj;
      src->sk[is_local] = sk;

      // time step, considering t0
      int  it_begin = (int) ( (wavelet_tstart[is] - t0) / dt);
      int  it_end   = it_begin + max_nt - 1;

      src->it_begin[is_local] = it_begin;
      src->it_end  [is_local] = it_end  ;
      // interp1d order
      int order = 3;     
      // import val source time, interp need used 
      for(int it_in=0; it_in<nt_in; it_in++)
      {
        t_in[it_in] = wavelet_tstart[is] + it_in*dt_in;
      }
      // for wavelet
      for (int it=it_begin; it<=it_end; it++)
      {
        int it_to_it1 = (it - it_begin);
        float t_shift = wavelet_tstart[is] - (it_begin * dt + t0);

        // order as: istage, is, it for better localization
        //    not work because it_begin may diff for diff source
        // int iptr_it = is_local * max_stage + it_to_it1 * max_stage * num_of_src_here;
        // use istage, it, is order
        int iptr_it = is_local * max_nt * max_stage + it_to_it1 * max_stage;

        for (int istage=0; istage<max_stage; istage++)
        {
          // time relative to start time of this source, considering diff from int conversion
          float t = it * dt + rk_stage_time[istage] * dt - t_shift;
          if ( is < in_num_force) {
          float this_Fx, this_Fy, this_Fz;
          // interp1d give t to get Force vector
          this_Fx = LagInterp_Piecewise_1d(t_in, &force_vector[is][0][0], nt_in, order, 
                                                               wavelet_tstart[is], dt_in, t);
          this_Fy = LagInterp_Piecewise_1d(t_in, &force_vector[is][1][0], nt_in, order, 
                                                               wavelet_tstart[is], dt_in, t);
          this_Fz = LagInterp_Piecewise_1d(t_in, &force_vector[is][2][0], nt_in, order, 
                                                               wavelet_tstart[is], dt_in, t);

          int iptr = iptr_it + istage;

          src->Fx[iptr]  = this_Fx;
          src->Fy[iptr]  = this_Fy;
          src->Fz[iptr]  = this_Fz;
          }

          if ( is >= in_num_force) {
          float this_Mxx, this_Myy, this_Mzz, this_Mxz, this_Myz, this_Mxy;
          // interp1d give t to get moment tensor
          this_Mxx = LagInterp_Piecewise_1d(t_in, &moment_tensor[is-in_num_force][0][0], nt_in, order,
                                                         wavelet_tstart[is], dt_in, t);
          this_Myy = LagInterp_Piecewise_1d(t_in, &moment_tensor[is-in_num_force][1][0], nt_in, order,
                                                         wavelet_tstart[is], dt_in, t);
          this_Mzz = LagInterp_Piecewise_1d(t_in, &moment_tensor[is-in_num_force][2][0], nt_in, order,
                                                         wavelet_tstart[is], dt_in, t);
          this_Mxz = LagInterp_Piecewise_1d(t_in, &moment_tensor[is-in_num_force][3][0], nt_in, order,
                                                         wavelet_tstart[is], dt_in, t);
          this_Myz = LagInterp_Piecewise_1d(t_in, &moment_tensor[is-in_num_force][4][0], nt_in, order,
                                                         wavelet_tstart[is], dt_in, t);
          this_Mxy = LagInterp_Piecewise_1d(t_in, &moment_tensor[is-in_num_force][5][0], nt_in, order,
                                                         wavelet_tstart[is], dt_in, t);

          int iptr = iptr_it + istage;
          src->Mxx[iptr] = this_Mxx;
          src->Myy[iptr] = this_Myy;
          src->Mzz[iptr] = this_Mzz;
          src->Mxz[iptr] = this_Mxz;
          src->Myz[iptr] = this_Myz;
          src->Mxy[iptr] = this_Mxy;
          }
        } // istage
      } // it

      // for extended points and coefs
      sx_inc = source_inc[is][0];
      sy_inc = source_inc[is][1];
      sz_inc = source_inc[is][2];
      float wid_gauss = npoint_half_ext / 2.0;
      float *this_ext_coef = src->ext_coef + is_local * max_ext;
      src_cal_norm_delt3d(this_ext_coef, sx_inc, sy_inc, sz_inc,
                          wid_gauss, wid_gauss, wid_gauss, npoint_half_ext);

      size_t iptr_ext = 0;
      for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
      {
        for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
        {
          for (int i=si-npoint_half_ext; i<=si+npoint_half_ext; i++)
          {
            if (gd_info_lindx_is_inner(i,j,k,gdinfo)==1)
            {
              // Note index need match coef
              int iptr_grid = i + j * siz_line + k * siz_slice;
              int iptr_coef =  (i-(si-npoint_half_ext))
                              + len_ext * (j-(sj-npoint_half_ext)) 
                              + len_ext * len_ext *(k-(sk-npoint_half_ext));
              src->ext_indx[iptr_ext + is_local * max_ext] = iptr_grid;
              src->ext_coef[iptr_ext + is_local * max_ext] = this_ext_coef[iptr_coef];
              iptr_ext++;
            }
          }
        }
      }
      // only count index inside phys region for this thread
      src->ext_num[is_local] = iptr_ext;

      is_local += 1;
    }
  }
  // free
  if(in_num_source > 0)
  {
    free(wrk3d);
    free(in_source_name);
    free(source_in_thread);
    free(wavelet_tstart);
    fdlib_mem_free_2l_float(source_coords, in_num_source, "free_source_coords");
    fdlib_mem_free_2l_float(source_inc, in_num_source, "free_source_inc");
    fdlib_mem_free_2l_int(source_index, in_num_source, "free_source_index");
  }
  if(in_num_force > 0)
  {
    fdlib_mem_free_3l_float(force_vector, in_num_force, CONST_NDIM, "free_force_vector");
  }
  if(in_num_moment > 0)
  {
    fdlib_mem_free_3l_float(moment_tensor, in_num_moment, CONST_NDIM_2, "free_moment_tensor");
    fdlib_mem_free_2l_char(moment_wavelet_mechism,in_num_moment,"free_wavelet_name");
  }
  return ierr;
}



int
src_read_locate_anasrc(gdinfo_t *gdinfo,
                       gdcurv_t *gdcurv,
                       src_t    *src,
                       char *pfilepath,
                       float t0,
                       float dt,
                       int   max_stage,
                       float *rk_stage_time,
                       int   npoint_half_ext,
                       MPI_Comm comm,
                       int myid,
                       int verbose)
{
  int ierr = 0;

  FILE* fp =NULL;
  char str[500];
  int in_num_force;
  int in_num_moment;
  int in_num_source;
  float wavelet_tlen;
  char *in_source_name = (char *)malloc(500 * sizeof(char));
  // read analysis source from inputfile
  if ((fp = fopen(pfilepath, "r"))==NULL) fprintf(stdout,"fail to open");
  fgets(str,500,fp);
  sscanf(str,"%s",in_source_name);
  fgets(str,500,fp);
  sscanf(str,"%d %d",&in_num_force,&in_num_moment);
  fgets(str,500,fp);
  sscanf(str,"%f",&wavelet_tlen);
  // set evtnm
  sprintf(src->evtnm,"%s",in_source_name);
  in_num_source = in_num_force + in_num_moment;
  float **source_coords = NULL;
  float **wavelet_coefs = NULL;
  float *wavelet_tstart = NULL;
  char  **wavelet_name = NULL;
  float **force_vector = NULL;
  float **moment_tensor = NULL;
  char  **moment_wavelet_mechism = NULL;
  int **source_index = NULL;
  float **source_inc = NULL;
  if(in_num_source > 0)
  {
    source_coords = (float **) fdlib_mem_malloc_2l_float(in_num_source,CONST_NDIM,"source_coords");
    wavelet_coefs = (float **) fdlib_mem_malloc_2l_float(in_num_source,2,"wavelet_coefs");
    wavelet_tstart = (float *) fdlib_mem_calloc_1d_float(in_num_source,0.0,"wavelet_tstart");
    wavelet_name = (char **) fdlib_mem_malloc_2l_char(in_num_source,100,"source_time_function type");
    source_index = (int **) fdlib_mem_malloc_2l_int(in_num_source,CONST_NDIM,"source_index");
    source_inc = (float **) fdlib_mem_malloc_2l_float(in_num_source,CONST_NDIM,"source_inc");
  }
  if(in_num_force > 0)
  {
    force_vector = (float **) fdlib_mem_malloc_2l_float(in_num_force,CONST_NDIM,"force_vector");
  }
  if(in_num_moment > 0)
  {
    moment_tensor = (float **) fdlib_mem_malloc_2l_float(in_num_moment,CONST_NDIM_2,"moment_tensor");
    moment_wavelet_mechism = (char **) fdlib_mem_malloc_2l_char(in_num_moment,100,"mechanism type:tensor or angle");
  }
  // read force
  if (in_num_force > 0)
  {
    for (int i=0;i<in_num_force;i++)
    {
      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&source_coords[i][0],&source_coords[i][1],&source_coords[i][2]);

      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&force_vector[i][0], &force_vector[i][1], &force_vector[i][2]);

      fgets(str,500,fp);
      sscanf(str,"%s",wavelet_name[i]);

      fgets(str,500,fp);
      sscanf(str,"%f %f",&wavelet_coefs[i][0], &wavelet_coefs[i][1]);

      fgets(str,500,fp);
      sscanf(str,"%f",&wavelet_tstart[i]);
    }
  }
  // read moment
  if (in_num_moment > 0)
  {
    for(int i=0;i<in_num_moment;i++)
    {
      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&source_coords[i+in_num_force][0],&source_coords[i+in_num_force][1],
                            &source_coords[i+in_num_force][2]);

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
          moment_tensor[i][j]=M0*temp_moment[j];
          //fprintf(stdout,"moment_tensor is %f \n",moment_tensor[i][j]);
        }
      } 
      if (strcmp("moment_tensor",moment_wavelet_mechism[i])==0)
      {  
        fgets(str,500,fp);
        sscanf(str,"%f %f %f %f %f %f",&moment_tensor[i][0],&moment_tensor[i][1],
            &moment_tensor[i][2],&moment_tensor[i][3],&moment_tensor[i][4],&moment_tensor[i][5]);
      } 
      fgets(str,500,fp);
      sscanf(str,"%s",wavelet_name[i+in_num_force]);

      fgets(str,500,fp);
      sscanf(str,"%f %f",&wavelet_coefs[i+in_num_force][0], &wavelet_coefs[i+in_num_force][1]);

      fgets(str,500,fp);
      sscanf(str,"%f",&wavelet_tstart[i+in_num_force]);
    }
  }
  fclose(fp);

  // get grid info from gdinfo
  int   ni1 = gdinfo->ni1;
  int   ni2 = gdinfo->ni2;
  int   nj1 = gdinfo->nj1;
  int   nj2 = gdinfo->nj2;
  int   nk1 = gdinfo->nk1;
  int   nk2 = gdinfo->nk2;
  int   nx  = gdinfo->nx ;
  int   ny  = gdinfo->ny ;
  int   nz  = gdinfo->nz ;
  int   npoint_ghosts = gdinfo->npoint_ghosts;
  size_t siz_line = gdinfo->siz_iy;
  size_t siz_slice= gdinfo->siz_iz;

  // get total elem of exted src region for a single point
  //    int max_ext = 7 * 7 * 7;
  int len_ext = 2*npoint_half_ext+1;
  int max_ext = len_ext * len_ext * len_ext;
  // local
  int si,sj,sk;
  int si_glob,sj_glob,sk_glob;
  float sx_inc, sy_inc, sz_inc;

  // workspace 3d var for distance calculation, only used for coord input
  float *wrk3d=NULL;
  wrk3d = (float *) fdlib_mem_calloc_1d_float(nx*ny*nz,0.0,"src_set_by_par");

  //
  // first run: loop all src to get info
  //
  int num_of_src_here = 0;
  int force_actived = 0;
  int moment_actived = 0;
  int   max_nt = 0;
  max_nt = (int) (wavelet_tlen / dt + 0.5) + 1; 
  // check if force and moment used
  if (in_num_force >= 1) force_actived = 1;
  if (in_num_moment >= 1) moment_actived = 1;

  for (int is=0; is < in_num_source; is++)
  {
    // count num of src in this thread
    
    // convert coord to glob index
      float sx = source_coords[is][0];
      float sy = source_coords[is][1];
      float sz = source_coords[is][2];
      fprintf(stdout,"locate source by coord (%f,%f,%f) ...\n",sx,sy,sz);
      fflush(stdout);
      src_coord_to_glob_indx(gdinfo,gdcurv,sx,sy,sz,comm,myid,
                             &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc,
                             wrk3d);
      // keep index to avoid duplicat run
      source_index[is][0] = si_glob;
      source_index[is][1] = sj_glob;
      source_index[is][2] = sk_glob;
      source_inc[is][0] = sx_inc;
      source_inc[is][1] = sy_inc;
      source_inc[is][2] = sz_inc;
    // check if in this thread using index
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
    {
      num_of_src_here += 1;
    }
  }

  // alloc src_t
  src_init(src,force_actived,moment_actived,num_of_src_here,max_nt,max_stage,max_ext);

  //
  // second run to set each src in this thread
  //
  int is_local = 0;
  for (int is=0; is < in_num_source; is++)
  {
    si_glob = source_index[is][0];
    sj_glob = source_index[is][1];
    sk_glob = source_index[is][2];

    // check if in this thread
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
    {
      // convert to local index
      si = gd_info_ind_glphy2lcext_i(si_glob, gdinfo);
      sj = gd_info_ind_glphy2lcext_j(sj_glob, gdinfo);
      sk = gd_info_ind_glphy2lcext_k(sk_glob, gdinfo);
      // keep
      src->si[is_local] = si;
      src->sj[is_local] = sj;
      src->sk[is_local] = sk;

      // time step, considering t0
      int  it_begin = (int) ( (wavelet_tstart[is] - t0) / dt);
      int  it_end   = it_begin + max_nt - 1;

      src->it_begin[is_local] = it_begin;
      src->it_end  [is_local] = it_end  ;

      // for wavelet
      for (int it=it_begin; it<=it_end; it++)
      {
        int it_to_it1 = (it - it_begin);
        float t_shift = wavelet_tstart[is] - (it_begin * dt + t0);

        // order as: istage, is, it for better localization
        //    not work because it_begin may diff for diff source
        // int iptr_it = is_local * max_stage + it_to_it1 * max_stage * num_of_src_here;
        // use istage, it, is order
        int iptr_it = is_local * max_nt * max_stage + it_to_it1 * max_stage;

        for (int istage=0; istage<max_stage; istage++)
        {
          // time relative to start time of this source, considering diff from int conversion
          float t = it_to_it1 * dt + rk_stage_time[istage] * dt - t_shift;

          float stf_val;
          if (strcmp(wavelet_name[is], "ricker")==0) {
            stf_val = fun_ricker(t, wavelet_coefs[is][0],wavelet_coefs[is][1]);
          } else if (strcmp(wavelet_name[is], "gaussian")==0) {
            stf_val = fun_gauss(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else if (strcmp(wavelet_name[is], "ricker_deriv")==0) {
            stf_val = fun_ricker_deriv(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else if (strcmp(wavelet_name[is], "gaussian_deriv")==0) {
            stf_val = fun_gauss_deriv(t, wavelet_coefs[is][0], wavelet_coefs[is][1]);
          } else{
            fprintf(stderr,"wavelet_name=%s\n", wavelet_name[is]); 
            fprintf(stderr,"   not implemented yet\n"); 
            fflush(stderr);
          }

          int iptr = iptr_it + istage;

          if ( is < in_num_force) {
            src->Fx[iptr]  = stf_val * force_vector[is][0];
            src->Fy[iptr]  = stf_val * force_vector[is][1];
            src->Fz[iptr]  = stf_val * force_vector[is][2];
          }

          if ( is >= in_num_force) {
            src->Mxx[iptr] = stf_val * moment_tensor[is-in_num_force][0];
            src->Myy[iptr] = stf_val * moment_tensor[is-in_num_force][1];
            src->Mzz[iptr] = stf_val * moment_tensor[is-in_num_force][2];
            src->Mxz[iptr] = stf_val * moment_tensor[is-in_num_force][3];
            src->Myz[iptr] = stf_val * moment_tensor[is-in_num_force][4];
            src->Mxy[iptr] = stf_val * moment_tensor[is-in_num_force][5];
          }
        } // istage
      } // it

      // for extended points and coefs
      sx_inc = source_inc[is][0];
      sy_inc = source_inc[is][1];
      sz_inc = source_inc[is][2];
      float wid_gauss = npoint_half_ext / 2.0;
      float *this_ext_coef = src->ext_coef + is_local * max_ext;
      src_cal_norm_delt3d(this_ext_coef, sx_inc, sy_inc, sz_inc,
                          wid_gauss, wid_gauss, wid_gauss, npoint_half_ext);

      size_t iptr_ext = 0;
      for (int k=sk-npoint_half_ext; k<=sk+npoint_half_ext; k++)
      {
        for (int j=sj-npoint_half_ext; j<=sj+npoint_half_ext; j++)
        {
          for (int i=si-npoint_half_ext; i<=si+npoint_half_ext; i++)
          {
            if (gd_info_lindx_is_inner(i,j,k,gdinfo)==1)
            {
              // Note index need match coef
              int iptr_grid = i + j * siz_line + k * siz_slice;
              int iptr_coef =  (i-(si-npoint_half_ext))
                              + len_ext * (j-(sj-npoint_half_ext)) 
                              + len_ext * len_ext *(k-(sk-npoint_half_ext));
              src->ext_indx[iptr_ext + is_local * max_ext] = iptr_grid;
              src->ext_coef[iptr_ext + is_local * max_ext] = this_ext_coef[iptr_coef];
              iptr_ext++;
            }
          }
        }
      }
      // only count index inside phys region for this thread
      src->ext_num[is_local] = iptr_ext;

      is_local += 1;
    }
  }

  // free
  free(wrk3d);
  free(wavelet_tstart);
  free(in_source_name);
  fdlib_mem_free_2l_float(source_coords, in_num_source, "free_source_coords");
  fdlib_mem_free_2l_float(wavelet_coefs, in_num_source, "free_wavelet_coefs");
  fdlib_mem_free_2l_float(force_vector, in_num_force, "free_force_vector");
  fdlib_mem_free_2l_float(moment_tensor, in_num_moment, "free_moment_tensor");
  fdlib_mem_free_2l_float(source_inc, in_num_source, "free_source_inc");
  fdlib_mem_free_2l_int(source_index, in_num_source, "free_source_index");
  fdlib_mem_free_2l_char(wavelet_name,in_num_source,"free_wavelet_name");
  fdlib_mem_free_2l_char(moment_wavelet_mechism,in_num_moment,"free_wavelet_name");

  return ierr;
}


/*
 * 3d spatial smoothing
 */

void
src_cal_norm_delt3d(float *delt, float x0, float y0, float z0,
                    float rx0, float ry0, float rz0, int LenDelt)
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
src_coord_to_local_indx(gdinfo_t *gdinfo,
                        gdcurv_t *gdcurv,
                        float sx, float sy, float sz,
                        int *si, int *sj, int *sk,
                        float *sx_inc, float *sy_inc, float *sz_inc,
                        float *restrict wrk3d)
{
  int is_here = 0; // default outside

  int nx = gdinfo->nx;
  int ny = gdinfo->ny;
  int nz = gdinfo->nz;
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line  = gdinfo->siz_iy;
  size_t siz_slice = gdinfo->siz_iz;

  float *restrict x3d = gdcurv->x3d;
  float *restrict y3d = gdcurv->y3d;
  float *restrict z3d = gdcurv->z3d;

  // outside coord range
  if ( sx < gdcurv->xmin || sx > gdcurv->xmax ||
       sy < gdcurv->ymin || sy > gdcurv->ymax ||
       sz < gdcurv->zmin || sz > gdcurv->zmax)
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

