/*
 * source term related processing
 */

// todo:

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_math.h"
#include "interp.h"
#include "fdlib_mem.h"
#include "src_t.h"

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
  src->Myz = NULL;
  src->Mxz = NULL;
  src->Mxy = NULL;

  if (moment_actived == 1) {
    src->Mxx= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Myy= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mzz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Myz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mxz= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    src->Mxy= (float *)malloc(max_stage * max_nt * num_of_src * sizeof(float));
    for (int iptr=0; iptr < max_stage * max_nt * num_of_src; iptr++) {
      src->Mxx[iptr] = 0.0;
      src->Myy[iptr] = 0.0;
      src->Mzz[iptr] = 0.0;
      src->Myz[iptr] = 0.0;
      src->Mxz[iptr] = 0.0;
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
               gd_t *gd,
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
      if (gd->type == GD_TYPE_CURV)
      {
        gd_curv_coord_to_glob_indx(gdinfo,gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
      else if (gd->type == GD_TYPE_CART)
      {
        gd_cart_coord_to_glob_indx(gdinfo,gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
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

  // print for QC
  if (myid == 0 && verbose > 2) {
    fprintf(stdout,"src located results:\n");
    for (int is=0; is < in_num_of_src; is++)
    {
      fprintf(stdout,"-- %d: coord=(%f,%f,%f), indx=(%d,%d,%d), inc=(%f,%f,%f)\n",
                    is,source_coords[is][0],source_coords[is][1],source_coords[is][2],
                    source_index[is][0],source_index[is][1],source_index[is][2],
                    source_inc[is][0],source_inc[is][1],source_inc[is][2]);
    }
    fflush(stdout);
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
          //moment_tensor: 0->Mxx 1>Myy 2->Mzz 3->Myz 4->Mxz 5->Mxy
          if (source_moment_actived[is]==1) {
            src->Mxx[iptr] = stf_val * moment_tensor[is][0];
            src->Myy[iptr] = stf_val * moment_tensor[is][1];
            src->Mzz[iptr] = stf_val * moment_tensor[is][2];
            src->Myz[iptr] = stf_val * moment_tensor[is][3];
            src->Mxz[iptr] = stf_val * moment_tensor[is][4];
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
                       gd_t *gd,
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
  if ((fp = fopen(pfilepath, "r"))==NULL)
  {
    fprintf(stdout,"fail to open valsrc file\n");
    fprintf(stdout,"must read only src file, code break now\n");
    free(in_source_name);
    exit(1);
  }
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

  //
  // first run: loop all src to get info
  //
  float *wavelet_tstart = NULL;
  int *source_in_thread = NULL;
  int **source_index = NULL;
  float **source_inc = NULL;
  if(in_num_source > 0)
  {
    wavelet_tstart = (float *) fdlib_mem_calloc_1d_float(in_num_source, 0.0, "wavelet_tstart");
    source_in_thread = (int *)fdlib_mem_calloc_1d_int(in_num_source,-1, "source_in_this_thread");
    source_index = (int **) fdlib_mem_malloc_2l_int(in_num_source,CONST_NDIM,"source_index");
    source_inc = (float **) fdlib_mem_malloc_2l_float(in_num_source,CONST_NDIM,"source_inc");
  }

  int num_of_src_here = 0;
  int force_actived = 0;
  int moment_actived = 0;
  int   max_nt = 0;
  max_nt = (int) (((nt_in-1)*dt_in / dt)+ 0.5) + 1; 

  for (int is=0; is < in_num_source; is++)
  {
    // count num of src in this thread
    
    // convert coord to glob index
      float sx = source_coords[is][0];
      float sy = source_coords[is][1];
      float sz = source_coords[is][2];

      if (gd->type == GD_TYPE_CURV)
      {
        gd_curv_coord_to_glob_indx(gdinfo,gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
      else if (gd->type == GD_TYPE_CART)
      {
        gd_cart_coord_to_glob_indx(gdinfo,gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
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

    //-- to notice user the progress by screen output
    if (myid == 0 && (is % 1000 ==0) ) {
      fprintf(stdout,"-- loc %d-th src index, finish %f\%\n",
                  is, (float)(is+1)/in_num_source*100.0);
      fflush(stdout);
    }
  }

  // print for QC
  if (myid == 0 && verbose > 2) {
    fprintf(stdout,"src located results:\n");
    fprintf(stdout,"in_num_force is %d, in_num_moment is %d\n",in_num_force,in_num_moment);
    for (int is=0; is < in_num_source; is++)
    {
      fprintf(stdout,"-- %d: coord=(%f,%f,%f), indx=(%d,%d,%d), inc=(%f,%f,%f)\n",
                    is,source_coords[is][0],source_coords[is][1],source_coords[is][2],
                    source_index[is][0],source_index[is][1],source_index[is][2],
                    source_inc[is][0],source_inc[is][1],source_inc[is][2]);
                    
      if(source_index[is][0] == -1000 || source_index[is][1] == -1000 || 
         source_index[is][2] == -1000)
        {
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"######### Warning ########\n");
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"source number %d physical coordinates are outside calculation area !\n",is);
        }
    }
    fflush(stdout);
  }

  // check if force and moment used
  if (in_num_force > 0 && num_of_src_here > 0) force_actived = 1;
  if (in_num_moment > 0 && num_of_src_here > 0) moment_actived = 1;
  fprintf(stdout,"force_actived is %d, moment_actived is %d\n",force_actived,moment_actived);
  fprintf(stdout,"num_of_src_here is %d, myid is %d\n",num_of_src_here,myid);
  fflush(stdout);

  float ***force_vector = NULL;
  float ***moment_tensor = NULL;
  char  **moment_wavelet_mechism = NULL;
  if(in_num_moment>0)
  {
    moment_wavelet_mechism = fdlib_mem_malloc_2l_char(in_num_moment, 500, "moment_wavelet_mechism");
  }
  if(num_of_src_here > 0)
  {
    force_vector = (float ***) fdlib_mem_calloc_3l_float(num_of_src_here, CONST_NDIM, nt_in, 0.0, "force_vector");
    moment_tensor = (float ***) fdlib_mem_calloc_3l_float(num_of_src_here, CONST_NDIM_2, nt_in, 0.0, "moment_tensor");
  }

  if (num_of_src_here>0)
  {
    int is_local = 0;
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
          sscanf(str,"%f %f %f",&force_vector[is_local][0][k],&force_vector[is_local][1][k],&force_vector[is_local][2][k]);
        }
        is_local += 1;
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
        //moment_tensor: 0->Mxx 1>Myy 2->Mzz 3->Myz 4->Mxz 5->Mxy
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
            moment_tensor[is_local][0][k] = M0*temp_moment[0];
            moment_tensor[is_local][1][k] = M0*temp_moment[1];
            moment_tensor[is_local][2][k] = M0*temp_moment[2];
            moment_tensor[is_local][3][k] = M0*temp_moment[3];
            moment_tensor[is_local][4][k] = M0*temp_moment[4];
            moment_tensor[is_local][5][k] = M0*temp_moment[5];
          }
        }
        if (strcmp("moment_tensor",moment_wavelet_mechism[index]) == 0)
        {
          for (int k=0;k<nt_in;k++)
          {
            fgets(str,500,fp);
            sscanf(str,"%f %f %f %f %f %f",&moment_tensor[is_local][0][k], &moment_tensor[is_local][1][k], 
                                           &moment_tensor[is_local][2][k], &moment_tensor[is_local][3][k],
                                           &moment_tensor[is_local][4][k], &moment_tensor[is_local][5][k]);
          }
        }
        is_local += 1;
      }

      if (source_in_thread[is] == -1)
      {
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
  if(num_of_src_here > 0)
  {
    for (int is=0; is < in_num_source; is++)
    {
      // check if in this thread
      if (source_in_thread[is] == 1)
      {
        si_glob = source_index[is][0];
        sj_glob = source_index[is][1];
        sk_glob = source_index[is][2];
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
            this_Fx = LagInterp_Piecewise_1d(t_in, &force_vector[is_local][0][0], nt_in, order, 
                                                                 wavelet_tstart[is], dt_in, t);
            this_Fy = LagInterp_Piecewise_1d(t_in, &force_vector[is_local][1][0], nt_in, order, 
                                                                 wavelet_tstart[is], dt_in, t);
            this_Fz = LagInterp_Piecewise_1d(t_in, &force_vector[is_local][2][0], nt_in, order, 
                                                                 wavelet_tstart[is], dt_in, t);

            int iptr = iptr_it + istage;

            src->Fx[iptr]  = this_Fx;
            src->Fy[iptr]  = this_Fy;
            src->Fz[iptr]  = this_Fz;
            }

            if ( is >= in_num_force) {
            float this_Mxx, this_Myy, this_Mzz, this_Mxz, this_Myz, this_Mxy;
            // interp1d give t to get moment tensor
            //moment_tensor: 0->Mxx 1>Myy 2->Mzz 3->Myz 4->Mxz 5->Mxy
            this_Mxx = LagInterp_Piecewise_1d(t_in, &moment_tensor[is_local][0][0], nt_in, order,
                                                           wavelet_tstart[is], dt_in, t);
            this_Myy = LagInterp_Piecewise_1d(t_in, &moment_tensor[is_local][1][0], nt_in, order,
                                                           wavelet_tstart[is], dt_in, t);
            this_Mzz = LagInterp_Piecewise_1d(t_in, &moment_tensor[is_local][2][0], nt_in, order,
                                                           wavelet_tstart[is], dt_in, t);
            this_Myz = LagInterp_Piecewise_1d(t_in, &moment_tensor[is_local][3][0], nt_in, order,
                                                           wavelet_tstart[is], dt_in, t);
            this_Mxz = LagInterp_Piecewise_1d(t_in, &moment_tensor[is_local][4][0], nt_in, order,
                                                           wavelet_tstart[is], dt_in, t);
            this_Mxy = LagInterp_Piecewise_1d(t_in, &moment_tensor[is_local][5][0], nt_in, order,
                                                           wavelet_tstart[is], dt_in, t);

            int iptr = iptr_it + istage;
            src->Mxx[iptr] = this_Mxx;
            src->Myy[iptr] = this_Myy;
            src->Mzz[iptr] = this_Mzz;
            src->Myz[iptr] = this_Myz;
            src->Mxz[iptr] = this_Mxz;
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
  }
  // free
  free(in_source_name);
  free(t_in);
  if(in_num_source > 0)
  {
    free(source_in_thread);
    free(wavelet_tstart);
    fdlib_mem_free_2l_float(source_coords, in_num_source, "free_source_coords");
    fdlib_mem_free_2l_float(source_inc, in_num_source, "free_source_inc");
    fdlib_mem_free_2l_int(source_index, in_num_source, "free_source_index");
  }
  if(num_of_src_here > 0)
  {
    fdlib_mem_free_3l_float(force_vector, num_of_src_here, CONST_NDIM, "free_force_vector");
    fdlib_mem_free_3l_float(moment_tensor, num_of_src_here, CONST_NDIM_2, "free_moment_tensor");
  }
  if(in_num_moment > 0)
  {
    fdlib_mem_free_2l_char(moment_wavelet_mechism,in_num_moment,"free_wavelet_name");
  }
  return ierr;
}



int
src_read_locate_anasrc(gdinfo_t *gdinfo,
                       gd_t *gd,
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
  if ((fp = fopen(pfilepath, "r"))==NULL) 
  {
    fprintf(stdout,"fail to open anasrc file\n");
    fprintf(stdout,"must read only src file, code break now\n");
    free(in_source_name);
    exit(1);
  }
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
  int *source_in_thread = NULL;
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
    source_in_thread = (int *) fdlib_mem_calloc_1d_int(in_num_source,-1, "source_in_this_thread");
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
    for (int is=0;is<in_num_force;is++)
    {
      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&source_coords[is][0],&source_coords[is][1],&source_coords[is][2]);

      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&force_vector[is][0], &force_vector[is][1], &force_vector[is][2]);

      fgets(str,500,fp);
      sscanf(str,"%s",wavelet_name[is]);

      fgets(str,500,fp);
      sscanf(str,"%f %f",&wavelet_coefs[is][0], &wavelet_coefs[is][1]);

      fgets(str,500,fp);
      sscanf(str,"%f",&wavelet_tstart[is]);
    }
  }
  // read moment
  if (in_num_moment > 0)
  {
    for(int is=0;is<in_num_moment;is++)
    {
      fgets(str,500,fp);
      sscanf(str,"%f %f %f",&source_coords[is+in_num_force][0],&source_coords[is+in_num_force][1],
                            &source_coords[is+in_num_force][2]);

      fgets(str,500,fp);
      sscanf(str, "%s", moment_wavelet_mechism[is]);

      //moment_tensor: 0->Mxx 1>Myy 2->Mzz 3->Myz 4->Mxz 5->Mxy
      if (strcmp("mechanism_angle",moment_wavelet_mechism[is])==0)
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
          moment_tensor[is][j]=M0*temp_moment[j];
        }
      } 
      if (strcmp("moment_tensor",moment_wavelet_mechism[is])==0)
      {  
        fgets(str,500,fp);
        sscanf(str,"%f %f %f %f %f %f",&moment_tensor[is][0],&moment_tensor[is][1],
            &moment_tensor[is][2],&moment_tensor[is][3],&moment_tensor[is][4],&moment_tensor[is][5]);
      } 
      fgets(str,500,fp);
      sscanf(str,"%s",wavelet_name[is+in_num_force]);

      fgets(str,500,fp);
      sscanf(str,"%f %f",&wavelet_coefs[is+in_num_force][0], &wavelet_coefs[is+in_num_force][1]);

      fgets(str,500,fp);
      sscanf(str,"%f",&wavelet_tstart[is+in_num_force]);
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

  //
  // first run: loop all src to get info
  //
  int num_of_src_here = 0;
  int force_actived = 0;
  int moment_actived = 0;
  int   max_nt = 0;
  max_nt = (int) (wavelet_tlen / dt + 0.5) + 1; 

  for (int is=0; is < in_num_source; is++)
  {
    // count num of src in this thread
    
    // convert coord to glob index
      float sx = source_coords[is][0];
      float sy = source_coords[is][1];
      float sz = source_coords[is][2];
      if (gd->type == GD_TYPE_CURV)
      {
        gd_curv_coord_to_glob_indx(gdinfo,gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
      else if (gd->type == GD_TYPE_CART)
      {
        gd_cart_coord_to_glob_indx(gdinfo,gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
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
    //-- to notice user the progress by screen output
    if (myid == 0 && (is % 1000 ==0) ) {
      fprintf(stdout,"-- loc %d-th src index, finish %f\%\n",
                  is, (float)(is+1)/in_num_source*100.0);
      fflush(stdout);
    }
  }

  // print for QC
  if (myid == 0 && verbose > 2) {
    fprintf(stdout,"src located results:\n");
    fprintf(stdout,"in_num_force is %d, in_num_moment is %d\n",in_num_force,in_num_moment);
    for (int is=0; is < in_num_source; is++)
    {
      fprintf(stdout,"-- %d: coord=(%f,%f,%f), indx=(%d,%d,%d), inc=(%f,%f,%f)\n",
                    is,source_coords[is][0],source_coords[is][1],source_coords[is][2],
                    source_index[is][0],source_index[is][1],source_index[is][2],
                    source_inc[is][0],source_inc[is][1],source_inc[is][2]);
      if(source_index[is][0] == -1000 || source_index[is][1] == -1000 || 
         source_index[is][2] == -1000)
        {
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"######### Warning ########\n");
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"source number %d physical coordinates are outside calculation area !\n",is);
        }
    }
    fflush(stdout);
  }
  // check if force and moment used
  if (in_num_force > 0 && num_of_src_here > 0) force_actived = 1;
  if (in_num_moment > 0 && num_of_src_here > 0) moment_actived = 1;
  fprintf(stdout,"force_actived is %d, moment_actived is %d\n",force_actived,moment_actived);
  fprintf(stdout,"num_of_src_here is %d, myid is %d\n",num_of_src_here,myid);
  fflush(stdout);
  // alloc src_t
  src_init(src,force_actived,moment_actived,num_of_src_here,max_nt,max_stage,max_ext);

  //
  // second run to set each src in this thread
  //
  int is_local = 0;
  if(num_of_src_here > 0)
  {
    for (int is=0; is < in_num_source; is++)
    {
      // check if in this thread
      if (source_in_thread[is] == 1)
      {
        si_glob = source_index[is][0];
        sj_glob = source_index[is][1];
        sk_glob = source_index[is][2];
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

            //moment_tensor: 0->Mxx 1>Myy 2->Mzz 3->Myz 4->Mxz 5->Mxy
            if ( is >= in_num_force) {
              src->Mxx[iptr] = stf_val * moment_tensor[is-in_num_force][0];
              src->Myy[iptr] = stf_val * moment_tensor[is-in_num_force][1];
              src->Mzz[iptr] = stf_val * moment_tensor[is-in_num_force][2];
              src->Myz[iptr] = stf_val * moment_tensor[is-in_num_force][3];
              src->Mxz[iptr] = stf_val * moment_tensor[is-in_num_force][4];
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
  }
  // free
  free(in_source_name);
  if(in_num_source > 0)
  {
    free(wavelet_tstart);
    free(source_in_thread);
    fdlib_mem_free_2l_float(source_coords, in_num_source, "free_source_coords");
    fdlib_mem_free_2l_float(wavelet_coefs, in_num_source, "free_wavelet_coefs");
    fdlib_mem_free_2l_float(source_inc, in_num_source, "free_source_inc");
    fdlib_mem_free_2l_int(source_index, in_num_source, "free_source_index");
    fdlib_mem_free_2l_char(wavelet_name,in_num_source,"free_wavelet_name");
  }
  if(in_num_force > 0)
  {
    fdlib_mem_free_2l_float(force_vector, in_num_force, "free_force_vector");
  }
  if(in_num_moment > 0)
  {
    fdlib_mem_free_2l_float(moment_tensor, in_num_moment, "free_moment_tensor");
    fdlib_mem_free_2l_char(moment_wavelet_mechism,in_num_moment,"free_wavelet_name");
  }

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
  //moment_tensor: 0->Mxx 1>Myy 2->Mzz 3->Myz 4->Mxz 5->Mxy
  // x->2 y->1 z->3
  source_moment_tensor[0] =  M22 ;  // Mxx
  source_moment_tensor[1] =  M11 ;  // Myy 
  source_moment_tensor[2] =  M33 ;  // Mzz
  source_moment_tensor[3] = -M13 ;  // Myz 
  source_moment_tensor[4] = -M23 ;  // Mxz
  source_moment_tensor[5] =  M12 ;  // Mxy 

  return;
}

