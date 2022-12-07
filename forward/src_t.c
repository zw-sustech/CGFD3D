/*
 */

// todo:

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "netcdf.h"

#include "fdlib_math.h"
#include "interp.h"
#include "fdlib_mem.h"
#include "io_funcs.h"
#include "src_t.h"

//#define DEBUG_SRC_T

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
src_glob_ext_ishere(int si, int sj, int sk, int half_ext, gd_t *gdinfo)
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
 * read .src file and convert into internal structure
 */

int
src_read_locate_file(
                     gd_t     *gd,
                     src_t    *src,
                     float    *restrict mu3d,
                     char     *in_src_file,
                     float     t0,
                     float     dt,
                     int       max_stage,
                     float    *rk_stage_time,
                     int       npoint_half_ext,
                     MPI_Comm  comm,
                     int       myid,
                     int       verbose)
{
  int ierr = 0;

  // get grid info from gdinfo
  int   ni1 = gd->ni1;
  int   ni2 = gd->ni2;
  int   nj1 = gd->nj1;
  int   nj2 = gd->nj2;
  int   nk1 = gd->nk1;
  int   nk2 = gd->nk2;
  int   nx  = gd->nx ;
  int   ny  = gd->ny ;
  int   nz  = gd->nz ;
  int   npoint_ghosts = gd->npoint_ghosts;
  size_t siz_line = gd->siz_iy;
  size_t siz_slice= gd->siz_iz;

  // get total elem of exted src region for a single point
  //    int max_ext = 7 * 7 * 7;
  int len_ext = 2*npoint_half_ext+1;
  int max_ext = len_ext * len_ext * len_ext;

  // local
  FILE *fp =NULL;
  char str[500];
  //char *in_source_name = (char *) malloc(sizeof(char)*500);

  // numbe of source, could be force and/or moment
  int in_num_source;
  // input location is grid index (0) or coordinate (1)
  int is_location_coord;
  // if above value is 1, the 3rd coord is coord (0) or depth (1)
  int in_3coord_meaning;
  // stf is specified by wavelet name (0) or values (1)
  int in_stf_given;
  float in_stf_length=0.0;
  float in_stf_dt=0.0;
  int   in_stf_nt=1;
  // which cmp is used, 1: force; 2: moment, 3: force + moment
  int in_cmp_type;
  // moment is given by tensor (0) or angle + mu D A (1)
  int in_mechanism_type;

  // open in_src_file
  if ((fp = fopen(in_src_file, "r"))==NULL) {
    fprintf(stderr,"ERROR: fail to open in_src_file=%s", in_src_file);
    fflush(stderr); exit(1);
  }

  // event name
  if (!io_get_nextline(fp, str,500)) {
    sprintf(src->evtnm,"%s",str);
  }
  // number of source
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d",&in_num_source);
  }
  if (in_num_source <= 0) {
    fprintf(stderr,"ERROR: in_num_source=%d <=0\n", in_num_source);
    fflush(stderr); exit(1);
  }

  // source time function is given by wavelet name or values
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d",&in_stf_given);
    if (in_stf_given == 0) { // by name
      sscanf(str,"%d %f",&in_stf_given, &in_stf_length);
    } else if (in_stf_given == 1) { // by value
      sscanf(str,"%d %f %d",&in_stf_given, &in_stf_dt, &in_stf_nt);
    } else {
      fprintf(stderr, "ERROR: in_stf_given=%d invalid (either 0 or 1)\n", in_stf_given);
      fflush(stderr); exit(1);
    }
  }

  // force and/or moment, and moment by tensor or angle + muDa
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d %d",&in_cmp_type, &in_mechanism_type);
  }

  // meaning of location and the 3rd input if location is given by coord
  if (!io_get_nextline(fp, str,500)) {
    sscanf(str,"%d %d",&is_location_coord, &in_3coord_meaning);
  }

  //
  // loop each source to get locate and index
  //

  int num_of_src_here = 0;

  float sx, sy, sz;
  int   si,sj,sk;
  int   si_glob,sj_glob,sk_glob;
  float sx_inc, sy_inc, sz_inc;

  float **all_inc    = fdlib_mem_calloc_2l_float(in_num_source,CONST_NDIM,0.0,"all_inc");
  int   **all_index = fdlib_mem_calloc_2l_int(in_num_source,CONST_NDIM,0,"all_index");
  int    *all_in_thread = fdlib_mem_calloc_1d_int(in_num_source,0,"all_in_thread");

  // read coords and determine if in this thread
  for (int is=0; is<in_num_source; is++)
  {
    // read in and get src global index
    if (!io_get_nextline(fp, str,500))
    {
      // read in as float value
      sscanf(str,"%f %f %f", &sx, &sy, &sz);

      if (is_location_coord == 1) // physical coord
      {
        // convert to global index
        //    todo: check if out of computational region
        if (gd->type == GD_TYPE_CURV)
        {
          // if sz is depth, convert to axis when it is in this thread
          if (in_3coord_meaning == 1) {
            gd_curv_depth_to_axis(gd,sx,sy,&sz,comm,myid);
          }
          gd_curv_coord_to_glob_indx(gd,sx,sy,sz,comm,myid,
                                 &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
        }
        else if (gd->type == GD_TYPE_CART)
        {
          // if sz is depth, convert to axis
          if (in_3coord_meaning == 1) {
            sz = gd->z1d[gd->nk2] - sz;
          }
          gd_cart_coord_to_glob_indx(gd,sx,sy,sz,comm,myid,
                                 &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
        }
        // keep index to avoid duplicat run
        all_index[is][0] = si_glob;
        all_index[is][1] = sj_glob;
        all_index[is][2] = sk_glob;
        all_inc  [is][0] = sx_inc;
        all_inc  [is][1] = sy_inc;
        all_inc  [is][2] = sz_inc;

        //-- to notice user the progress using screen output for large input
        if (myid == 0 && (is % 1000 ==0) && verbose>99) {
          fprintf(stdout,"-- loc %d-th src index, finish %2.0f%%\n",
                      is, (float)(is+1)/in_num_source*100.0);
          fflush(stdout);
        }
      }
      else // computational coordinate or grid index
      {
        // if sz is relative to surface, convert to normal index
        if (in_3coord_meaning == 1) {
          sz = gd->gnk2 - sz;
        }

        // nearest integer index
        si_glob = (int) (sx + 0.5);
        sj_glob = (int) (sy + 0.5);
        sk_glob = (int) (sz + 0.5);
        // relative shift
        sx_inc = sx - si_glob;
        sy_inc = sy - sj_glob;
        sz_inc = sz - sk_glob;

        all_index[is][0] = si_glob;
        all_index[is][1] = sj_glob;
        all_index[is][2] = sk_glob;
        all_inc  [is][0] = sx_inc;
        all_inc  [is][1] = sy_inc;
        all_inc  [is][2] = sz_inc;
      }
    }

    // check if in this thread using index
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gd)==1)
    {
      all_in_thread[is] = 1;
      num_of_src_here += 1;
    }
  } // is loop

  // print for QC
  if (myid == 0 && verbose > 99) {
    fprintf(stdout,"src located results:\n");
    for (int is=0; is < in_num_source; is++)
    {
      fprintf(stdout,"-- %d: indx=(%d,%d,%d), inc=(%f,%f,%f)\n",
                    is, all_index[is][0],all_index[is][1],all_index[is][2],
                        all_inc[is][0],all_inc[is][1],all_inc[is][2]);
    }
    fflush(stdout);
  }

  if (myid == 0) {
    for (int is=0; is < in_num_source; is++)
    {
      if(all_index[is][0] == -1000 || all_index[is][1] == -1000 || 
         all_index[is][2] == -1000)
        {
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"######### Warning ########\n");
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"source number %d physical coordinates are outside calculation area !\n",is);
        }
    }
    fflush(stdout);
  }
  //
  // alloc src_t struct for this thread
  //

  // check if force and moment used
  int force_actived  = 0;
  int moment_actived = 0;
  if (num_of_src_here > 0)
  {
    if (in_cmp_type == 1 || in_cmp_type == 3) {
      force_actived = 1;
    }

    if (in_cmp_type == 2 || in_cmp_type == 3) {
      moment_actived = 1;
    } 
  }

  // get number of sample for src_t
  int max_nt = 0;
  if (in_stf_given == 0) { // by name
    max_nt = (int) (in_stf_length / dt + 0.5);
  } else { // by value
    max_nt = (int) (((in_stf_nt-1)*in_stf_dt / dt)+ 0.5) + 1; 
  }

  // alloc src_t
  src_init(src,force_actived,moment_actived,num_of_src_here,max_nt,max_stage,max_ext);

  //
  // loop all source and only keep those in this thread
  //

  float wavelet_tstart;
  char  wavelet_name[50]; // assume max size of name is <= 50
  float wavelet_coefs[10]; // assume max number of coef <= 10
  int  it_begin;
  int  it_end;

  float fx,fy,fz;
  float mxx,myy,mzz,myz,mxz,mxy;
  float *f1 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"f1");
  float *f2 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"f2");
  float *f3 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"f3");
  float *m11 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m11");
  float *m22 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m22");
  float *m33 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m33");
  float *m23 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m23");
  float *m13 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m13");
  float *m12 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"m12");
  float *t_in = (float *)malloc(in_stf_nt*sizeof(float));

  int is_local = 0;
  float M0 = 0.0;
  for (int is=0; is<in_num_source; is++)
  {
    // read stf and cmp of each source
    if (in_stf_given == 0) // wavelet name
    {
      // read in stf
      if (!io_get_nextline(fp, str,500))
      {
        if (all_in_thread[is] == 1) // in in this thread
        {
          // read up to 10 coefs, may be less than 10
          sscanf(str,"%f %s %f %f %f %f %f %f %f %f %f %f",
                  &wavelet_tstart, wavelet_name, wavelet_coefs+0,
                  wavelet_coefs+1, wavelet_coefs+2, wavelet_coefs+3,
                  wavelet_coefs+4, wavelet_coefs+5, wavelet_coefs+6,
                  wavelet_coefs+7, wavelet_coefs+8, wavelet_coefs+9);
          //fprintf(stdout,"-- wavelet_tstart=%g, wavelet_name=%s\n",
          //                    wavelet_tstart, wavelet_name);
          //fprintf(stdout,"---- coef[0]=%g, coef[1]=%g, coef[3]=%g\n",
          //                    wavelet_coefs[0],wavelet_coefs[1],wavelet_coefs[2]);
          //fflush(stdout);
        }
      }
      // read in cmp
      if (!io_get_nextline(fp, str,500))
      {
        if (all_in_thread[is] == 1) // in in this thread
        {
          if (in_cmp_type == 1) { // force
            sscanf(str,"%f %f %f",&fx,&fy,&fz);
          } else if (in_cmp_type == 2) { // moment
            sscanf(str,"%f %f %f %f %f %f",&mxx,&myy,&mzz,&myz,&mxz,&mxy);
          } else { // force + moment
            sscanf(str,"%f %f %f %f %f %f %f %f %f",
                             &fx,&fy,&fz,&mxx,&myy,&mzz,&myz,&mxz,&mxy);
          }

          // convert uDA input into moment tensor
          if (moment_actived==1 && in_mechanism_type ==1)
          {
            si = gd_indx_glphy2lcext_i(all_index[is][0], gd);
            sj = gd_indx_glphy2lcext_j(all_index[is][1], gd);
            sk = gd_indx_glphy2lcext_k(all_index[is][2], gd);
            size_t iptr = si + sj * siz_line + sk * siz_slice;
            // mu < 0 means to use internal model mu value
            float mu =  myz;
            if (mu < 0.0) { mu =  mu3d[iptr]; }
            //mxz is D, mxy[it] is A,
            M0 += mu*mxz*mxy;
            src_muDA_to_moment(mxx,myy,mzz,mu,mxz,mxy,
                      &mxx,&myy,&mzz,&myz,&mxz,&mxy);
          }
        }
      }
    }
    else // by values
    {
      // read t0
      if (!io_get_nextline(fp, str,500)) {
        sscanf(str,"%f",&wavelet_tstart);
      }

      // read cmp in number of in_stf_nt no matter in_thread or not
      for (int it=0; it<in_stf_nt; it++)
      {
        if (!io_get_nextline(fp, str,500))
        {
          if (all_in_thread[is] == 1) // in in this thread
          {
            if (in_cmp_type == 1) { // force
              sscanf(str,"%f %f %f",f1+it,f2+it,f3+it);
            } else if (in_cmp_type == 2) { // moment
              sscanf(str,"%f %f %f %f %f %f",m11+it,m22+it,m33+it,m23+it,m13+it,m12+it);
            } else { // force + moment
              sscanf(str,"%f %f %f %f %f %f %f %f %f",
                  f1+it,f2+it,f3+it,m11+it,m22+it,m33+it,m23+it,m13+it,m12+it);
            }

            // convert uDA input into moment tensor
            if (moment_actived==1 && in_mechanism_type ==1)
            {
              si = gd_indx_glphy2lcext_i(all_index[is][0], gd);
              sj = gd_indx_glphy2lcext_j(all_index[is][1], gd);
              sk = gd_indx_glphy2lcext_k(all_index[is][2], gd);
              size_t iptr = si + sj * siz_line + sk * siz_slice;

              // mu < 0 means to use internal model mu value
              float mu =  m23[it];
              if (mu < 0.0) { mu =  mu3d[iptr]; }

              //m13[it] is v, m12[it] is A,
              M0 += mu*m13[it]*in_stf_dt*m12[it];
              src_muDA_to_moment(m11[it],m22[it],m33[it],mu,m13[it],m12[it],
                                 m11+it ,m22+it ,m33+it ,m23+it ,m13+it ,m12+it);
            }
          } // in this thread
        } // get next line
      } // it
    } // read in stf for one is

    // push into src_t if in this thread
    if (all_in_thread[is] == 1)
    {
      // convert global index to local index
      si = gd_indx_glphy2lcext_i(all_index[is][0], gd);
      sj = gd_indx_glphy2lcext_j(all_index[is][1], gd);
      sk = gd_indx_glphy2lcext_k(all_index[is][2], gd);
  
      // keep into src_t
      src->si[is_local] = si;
      src->sj[is_local] = sj;
      src->sk[is_local] = sk;

      // for extended points and coefs
      sx_inc = all_inc[is][0];
      sy_inc = all_inc[is][1];
      sz_inc = all_inc[is][2];
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
              if (gd_lindx_is_inner(i,j,k,gd)==1)
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

      //
      // wavelet
      //

      // time step, considering t0
      it_begin = (int) ( (wavelet_tstart - t0) / dt);
      it_end   = it_begin + max_nt - 1;

      src->it_begin[is_local] = it_begin;
      src->it_end  [is_local] = it_end  ;

      // setp input t vector for interp 
      for(int it=0; it<in_stf_nt; it++)
      {
        t_in[it] = wavelet_tstart + it*in_stf_dt;
      }

      for (int it=it_begin; it<=it_end; it++)
      {
        int it_to_it1 = (it - it_begin);
        int iptr_it = is_local * max_nt * max_stage + it_to_it1 * max_stage;
        // need to explain purpose
        float t_shift = wavelet_tstart - (it_begin * dt + t0);

        for (int istage=0; istage<max_stage; istage++)
        {
          int iptr = iptr_it + istage;

          // cal stf for given wavelet name
          if (in_stf_given==0)
          {
            // time relative to start time of this source, considering diff from int conversion
            float t = it_to_it1 * dt + rk_stage_time[istage] * dt - t_shift;

            float stf_val = src_cal_wavelet(t,dt,wavelet_name,wavelet_coefs);
            if (force_actived==1) {
              src->Fx[iptr]  = stf_val * fx;
              src->Fy[iptr]  = stf_val * fy;
              src->Fz[iptr]  = stf_val * fz;
            }
            if (moment_actived==1) {
              src->Mxx[iptr] = stf_val * mxx;
              src->Myy[iptr] = stf_val * myy;
              src->Mzz[iptr] = stf_val * mzz;
              src->Myz[iptr] = stf_val * myz;
              src->Mxz[iptr] = stf_val * mxz;
              src->Mxy[iptr] = stf_val * mxy;
            }
          }
          // interp for input values
          else
          {
            // time relative to start time of this source, considering diff from int conversion
            float t = it * dt + rk_stage_time[istage] * dt - t_shift;

            // interp1d order
            int order = 3;     

            if (force_actived==1)
            {
              fx = LagInterp_Piecewise_1d(t_in, f1, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              fy = LagInterp_Piecewise_1d(t_in, f2, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              fz = LagInterp_Piecewise_1d(t_in, f3, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);

              src->Fx[iptr]  = fx;
              src->Fy[iptr]  = fy;
              src->Fz[iptr]  = fz;
            }

            if (moment_actived==1)
            {
              mxx = LagInterp_Piecewise_1d(t_in, m11, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              myy = LagInterp_Piecewise_1d(t_in, m22, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              mzz = LagInterp_Piecewise_1d(t_in, m33, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              myz = LagInterp_Piecewise_1d(t_in, m23, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              mxz = LagInterp_Piecewise_1d(t_in, m13, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              mxy = LagInterp_Piecewise_1d(t_in, m12, in_stf_nt, order,
                      wavelet_tstart, in_stf_dt, t);
              src->Mxx[iptr] = mxx;
              src->Myy[iptr] = myy;
              src->Mzz[iptr] = mzz;
              src->Myz[iptr] = myz;
              src->Mxz[iptr] = mxz;
              src->Mxy[iptr] = mxy;
            }
          }
        } // istage
      } // it
      
      // local is increase
      is_local += 1;

    } // if in_thread
  } // is
  if (in_cmp_type == 2 && in_mechanism_type ==1)
  {
    float sendbuf = M0;
    MPI_Allreduce(&sendbuf, &M0, 1, MPI_FLOAT, MPI_SUM, comm);
    float Mw = 2.0/3.0*log10(M0)-6.06;
    fprintf(stdout,"Mw is %f\n",Mw);
  }

  // close file and free local pointer

  fclose(fp); 

  free(f1);
  free(f2);
  free(f3);
  free(m11);
  free(m22);
  free(m33);
  free(m23);
  free(m13);
  free(m12);
  free(t_in);
  free(all_in_thread);

  fdlib_mem_free_2l_float(all_inc, in_num_source, "free_all_inc");
  fdlib_mem_free_2l_int  (all_index, in_num_source, "free_all_index");

  return ierr;
}

/*
 * read src nc file ordered as [t,n] of each cmp, then interp and convert
 *  into internal structure [cmp,stage,is,it]
 */

int
src_dd_read2local(
                  gd_t     *gd,
                  src_t    *src,
                  char     *in_file,
                  char     *tmp_dir,
                  int       dd_is_add_at_point,
                  int       dd_nt_per_read,
                  float     t0,
                  float     dt,
                  int       nt_total,
                  int       max_stage,
                  float    *rk_stage_time,
                  int       npoint_half_ext,
                  MPI_Comm  comm,
                  int       myid,
                  int*      topoid,
                  int       verbose)
{
  int ierr = 0;

  // get grid info from gd
  int   ni1 = gd->ni1;
  int   ni2 = gd->ni2;
  int   nj1 = gd->nj1;
  int   nj2 = gd->nj2;
  int   nk1 = gd->nk1;
  int   nk2 = gd->nk2;
  int   nx  = gd->nx ;
  int   ny  = gd->ny ;
  int   nz  = gd->nz ;
  size_t siz_line = gd->siz_iy;
  size_t siz_slice= gd->siz_iz;

  // set default values
  src->fp_vi  = NULL;
  src->fp_mij = NULL;
  src->dd_is_valid = 0;

  // return if input dd file is empty
  if (strlen(in_file) == 0) {
    return ierr;
  }

  // get values from input
  src->dd_is_add_at_point = dd_is_add_at_point;
  src->dd_nt_per_read     = dd_nt_per_read    ;
  src->max_stage          = max_stage;
  src->dd_it_here         = 0;

  if (src->dd_is_add_at_point==1) {
    src->dd_smo_hlen = 0;
  } else {
    src->dd_smo_hlen = npoint_half_ext;
  }

  // local
  // numbe of source, could be force and/or moment
  int in_num_source;
  // input location is grid index (0) or coordinate (1)
  int is_location_coord;
  // if above value is 1, the 3rd coord is coord (0) or depth (1)
  int in_3coord_meaning;
  // stf is specified by wavelet name (0) or values (1)
  float in_stf_length=0.0;
  float in_stf_dt=0.0;
  int   in_stf_nt=1;

  char filenm[CONST_MAX_STRLEN];

  // which cmp is used, 1: force; 2: moment, 3: force + moment
  int in_cmp_type;
  // moment is given by tensor (0) or angle + mu D A (1)
  int in_mechanism_type;

  int ncid,tdimid,ndimid;
  int xid,yid,zid,tid;
  
  // read in nc
  ierr = nc_open(in_file, NC_NOWRITE, &ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s for file %s\n", nc_strerror(ierr), in_file);
    MPI_Abort(comm,ierr);
  }

  // get att for flags
  nc_get_att_int(ncid,NC_GLOBAL,"location_is_axis",&is_location_coord);
  nc_get_att_int(ncid,NC_GLOBAL,"z_is_depth",&in_3coord_meaning);
  if (myid == 0 && verbose>999) {
   fprintf(stdout,"-- location_is_axis=%d\n", is_location_coord);
   fprintf(stdout,"-- z_is_depth=%d\n", in_3coord_meaning);
   fflush(stdout);
  }

  // get dim size
  size_t tdim_len;
  ierr = nc_inq_dimid(ncid, "time", &tdimid);
  ierr = nc_inq_dimlen(ncid, tdimid, &tdim_len);
  in_stf_nt = tdim_len;

  size_t ndim_len;
  ierr = nc_inq_dimid(ncid, "number", &ndimid);
  ierr = nc_inq_dimlen(ncid, ndimid, &ndim_len);
  in_num_source = ndim_len;

  if (myid == 0 && verbose>999) {
    fprintf(stdout,"-- tdim_len=%d,ndim_len=%d\n", tdim_len,ndim_len);
    fflush(stdout);
  }

  float *t_in = (float *)malloc(in_stf_nt*sizeof(float));
  float *all_x = fdlib_mem_calloc_1d_float(in_num_source,0.0,"all_x");
  float *all_y = fdlib_mem_calloc_1d_float(in_num_source,0.0,"all_y");
  float *all_z = fdlib_mem_calloc_1d_float(in_num_source,0.0,"all_z");

  // get time var
  ierr = nc_inq_varid(ncid, "time", &tid);
  ierr = nc_get_var(ncid,tid,t_in);

  // currently, t_in[0] should be equal to t0
  if (fabs(t_in[0] - t0) > 1e-10)
  {
    fprintf(stderr,"ERROR: fd t0=%g is not equal to t0=%g in %s\n",
            t0, t_in[0], in_file);
    fflush(stderr);
    MPI_Abort(comm,ierr);
  }

  // stf dt
  in_stf_dt = t_in[1] - t_in[0];
  if (myid == 0 && verbose>999) {
    fprintf(stdout,"-- in_stf_dt=%g\n", in_stf_dt);
    fflush(stdout);
  }

  // get number of sample for src_t
  int max_nt = (int) (( (t_in[in_stf_nt-1]-t0) / dt)+ 0.5); 
  if (myid == 0 && verbose>999) {
    fprintf(stdout,"-- max_nt converted from iput nc=%d, nt_total=%d\n",
                    max_nt, nt_total);
    fflush(stdout);
  }

  // should be less than fd total nt
  if (max_nt > nt_total) {
    max_nt = nt_total;
  }
  src->dd_max_nt = max_nt;

  // reset per read if larger than max_nt
  if (src->dd_nt_per_read > max_nt) {
    src->dd_nt_per_read = max_nt;
  }

  // get location vars
  ierr = nc_inq_varid(ncid, "x", &xid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    MPI_Abort(comm,ierr);
  }
  ierr = nc_inq_varid(ncid, "y", &yid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    MPI_Abort(comm,ierr);
  }
  ierr = nc_inq_varid(ncid, "z", &zid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    MPI_Abort(comm,ierr);
  }

  ierr = nc_get_var_float(ncid,xid,all_x);
  ierr = nc_get_var_float(ncid,yid,all_y);
  ierr = nc_get_var_float(ncid,zid,all_z);

  // allocate other loc vars
  float **all_inc    = fdlib_mem_calloc_2l_float(in_num_source,CONST_NDIM,0.0,"all_inc");
  int   **all_index = fdlib_mem_calloc_2l_int(in_num_source,CONST_NDIM,0,"all_index");
  int    *all_in_thread = fdlib_mem_calloc_1d_int(in_num_source,0,"all_in_thread");

  int num_of_src_here = 0;
  int   si,sj,sk;
  int   si_glob,sj_glob,sk_glob;
  float sx_inc, sy_inc, sz_inc;

  //
  // locate
  //

  for (int is=0; is<in_num_source; is++)
  {
    float sx = all_x[is];
    float sy = all_y[is];
    float sz = all_z[is];

    if (is_location_coord == 1) // physical coord
    {
      // convert to global index
      //    todo: check if out of computational region
      if (gd->type == GD_TYPE_CURV)
      {
        // if sz is depth, convert to axis when it is in this thread
        if (in_3coord_meaning == 1) {
          gd_curv_depth_to_axis(gd,sx,sy,&sz,comm,myid);
        }
        gd_curv_coord_to_glob_indx(gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
      else if (gd->type == GD_TYPE_CART)
      {
        // if sz is depth, convert to axis
        if (in_3coord_meaning == 1) {
          sz = gd->z1d[gd->nk2] - sz;
        }
        gd_cart_coord_to_glob_indx(gd,sx,sy,sz,comm,myid,
                               &si_glob,&sj_glob,&sk_glob,&sx_inc,&sy_inc,&sz_inc);
      }
      // keep index to avoid duplicat run
      all_index[is][0] = si_glob;
      all_index[is][1] = sj_glob;
      all_index[is][2] = sk_glob;
      all_inc  [is][0] = sx_inc;
      all_inc  [is][1] = sy_inc;
      all_inc  [is][2] = sz_inc;

      //-- to notice user the progress using screen output for large input
      if (myid == 0 && (is % 1000 ==0) && verbose>1) {
          fprintf(stdout,"-- loc %d-th src index, finish %2.0f%%\n",
                      is, (float)(is+1)/in_num_source*100.0);
        fflush(stdout);
      }
    }
    else // computational coordinate or grid index
    {
      // if sz is relative to surface, convert to normal index
      if (in_3coord_meaning == 1) {
        sz = gd->gnk2 - sz;
      }

      // nearest integer index
      si_glob = (int) (sx + 0.5);
      sj_glob = (int) (sy + 0.5);
      sk_glob = (int) (sz + 0.5);
      // relative shift
      sx_inc = sx - si_glob;
      sy_inc = sy - sj_glob;
      sz_inc = sz - sk_glob;

      all_index[is][0] = si_glob;
      all_index[is][1] = sj_glob;
      all_index[is][2] = sk_glob;
      all_inc  [is][0] = sx_inc;
      all_inc  [is][1] = sy_inc;
      all_inc  [is][2] = sz_inc;
    }

    // check if in this thread using index
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,src->dd_smo_hlen,gd)==1)
    {
      all_in_thread[is] = 1;
      num_of_src_here += 1;
    }
  } // is loop

  // set related values in src_t
  src->dd_total_number = num_of_src_here;

  // print for QC
  if (myid == 0 && verbose > 99) {
    fprintf(stdout,"ddsrc located results:\n");
    for (int is=0; is < in_num_source; is++)
    {
      fprintf(stdout,"-- %d: indx=(%d,%d,%d), inc=(%f,%f,%f)\n",
                    is, all_index[is][0],all_index[is][1],all_index[is][2],
                        all_inc[is][0],all_inc[is][1],all_inc[is][2]);
                    
      if(all_index[is][0] == -1000 || all_index[is][1] == -1000 || 
         all_index[is][2] == -1000)
        {
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"######### Warning ########\n");
          fprintf(stdout,"#########         ########\n");
          fprintf(stdout,"source number %d physical coordinates are outside calculation area !\n",is);
        }
    }
    fflush(stdout);
  }

  // set valid if num_of_src_here > 0
  if (num_of_src_here > 0) {
    src->dd_is_valid = 1;
  }
  if (verbose>999)
  {
    MPI_Barrier(MPI_COMM_WORLD);
    fprintf(stdout,"-- %d,%d: num_of_src_here=%d, dd_is_valid=%d\n", 
                  topoid[0],topoid[1],num_of_src_here,src->dd_is_valid);
    fflush(stdout);
  }
  
  // set index in src_t

  if (src->dd_is_valid==1)
  {
    src->dd_indx = (size_t *)malloc(sizeof(size_t)*num_of_src_here);
    src->dd_si = (int *)malloc(num_of_src_here*sizeof(int));
    src->dd_sj = (int *)malloc(num_of_src_here*sizeof(int));
    src->dd_sk = (int *)malloc(num_of_src_here*sizeof(int));
    src->dd_si_inc = (float *)malloc(num_of_src_here*sizeof(float));
    src->dd_sj_inc = (float *)malloc(num_of_src_here*sizeof(float));
    src->dd_sk_inc = (float *)malloc(num_of_src_here*sizeof(float));
  }

  int is_local = 0;
  for (int is=0; is<in_num_source; is++)
  {
    if (all_in_thread[is] == 1)
    {
      if (verbose > 999) {
        fprintf(stdout,"--- %d %d is=%d : set index\n", 
                            topoid[0],topoid[1],is);
        fflush(stdout);
      }

      // convert global index to local index
      si = gd_indx_glphy2lcext_i(all_index[is][0], gd);
      sj = gd_indx_glphy2lcext_j(all_index[is][1], gd);
      sk = gd_indx_glphy2lcext_k(all_index[is][2], gd);
  
      // keep into src_t
      src->dd_si[is_local] = si;
      src->dd_sj[is_local] = sj;
      src->dd_sk[is_local] = sk;

      // for extended points and coefs
      src->dd_si_inc[is_local] = all_inc[is][0];
      src->dd_sj_inc[is_local] = all_inc[is][1];
      src->dd_sk_inc[is_local] = all_inc[is][2];

      // indx
      src->dd_indx[is_local] = si + sj * siz_line + sk * siz_slice;

      is_local+=1;
    }
  }

  //
  // read stf and convert to local file
  //

  int vxid,vyid,vzid;
  int mxxid,myyid,mzzid;
  int myzid,mxzid,mxyid;

  float v3e[CONST_NDIM];
  float v6e[CONST_NDIM_2];

  // check if force and moment used
  //  check if there is "Fx", "Mxx" etc
  ierr = nc_inq_varid(ncid, "Fx", &vxid);
  if ((ierr == NC_NOERR) && (src->dd_is_valid==1))
  {
    src->dd_vi_actived = 1;
    ierr = nc_inq_varid(ncid, "Fy", &vyid);
    ierr = nc_inq_varid(ncid, "Fz", &vzid);
  }
  else
  {
    src->dd_vi_actived = 0;
  }

  ierr = nc_inq_varid(ncid, "Mxx_rate", &mxxid);
  if ((ierr == NC_NOERR) && (src->dd_is_valid==1))
  {
    src->dd_mij_actived = 1;

    ierr = nc_inq_varid(ncid, "Myy_rate", &myyid);
    ierr = nc_inq_varid(ncid, "Mzz_rate", &mzzid);
    ierr = nc_inq_varid(ncid, "Myz_rate", &myzid);
    ierr = nc_inq_varid(ncid, "Mxz_rate", &mxzid);
    ierr = nc_inq_varid(ncid, "Mxy_rate", &mxyid);
  }
  else
  {
    src->dd_mij_actived = 0;
  }

  // alloc vars
  if (src->dd_vi_actived==1) {
    src->dd_vi = (float *)malloc( src->dd_nt_per_read * max_stage 
                                  * num_of_src_here * CONST_NDIM * sizeof(float));
  }
  if (src->dd_mij_actived==1) {
    src->dd_mij = (float *)malloc( src->dd_nt_per_read * max_stage 
                                  * num_of_src_here * CONST_NDIM_2 * sizeof(float));
  }

  // create local file
  if (src->dd_vi_actived==1)
  {
    io_build_fname(tmp_dir,"src_ddvi",".bin",topoid,filenm);
    if ((src->fp_vi = fopen(filenm,"w+b"))==NULL) {
      fprintf(stderr,"ERROR: fail to create temp dd src file=%s", filenm);
      fflush(stderr);
      MPI_Abort(comm,ierr);
    } else {
      if (verbose > 99) {
        fprintf(stdout,"--- %d %d write srcdd_vi to file : %s\n", 
                            topoid[0],topoid[1],filenm);
      }
    }
  }
  if (src->dd_mij_actived==1)
  {
    io_build_fname(tmp_dir,"src_ddmij",".bin",topoid,filenm);
    if ((src->fp_mij = fopen(filenm,"w+b"))==NULL) {
      fprintf(stderr,"ERROR: fail to create temp dd src file=%s", filenm);
      fflush(stderr);
      MPI_Abort(comm,ierr);
    } else {
      if (verbose > 99) {
        fprintf(stdout,"--- %d %d write srcdd_mij to file : %s\n", 
                            topoid[0],topoid[1],filenm);
      }
    }
  }

  //
  // loop all source and only keep those in this thread
  //

  float *in_f1  = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_f1");
  float *in_f2  = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_f2");
  float *in_f3  = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_f3");
  float *in_m11 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_m11");
  float *in_m22 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_m22");
  float *in_m33 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_m33");
  float *in_m23 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_m23");
  float *in_m13 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_m13");
  float *in_m12 = fdlib_mem_calloc_1d_float(in_stf_nt,0.0,"in_m12");

  float **ou_f1  = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_f1");
  float **ou_f2  = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_f2");
  float **ou_f3  = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_f3");
  float **ou_m11 = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_m11");
  float **ou_m22 = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_m22");
  float **ou_m33 = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_m33");
  float **ou_m23 = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_m23");
  float **ou_m13 = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_m13");
  float **ou_m12 = fdlib_mem_calloc_2l_float(max_nt,max_stage,0.0,"ou_m12");

  is_local = 0;
  for (int is=0; is<in_num_source; is++)
  {
    if (all_in_thread[is] == 1)
    {
      if (verbose > 999) {
        fprintf(stdout,"--- %d %d is=%d : read vi/mij\n", 
                            topoid[0],topoid[1],is);
        fflush(stdout);
      }

      // get data

      if (src->dd_vi_actived == 1)
      {
        //size_t subs[] = {0, is};
        //size_t subc[] = {in_stf_nt, 1};
        size_t subs[] = {is, 0};
        size_t subc[] = {1, in_stf_nt};
        nc_get_vara_float(ncid,vxid,subs,subc,in_f1);
        nc_get_vara_float(ncid,vyid,subs,subc,in_f2);
        nc_get_vara_float(ncid,vzid,subs,subc,in_f3);
      }

      if (src->dd_mij_actived == 1)
      {
        //size_t subs[] = {0, is};
        //size_t subc[] = {in_stf_nt, 1};
        size_t subs[] = {is, 0};
        size_t subc[] = {1, in_stf_nt};
        nc_get_vara_float(ncid,mxxid,subs,subc,in_m11);
        nc_get_vara_float(ncid,myyid,subs,subc,in_m22);
        nc_get_vara_float(ncid,mzzid,subs,subc,in_m33);
        nc_get_vara_float(ncid,myzid,subs,subc,in_m23);
        nc_get_vara_float(ncid,mxzid,subs,subc,in_m13);
        nc_get_vara_float(ncid,mxyid,subs,subc,in_m12);
      }

      if (verbose > 999) {
        fprintf(stdout,"--- %d %d is=%d : read in stf_dt=%g\n", 
                            topoid[0],topoid[1], is, in_stf_dt);
        for (int it=0; it < in_stf_nt; it++) {
          fprintf(stdout,"--- %d %d is=%d : in_it=%d,in_t=%g,mxx=%g,myy=%g,mzz=%g,myz=%g,mxz=%g,mxy=%g\n",
                              topoid[0],topoid[1], is, 
                              it,it*in_stf_dt, 
                              in_m11[it], in_m22[it], in_m33[it], in_m23[it], in_m13[it], in_m12[it]);
        }
        fflush(stdout);
      }

      if (verbose > 999) {
        fprintf(stdout,"--- %d %d is=%d : interp stf from stf_dt=%g to dt=%g\n", 
                            topoid[0],topoid[1], is, in_stf_dt, dt);
        fflush(stdout);
      }


      // interp to fd dt 
      for (int it=0; it<max_nt; it++)
      {
        float t_start = it * dt + t0;

        for (int istage=0; istage<max_stage; istage++)
        {
          float t = t_start + rk_stage_time[istage] * dt;

          // interp1d order
          int order = 3;     

          if (src->dd_vi_actived==1)
          {
            ou_f1[it][istage] = LagInterp_Piecewise_1d(t_in, in_f1, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_f2[it][istage] = LagInterp_Piecewise_1d(t_in, in_f2, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_f3[it][istage] = LagInterp_Piecewise_1d(t_in, in_f3, in_stf_nt, order,
                    t0, in_stf_dt, t);
          }

          if (src->dd_mij_actived==1)
          {
            ou_m11[it][istage] = LagInterp_Piecewise_1d(t_in, in_m11, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_m22[it][istage] = LagInterp_Piecewise_1d(t_in, in_m22, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_m33[it][istage] = LagInterp_Piecewise_1d(t_in, in_m33, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_m23[it][istage] = LagInterp_Piecewise_1d(t_in, in_m23, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_m13[it][istage] = LagInterp_Piecewise_1d(t_in, in_m13, in_stf_nt, order,
                    t0, in_stf_dt, t);
            ou_m12[it][istage] = LagInterp_Piecewise_1d(t_in, in_m12, in_stf_nt, order,
                    t0, in_stf_dt, t);

            if (verbose > 999) {
              fprintf(stdout,"--- %d %d is=%d,it=%d,istage=%d,t=%g,mxx=%g,myy=%g,mzz=%g,myz=%g,mxz=%g,mxy=%g\n",
                                    topoid[0],topoid[1], is, it, istage, t,
                                    ou_m11[it][istage], 
                                    ou_m22[it][istage], 
                                    ou_m33[it][istage], 
                                    ou_m23[it][istage], 
                                    ou_m13[it][istage], 
                                    ou_m12[it][istage]);
              fflush(stdout);
            }

          }
        } // istage
      } // it

      // save to local file, vi
      if (src->dd_vi_actived==1)
      {
        if (verbose>999) {
          fprintf(stdout,"-- write vi to local file\n");
          fflush(stdout);
        }

        for (int it=0; it<max_nt; it++)
        {
          // file pos of is=0 for it
          long int offset = it * max_stage 
                               * num_of_src_here * CONST_NDIM * sizeof(float)
                               +        is_local * CONST_NDIM * sizeof(float);

          fseek(src->fp_vi, offset , SEEK_SET);

          for (int istage=0; istage<max_stage; istage++)
          {
            v3e[0] = ou_f1[it][istage];
            v3e[1] = ou_f2[it][istage];
            v3e[2] = ou_f3[it][istage];

            // write
            fwrite(v3e, sizeof(float), CONST_NDIM, src->fp_vi);

            fseek(src->fp_vi, CONST_NDIM*(num_of_src_here-1)*sizeof(float), SEEK_CUR);
          }
        }
      }

      // save to local file, mij
      if (src->dd_mij_actived==1)
      {
        if (verbose>999) {
          fprintf(stdout,"-- write is=%d,is_loca=%d mij to local file\n",is,is_local);
          fflush(stdout);
        }

        for (int it=0; it<max_nt; it++)
        {
          long int offset = it * max_stage 
                               * num_of_src_here * CONST_NDIM_2 * sizeof(float)
                               +        is_local * CONST_NDIM_2 * sizeof(float);
          fseek(src->fp_mij, offset , SEEK_SET);

          for (int istage=0; istage<max_stage; istage++)
          {
            v6e[0] = ou_m11[it][istage];
            v6e[1] = ou_m22[it][istage];
            v6e[2] = ou_m33[it][istage];
            v6e[3] = ou_m23[it][istage];
            v6e[4] = ou_m13[it][istage];
            v6e[5] = ou_m12[it][istage];

            // write
            fwrite(v6e, sizeof(float), CONST_NDIM_2, src->fp_mij);

            fseek(src->fp_mij, CONST_NDIM_2*(num_of_src_here-1)*sizeof(float), SEEK_CUR);
          }
        }
      }
      
      // local is increase
      is_local += 1;

    } // if in_thread
  } // is

  // close file
  ierr = nc_close(ncid);
  if (ierr != NC_NOERR){
    fprintf(stderr,"nc error: %s\n", nc_strerror(ierr));
    MPI_Abort(comm,ierr);
  }

  if (verbose>999) {
  if (src->dd_is_valid==1)
  {
    fprintf(stdout,"-- summarize of dd src:\n");
    fprintf(stdout,"---- %d %d dd_is_valid=%d,dd_is_add_at_point=%d\n",
                            topoid[0],topoid[1],src->dd_is_valid,src->dd_is_add_at_point);
    fprintf(stdout,"---- %d %d dd_vi_actived=%d,dd_mij_actived=%d\n",
                            topoid[0],topoid[1],src->dd_vi_actived,src->dd_mij_actived);
    fprintf(stdout,"---- %d %d dd_total_number=%d,dd_max_nt=%d,dd_nt_per_read=%d\n",
                            topoid[0],topoid[1],
                            src->dd_total_number,src->dd_max_nt,src->dd_nt_per_read);
    fflush(stdout);
  }
  }

  // do not close dd src tmp file, will read during calculation
  if (src->fp_vi)  fseek(src->fp_vi , 0, SEEK_SET);
  if (src->fp_mij) fseek(src->fp_mij, 0, SEEK_SET);

  // read first part
  // read Vi
  if (src->fp_vi)
  {
    size_t num_this_read =    src->dd_nt_per_read
                            * src->max_stage
                            * src->dd_total_number
                            * CONST_NDIM;
                            
    size_t num_real_read = fread(src->dd_vi, sizeof(float), num_this_read, src->fp_vi);

#ifdef DEBUG_SRC_T
    fprintf(stdout,"-- read srcdd mi: after read, num_real_read=%ld,cur pos=%ld\n",
                                   num_real_read, ftell(src->fp_vi));
    for (int it=0; it < num_this_read; it++)
    {
      fprintf(stdout,"--- %d %d indx=%d,dd_mij=%g\n",
                                      topoid[0],topoid[1], it,
                                      src->dd_vi[it]);
    }
    fflush(stdout);
#endif
  }
  // read Mij
  if (src->fp_mij)
  {
    size_t num_this_read =    src->dd_nt_per_read
                            * src->max_stage
                            * src->dd_total_number
                            * CONST_NDIM_2;
#ifdef DEBUG_SRC_T
    fprintf(stdout,"-- read srcdd mij: num_this_read=%zu,cur pos=%ld\n",
                      num_this_read,ftell(src->fp_mij));
#endif

    size_t num_real_read = fread(src->dd_mij, sizeof(float), num_this_read, src->fp_mij);

#ifdef DEBUG_SRC_T
    fprintf(stdout,"-- read srcdd mij: after read, num_real_read=%ld,cur pos=%ld\n",
                                   num_real_read, ftell(src->fp_mij));
    for (int it=0; it < num_this_read; it++)
    {
      fprintf(stdout,"--- %d %d indx=%d,dd_mij=%g\n",
                                      topoid[0],topoid[1], it,
                                      src->dd_mij[it]);
    }
    fflush(stdout);
#endif
  }

  // free tmp vars
  free(in_f1);
  free(in_f2);
  free(in_f3);
  free(in_m11);
  free(in_m22);
  free(in_m33);
  free(in_m23);
  free(in_m13);
  free(in_m12);

  free(t_in);
  free(all_in_thread);

  fdlib_mem_free_2l_float(all_inc, in_num_source, "free_all_inc");
  fdlib_mem_free_2l_int  (all_index, in_num_source, "free_all_index");

  fdlib_mem_free_2l_float(ou_f1 , max_nt, "ou_f1 ");
  fdlib_mem_free_2l_float(ou_f2 , max_nt, "ou_f2 ");
  fdlib_mem_free_2l_float(ou_f3 , max_nt, "ou_f3 ");
  fdlib_mem_free_2l_float(ou_m11, max_nt, "ou_m11");
  fdlib_mem_free_2l_float(ou_m22, max_nt, "ou_m22");
  fdlib_mem_free_2l_float(ou_m33, max_nt, "ou_m33");
  fdlib_mem_free_2l_float(ou_m23, max_nt, "ou_m23");
  fdlib_mem_free_2l_float(ou_m13, max_nt, "ou_m13");
  fdlib_mem_free_2l_float(ou_m12, max_nt, "ou_m12");

  return ierr;
}

int
src_dd_free(src_t *src)
{
  if (src->fp_vi ) fclose(src->fp_vi );
  if (src->fp_mij) fclose(src->fp_mij);
}

/*
 * increase it, reload stf if required
 */

int
src_dd_accit_loadstf(src_t *src, int it, int myid)
{
  int ierr = 0;

#ifdef DEBUG_SRC_T
  fprintf(stdout,"-- reload srcdd: it=%d,dd_is_valid=%d,dd_max_nt=%d\n",
                      it,src->dd_is_valid,src->dd_max_nt);
  fflush(stdout);
#endif

  // set invlid if out of range
  if (src->dd_is_valid==1 && it >= src->dd_max_nt) {
    src->dd_is_valid = 0;
    return ierr;
  }

  // increase src it
  if (it>0) {
    src->dd_it_here += 1;
  } else {
    src->dd_it_here  = 0;
  }
  // reload data if required
  if (src->dd_it_here >= src->dd_nt_per_read)
  {
#ifdef DEBUG_SRC_T
    fprintf(stdout,"-- reload srcdd: myid=%d,it=%d,dd_it_here=%d,dd_nt_per_read=%d\n",
                        myid,it,src->dd_it_here,src->dd_nt_per_read);
    fflush(stdout);
#endif
    // if read out of range, reduce
    if ( (src->dd_nt_per_read + it) > src->dd_max_nt )
    {
      src->dd_nt_this_read = src->dd_max_nt - it;
    } else {
      src->dd_nt_this_read = src->dd_nt_per_read;
    }

    // read Vi
    if (src->dd_vi_actived == 1)
    {
      size_t num_this_read =    src->dd_nt_this_read
                              * src->max_stage
                              * src->dd_total_number
                              * CONST_NDIM;

      size_t num_real_read = fread(src->dd_vi, sizeof(float), num_this_read, src->fp_vi);

#ifdef DEBUG_SRC_T
      fprintf(stdout,"-- reload srcdd vi: myid=%d,it=%d,it_here=%d,dd_nt_this_read=%d,num_this_read=%zu,num_real_read=%zu\n",
                      myid,it,src->dd_it_here,src->dd_nt_this_read,num_this_read,num_real_read);
      fflush(stdout);
#endif
    }

    // read Mij
    if (src->dd_mij_actived == 1)
    {
      size_t num_this_read =    src->dd_nt_this_read
                              * src->max_stage
                              * src->dd_total_number
                              * CONST_NDIM_2;

      size_t num_real_read = fread(src->dd_mij, sizeof(float), num_this_read, src->fp_mij);

#ifdef DEBUG_SRC_T
      fprintf(stdout,"-- reload srcdd mij: myid=%d,it=%d,it_here=%d,dd_nt_this_read=%d,num_this_read=%zu,num_real_read=%zu\n",
                      myid,it,src->dd_it_here,src->dd_nt_this_read,num_this_read,num_real_read);
      fflush(stdout);
#endif
    }

    // reset dd_it_here
    src->dd_it_here = 0;
  }

  return 0;
}

/*
 * get stf value at a given t
 */

float
src_cal_wavelet(float t, float dt, char *wavelet_name, float *wavelet_coefs)
{
  float stf_val;

  if (strcmp(wavelet_name, "ricker")==0) {
    stf_val = fun_ricker(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "gaussian")==0) {
    stf_val = fun_gauss(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "ricker_deriv")==0) {
    stf_val = fun_ricker_deriv(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "gaussian_deriv")==0) {
    stf_val = fun_gauss_deriv(t, wavelet_coefs[0], wavelet_coefs[1]);
  } else if (strcmp(wavelet_name, "bell")==0) {
    stf_val = fun_bell(t, wavelet_coefs[0]);
  } else if (strcmp(wavelet_name, "bell_deriv")==0) {
    stf_val = fun_bell_deriv(t, wavelet_coefs[0]);
  } else if (strcmp(wavelet_name, "step")==0) {
    stf_val = fun_step(t);
  } else if (strcmp(wavelet_name, "delta")==0) {
    stf_val = fun_delta(t, dt);
  } else if (strcmp(wavelet_name, "liuetal2006")==0) {
    stf_val = fun_liuetal2006(t, wavelet_coefs[0]);
  } else if (strcmp(wavelet_name, "klauder")==0) {
    stf_val = fun_klauder(t, wavelet_coefs[0], wavelet_coefs[1], wavelet_coefs[2], wavelet_coefs[3]);
  } else if (strcmp(wavelet_name, "klauder_blackman")==0) {
    stf_val = fun_klauder_blackman(t, wavelet_coefs[0], wavelet_coefs[1], 
                                      wavelet_coefs[2], wavelet_coefs[3], dt);
  } else{
    fprintf(stderr,"wavelet_name=%s\n", wavelet_name); 
    fprintf(stderr,"   not implemented yet\n"); 
    fflush(stderr); exit(1);
  }

  return stf_val;
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
  float u = (t-t0)*2.0*PI*fc;
  float v = (1-u*u/2.0)*exp(-u*u/4.0);

  return v;
}

float 
fun_ricker_deriv(float t, float fc, float t0)
{
  //float pi = acos(-1.0);
  float u = (t-t0)*2.0*PI*fc;
  float v = u*(-3+1/2*u*u)*exp(-u*u/4)*PI*fc;

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

float
fun_bell(float t, float riset)
{
  float f;

  if (t>=0.0 & t <= riset) {
    f = 2.0 * PI * sin(2*PI*t/riset) / riset;
  } else {
    f = 0.0;
  }

  return f;
}

float
fun_bell_deriv(float t, float riset)
{
  float f;

  if (t < 0.0) {
    f = 0.0;
  } else if (t <= riset) {
    f = t / riset - sin(2.0*PI*t/riset) / (2.0*PI);
  } else {
    f = 1.0;
  }

  return f;
}

float
fun_step(float t)
{
  float f;

  if (t<0.0) {
    f = 0.0;
  } else {
    f = 1.0;
  }

  return f;
}

float
fun_delta(float t, float dt)
{
  float f;

  if (t<0.0 | t > dt) {
    f = 0.0;
  } else {
    f = 1.0 / dt;
  }

  return f;
}

float
fun_klauder(float t, float tshift, float f1, float f2, float T)
{
  float f;

  float K  = (f2-f1)/2.0;
  float fM = (f2+f1)/2.0;

  float t_zero = t - tshift;

  // need to verify the exact eqn
  // need to check if specil treatment is needed when t_zero close to 0
  f = (sin( PI*K*t_zero * (T-t_zero)) * cos(2*PI*fM*t_zero)) / (PI*K*t_zero);

  return f;
}

float
fun_klauder_blackman(float t, float tshift, float f1, float f2, float T, float dt)
{
  float f;

  f = fun_klauder(f,tshift,f1,f2,T);

  float B = blackman_window(t,dt,tshift);

  f = f * B;

  return f;
}

float
blackman_window(float t, float dt, float tshift)
{
    float B;
    float i = t/dt;
    float n = tshift/dt;

    if (i<=2*n) {
      B = 0.42-0.5*cos(2*PI*(i-1)/(2*n-1)) + 0.08*cos(4*PI*(i-1)/(2*n-1));
    } else {
      B = 0.0;
    }

    return B;
}


// eqn 7a in liu et al., 2006
float
fun_liuetal2006(float t, float tau)
{
  float f;

  float tau_1 = 0.13 * tau;
  float tau_2 = tau - tau_1;
  float CN = PI / (1.4*PI*tau_1 + 1.2*tau_1 + 0.3*PI*tau_2);

  if (t >=0.0 & t < tau_1) {
    f = CN * (0.7 - 0.7 * cos(PI*t/tau_1) + 0.6 * sin(0.5*PI*t/tau_1));
  } else if (t < 2.0 * tau_1) {
    f = CN * (1.0 - 0.7 * cos(PI*t/tau_1) + 0.3 * cos(PI*(t-tau_1)/tau_2));
  } else if (t < tau) {
    f = CN * (0.3 + 0.3 * cos(PI*(t-tau_1) / tau_2));
  } else {
    f = 0.0;
  }

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

/*
 *
 */

int
src_muDA_to_moment(float strike, float dip, float rake, float mu, float D, float A,
          float *mxx, float *myy, float *mzz, float *myz, float *mxz, float *mxy)
{
  int ierr = 0;

  // convert to moment tensor
  float temp_moment[6];
  angle2moment(strike, dip, rake, temp_moment);

  float M0 = mu * D * A;
  
  *mxx = M0*temp_moment[0];
  *myy = M0*temp_moment[1];
  *mzz = M0*temp_moment[2];
  *myz = M0*temp_moment[3];
  *mxz = M0*temp_moment[4];
  *mxy = M0*temp_moment[5];

  return ierr;
}

/*
 *
 */

int
src_print(src_t *src, int verbose)
{
  int ierr = 0;

  // evtnm has a terminal char at end
  fprintf(stdout,"-- evtnm=%s\n", src->evtnm);

  fprintf(stdout,"-- total_number=%d\n", src->total_number);
  fprintf(stdout,"-- force_actived=%d, moment_actived=%d\n",
          src->force_actived, src->moment_actived);
  fprintf(stdout,"-- max_nt=%d,max_stage=%d,max_ext=%d\n",
          src->max_nt,src->max_stage,src->max_ext);
  
  // only print for large verbose
  if (verbose > 99)
  {
    for (int is=0; is<src->total_number; is++)
    {
      fprintf(stdout,"--- is=%d, si=%d,sj=%d,sk=%d,ext_num=%d,it_begin=%d,it_end=%d\n",
            is,src->si[is],src->sj[is],src->sk[is],src->ext_num[is],
            src->it_begin[is],src->it_end[is]);
      // should not print time series for normal usage
      if (verbose > 999)
      {
        for (int it = src->it_begin[is]; it <= src->it_end[is]; it++)
        {
          int it_to_it_start = it - src->it_begin[is];
          // print 0 stage
          size_t iptr = is * src->max_nt * src->max_stage
                        + it_to_it_start * src->max_stage;
          fprintf(stdout, "---- it=%d",it);
          if (src->force_actived==1) {
            fprintf(stdout, ",fx=%g,fy=%g,fz=%g",
                      src->Fx[iptr], src->Fy[iptr], src->Fz[iptr]);
          }
          if (src->moment_actived==1) {
            fprintf(stdout, ",Mxx=%g,Myy=%g,Mzz=%g,Myz=%g,Mxz=%g,Mxy=%g",
                      src->Mxx[iptr], src->Myy[iptr], src->Mzz[iptr],
                      src->Myz[iptr], src->Mxz[iptr], src->Mxy[iptr]);
          }
          fprintf(stdout, "\n");
        }
      }
    }
  }

  return ierr;
}

