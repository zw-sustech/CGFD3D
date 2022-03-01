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
#include "io_funcs.h"
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
 * read .src file and convert into internal structure
 */

int
src_read_locate_file(gdinfo_t *gdinfo,
                     gd_t     *gd,
                     src_t    *src,
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
  FILE *fp =NULL;
  char str[500];
  //char *in_source_name = (char *) malloc(sizeof(char)*500);

  // numbe of source, could be force and/or moment
  int in_num_source;
  // input location is grid index (0) or coordinate (1)
  int in_location_meaning;
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
    sscanf(str,"%s",str);
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
    sscanf(str,"%d %d",&in_location_meaning, &in_3coord_meaning);
  }

  //
  // loop each source to get locate and index
  //

  int num_of_src_here = 0;

  float sx, sy, sz;
  int   si,sj,sk;
  int   si_glob,sj_glob,sk_glob;
  float sx_inc, sy_inc, sz_inc;

  float **all_coords = fdlib_mem_calloc_2l_float(in_num_source,CONST_NDIM,0.0,"all_coords");
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

      if (in_location_meaning == 1) // physical coord
      {
        // convert to global index
        //    todo: check if out of computational region
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
        all_index[is][0] = si_glob;
        all_index[is][1] = sj_glob;
        all_index[is][2] = sk_glob;
        all_inc  [is][0] = sx_inc;
        all_inc  [is][1] = sy_inc;
        all_inc  [is][2] = sz_inc;

        //-- to notice user the progress using screen output for large input
        if (myid == 0 && (is % 1000 ==0) && verbose>99) {
          fprintf(stdout,"-- loc %d-th src index, finish %f\%\n",
                      is, (float)(is+1)/in_num_source*100.0);
          fflush(stdout);
        }
      }
      else // computational coordinate or grid index
      {
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
    if (src_glob_ext_ishere(si_glob,sj_glob,sk_glob,npoint_half_ext,gdinfo)==1)
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
    } else if (in_cmp_type == 2 || in_cmp_type == 3) {
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
            src_muDA_to_moment(mxx,myy,mzz,myz,mxz,mxy,
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
              src_muDA_to_moment(m11[it],m22[it],m33[it],m23[it],m13[it],m12[it],
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
      si = gd_info_ind_glphy2lcext_i(all_index[is][0], gdinfo);
      sj = gd_info_ind_glphy2lcext_j(all_index[is][1], gdinfo);
      sk = gd_info_ind_glphy2lcext_k(all_index[is][2], gdinfo);
  
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

            float stf_val = src_cal_wavelet(t,wavelet_name,wavelet_coefs);
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

  fdlib_mem_free_2l_float(all_coords, in_num_source, "free_all_coords");
  fdlib_mem_free_2l_float(all_inc, in_num_source, "free_all_inc");
  fdlib_mem_free_2l_int  (all_index, in_num_source, "free_all_index");

  return ierr;
}

/*
 * get stf value at a given t
 */

float
src_cal_wavelet(float t, char *wavelet_name, float *wavelet_coefs)
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
  fprintf(stdout,"-- evtnm=%s", src->evtnm);

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

