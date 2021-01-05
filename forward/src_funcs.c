/*
 * source term related processing
 */

// todo:

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "src_funcs.h"

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
 * test
 */

void
src_gen_test(
    float t0,
    float dt,
    int   nt_total,
    int   num_of_stages,
    int  *num_of_force,
    int **restrict p_force_info,
    float  **restrict p_force_vec_stf,
    int  *num_of_moment,
    int **restrict p_moment_info,
    float  **restrict p_moment_ten_rate,
    int verbose)
{
  float stf_fc = 2.0;
  float stf_dur = 1.0;
  //float stf_fc = 1.0;
  //float stf_dur = 2.0;

  int nforce = 0;
  
  int nmoment = 0;
  int *moment_info = (int *)fdlib_mem_calloc_1d_int(nmoment*M_SRC_INFO_NVAL, 1, "src_gen_test");

//  moment_info[M_SRC_INFO_SEQ_SI   ] = 52; //si
 // moment_info[M_SRC_INFO_SEQ_SJ   ] = 22; //sj
 // moment_info[M_SRC_INFO_SEQ_SK   ] = 22; //sk
  moment_info[M_SRC_INFO_SEQ_SJ   ] = 52; //sj
  moment_info[M_SRC_INFO_SEQ_SK   ] = 58; //sk

  moment_info[M_SRC_INFO_SEQ_SK   ] = 32; //sk
  moment_info[M_SRC_INFO_SEQ_POS  ] = 0;
  moment_info[M_SRC_INFO_SEQ_ITBEG] = 0;
  moment_info[M_SRC_INFO_SEQ_ITEND] = (int) (stf_dur / dt);

  int ipos = moment_info[M_SRC_INFO_SEQ_POS];
  int it1  = moment_info[M_SRC_INFO_SEQ_ITBEG];
  int it2  = moment_info[M_SRC_INFO_SEQ_ITEND];
  int nt_moment = it2 - it1 + 1;

  float *moment_ten_rate = (float *)fdlib_mem_calloc_1d_float(
                                      nt_moment*6*num_of_stages,0.0,"src_gen_test");

  float rk_tinc[4] = {0.0,0.5,0.5,1.0};

  for (int icmp=0; icmp<6; icmp++)
  {
    float Mij=0.0;
    // explosive
    if (icmp<3) Mij=1.0e9;

    for (int it=it1; it<=it2; it++)
    {
      int it_to_it1 = (it - it1);
      for (int istage=0; istage<num_of_stages; istage++)
      {
        int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_moment,num_of_stages);
        float t = it * dt + t0 + rk_tinc[istage] * dt;
        float stf_val = fun_ricker(t, stf_fc, stf_dur/2.0);
        moment_ten_rate[iptr] = stf_val * Mij;
      }
    }
  }

  *num_of_force = nforce;
  *num_of_moment = nmoment;
  *p_moment_info = moment_info;
  *p_moment_ten_rate = moment_ten_rate;
}

void
src_gen_test_gauss(
    size_t siz_line,
    size_t siz_slice,
    float t0,
    float dt,
    int   nt_total,
    int   num_of_stages,
    int  *num_of_force,
    int **restrict p_force_info,
    float  **restrict p_force_vec_stf,
    int  *num_of_moment,
    int **restrict p_moment_info,
    float  **restrict p_moment_ten_rate,
    int **restrict p_moment_ext_indx,
    float  **restrict p_moment_ext_coef,
    int verbose)
{
  float stf_fc = 2.0;
  float stf_dur = 1.0;
  //float stf_fc = 1.0;
  //float stf_dur = 2.0;

  int nforce = 0;
  int nmoment = 1;

  int siz_ext = 7 * 7 * 7;
  int si = 55;
  int sj = 55;
  int sk = 32;
  //int sk = 58; // for free surface
  float sx_inc = 0.0;
  float sy_inc = 0.0;
  float sz_inc = 0.0;
  //float sx_inc = -0.3;
  //float sy_inc = 0.4;
  //float sz_inc = 0.2;

  // info
  int *moment_info = (int *)fdlib_mem_calloc_1d_int(nmoment*M_SRC_INFO_NVAL, 1, "src_gen_test");

  moment_info[M_SRC_INFO_SEQ_SI   ] = si;
  moment_info[M_SRC_INFO_SEQ_SJ   ] = sj;
  moment_info[M_SRC_INFO_SEQ_SK   ] = sk;
  moment_info[M_SRC_INFO_SEQ_POS  ] = 0;
  moment_info[M_SRC_INFO_SEQ_ITBEG] = 0;
  moment_info[M_SRC_INFO_SEQ_ITEND] = (int) (stf_dur / dt);
  moment_info[M_SRC_INFO_SEQ_NEXT ] = siz_ext;

  int ipos = moment_info[M_SRC_INFO_SEQ_POS];
  int it1  = moment_info[M_SRC_INFO_SEQ_ITBEG];
  int it2  = moment_info[M_SRC_INFO_SEQ_ITEND];
  int nt_moment = it2 - it1 + 1;

  // stf
  float *moment_ten_rate = (float *)fdlib_mem_calloc_1d_float(
                                      nt_moment*6*num_of_stages,0.0,"src_gen_test");

  float rk_tinc[4] = {0.0,0.5,0.5,1.0};

  for (int icmp=0; icmp<6; icmp++)
  {
    float Mij=0.0;
    // explosive
    if (icmp<3) Mij=1.0e9;

    for (int it=it1; it<=it2; it++)
    {
      int it_to_it1 = (it - it1);
      for (int istage=0; istage<num_of_stages; istage++)
      {
        int iptr = M_SRC_IND(icmp,it_to_it1,istage,nt_moment,num_of_stages);
        float t = it * dt + t0 + rk_tinc[istage] * dt;
        float stf_val = fun_ricker(t, stf_fc, stf_dur/2.0);
        moment_ten_rate[iptr] = stf_val * Mij;
      }
    }
  }

  // ext
  float *moment_ext_coef = (float *)malloc(siz_ext * sizeof(float));
  int   *moment_ext_indx = (int   *)malloc(siz_ext * sizeof(int  ));

  cal_norm_delt3d(moment_ext_coef, sx_inc, sy_inc, sz_inc, 1.5, 1.5, 1.5, 3);

  size_t iptr_s = 0;
  for (int k=sk-3; k<=sk+3; k++) {
    for (int j=sj-3; j<=sj+3; j++) {
      for (int i=si-3; i<=si+3; i++) {
        int iptr = i + j * siz_line + k * siz_slice;
        moment_ext_indx[iptr_s] = iptr;
        iptr_s++;
      }
    }
  }

  *num_of_force = nforce;
  *num_of_moment = nmoment;
  *p_moment_info = moment_info;
  *p_moment_ten_rate = moment_ten_rate;
  *p_moment_ext_indx = moment_ext_indx;
  *p_moment_ext_coef = moment_ext_coef;
}

// ricker and it deriv.
float 
fun_ricker(float t, float fc, float t0)
{
    float pi = acos(-1.0);
    float f0 = sqrtf(pi)/2.0;
    float u = (t-t0)*2.0*pi*fc;
    float v = (u*u/4-0.5)*exp(-u*u/4)*f0;

    return v;
}

float
fun_gauss(float t, float a, float t0)
{
    float f;
    f = exp(-(t-t0)*(t-t0)/(a*a))/(sqrtf(PI)*a);
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
    int ipos = force_info[M_SRC_INFO_SEQ_POS];
    int it1  = force_info[M_SRC_INFO_SEQ_ITBEG];
    int it2  = force_info[M_SRC_INFO_SEQ_ITEND];
    int nt_force = it2 - it1 + 1;

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
    int ipos = moment_info[M_SRC_INFO_SEQ_POS];
    int it1  = moment_info[M_SRC_INFO_SEQ_ITBEG];
    int it2  = moment_info[M_SRC_INFO_SEQ_ITEND];
    int nt_moment = it2 - it1 + 1;

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
