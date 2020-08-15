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
  float stf_fc = 1.0;
  float stf_dur = 4.0;

  int nforce = 0;
  
  int nmoment = 1;
  int *moment_info = (int *)fdlib_mem_calloc_1d_int(nmoment*6, 1, "src_gen_test");

  moment_info[0] = 50; //si
  moment_info[1] = 50; //sj
  moment_info[2] = 50; //sk
  moment_info[3] = 0;
  moment_info[4] = 0;
  moment_info[5] = (int) stf_dur / dt;

  int ipos = moment_info[3];
  int it1 = moment_info[4];
  int it2 = moment_info[5];
  int nt_moment = it2 - it1 + 1;

  float *moment_ten_rate = (float *)fdlib_mem_calloc_1d_float(
                                      nt_moment*6*num_of_stages,0.0,"src_gen_test");

  float rk_a[4] = {0.0,0.5,0.5,1.0};

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
        float t = it * dt + t0 + rk_a[istage] * dt;
        float stf_val = fun_ricker(t, stf_fc, 2.0);
        moment_ten_rate[iptr] = stf_val * Mij;
      }
    }
  }

  *num_of_force = nforce;
  *num_of_moment = nmoment;
  *p_moment_info = moment_info;
  *p_moment_ten_rate = moment_ten_rate;
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
    int ipos = force_info[3];
    int it1  = force_info[4];
    int it2  = force_info[5];
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
    int ipos = moment_info[3];
    int it1 = moment_info[4];
    int it2 = moment_info[5];
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
}
