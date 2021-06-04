/*
 * source term related processing
 */

// todo:

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "src_funcs.h"

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
                           int   *source_gridindex,
                           float *source_coords,
                           float *force_vector,
                           float *moment_tensor,
                           char  *wavelet_name,
                           float *wavelet_coefs,
                           float wavelet_tstart,
                           float wavelet_tend,
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
  int  it_end   = (int) ((wavelet_tend + 0.5) / dt);
  int  nt_total_wavelet = it_end - it_begin + 1;
  float *wavelet_values = NULL;

  // default local index and relative shift to -1 and 0
  int si=-1; int sj=-1; int sk=-1;
  float sx_inc = 0.0; float sy_inc = 0.0; float sz_inc = 0.0;

  // locate into this thread
  if (sgpi < 0 || sgpj < 0 || sgpk < 0) {
    // use coord to locate local index
    fprintf(stdout,"locate source by coord ...\n"); 
    fprintf(stdout,"   not implemented yet\n"); 
    fflush(stdout);
  }
  else // use grid index
  {
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
        } else {
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

          int iptr = i + j * siz_line + k * siz_slice;
          moment_ext_indx[iptr_s] = iptr;
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

          int iptr = i + j * siz_line + k * siz_slice;
          force_ext_indx[iptr_s] = iptr;
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
}

/*
 * Locate 
 */
void 
eighth_locate(float *location, 
              float *vx, float *vy, float *vz, 
              float sx_inc, int sy_inc, int sz_inc)
{
  float x0 = location[0];
  float y0 = location[1];
  float z0 = location[2];
  
  float x[125] ;
  float y[125] ;
  float z[125] ;
 

  for( int i=0; i<5; i=i+2 )
  {
    for( int j=0; j<5; j=j+2) 
    {
     for ( int k=0; k<5; k=k+2)
      {
       int indx_01  = k * 25 + j * 5 + i ;
       int indx_02  = ((int)(k / 2)) * 9 + ((int)(j / 2)) * 3 + ((int)(i / 2));

       
      
      }

    }

  }



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
    int ipos = moment_info[M_SRC_INFO_SEQ_POS];
    int it1  = moment_info[M_SRC_INFO_SEQ_ITBEG];
    int it2  = moment_info[M_SRC_INFO_SEQ_ITEND];
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


void 
angle2moment(float strike, float dip, float rake, float* source_moment_tensor)
{
  float strike_pi,dip_pi,rake_pi; 
  float M11,M22,M33,M12,M13,M23;

  dip_pi    = dip    / 180.0 * PI; 
  strike_pi = strike / 180.0 * PI;
  rake_pi   = rake   / 180.0 * PI;

 // in Aki and Richard's
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
 
  source_moment_tensor[0] = M11 ; 
  source_moment_tensor[1] = M22 ;   
  source_moment_tensor[2] = M33 ;
  source_moment_tensor[3] = M12 ;  
  source_moment_tensor[4] = M13 ;
  source_moment_tensor[5] = M23 ;  
 
}

/* structure for keep vertex of cubic*/

/*struct CubicPt {

float xcoord; 
float ycoord; 
float zcoord;

int xindx; 
int yindx;
int zindx;

} ;*/


/* search location of source
 *
 *
 */
void
LocaSrc(float sx, float sy, float sz,
        int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
        size_t siz_line, size_t siz_slice, size_t siz_volume, 
        float *restrict c3d,
        size_t *restrict c3d_pos,
        struct 
        )
{

  int indx, NearIndx;
  float Dist = 0.0 ; 
  float DistInt = 0.0 ;

  float *restrict xcoord = c3d + c3d_pos[0];
  float *restrict ycoord = c3d + c3d_pos[1];
  float *restrict zcoord = c3d + c3d_pos[2];
   
  /* search minimum distance  */ 
  DistInt =  (sx - xcoord[0]) * (sx - xcoord[0])
           + (sy - ycoord[0]) * (sy - ycoord[0])
           + (sz - zcoord[0]) * (sz - zcoord[0]);

  for(int k=nk1; k<nk2; k++)
   {
     for(int j=nj1; j<nj2; j++)
      {
        for(int i=ni1; i<ni2; i++)
         {
          indx = i + j * siz_line + k * siz_slice ;
         
          Dist =  (sx - xcoord[indx]) * (sx - xcoord[indx])
                + (sy - ycoord[indx]) * (sy - ycoord[indx])
                + (sz - zcoord[indx]) * (sz - zcoord[indx]);
                  
          if (Dist < DistInt)
           {
             /* Keep indx for minimum distance */
             NearSi = i ;  NearSj = j ; NearSk = k ;
             NearIndx = Indx;
           }
         }
      }
   }

  float NearPtX = xcoord[NearIndx]; 
  float NearPtY = ycoord[NearIndx]; 
  float NearPtZ = zcoord[NearIndx]; 

 /* eight situtations 
  * Always save the nearest point at the first
  * */ 
 
  if ( sx > NearPtX && sy > NearPtY && sz > NearPtZ  )
     {
       

     }
  else if ( sx > NearPtX && sy > NearPtY && sz < NearPtZ)
     {
     
     }
  else if ( sx > NearPtX && sy < NearPtY && sz > NearPtZ) 
     {
     
     }
  else if ( sx > NearPtX && sy < NearPtY && sz < NearPtZ) 
     {
     
     }
  else if ( sx < NearPtX && sy > NearPtY && sz > NearPtZ  )
     {
     
     }
  else if ( sx < NearPtX && sy > NearPtY && sz < NearPtZ)
     {
     
     }
  else if ( sx < NearPtX && sy < NearPtY && sz > NearPtZ) 
     {
     
     }
  else if ( sx < NearPtX && sy < NearPtY && sz < NearPtZ) 
     {
     
     }

}



/* Coordinate mapping using  Inverse Distance Weight
 * sx, sy, sz: physical cooedinate for source  
 * struct cubcoord : structure cubic point coordinate coordx, coordy, coordz (float)
 * struct cubindx  : structure cubic point indx       indxx,  indxy,  indxz  (int)
 *
 * 
 */
float 
CoorMap(float sx, float sy, float sz, 
        struct cubcoord *Pcoord,
        struct cubindx  *Pindx)

{

 float Dist[8];
 float SUM;
 float Weight[8];
 float Mapsi, Mapsj, Mapsk;
 float shift;
 int si, sj, sk; 
 
 SUM = 0.0 ;
 
 for (int i=0; i<8; i++)
  {
   Dist[i] = sqrt ((sx - Pcoord[i].coordx) * (sx - Pcoord[i].coordx)
                 + (sy - Pcoord[i].coordy) * (sy - Pcoord[i].coordy)
                 + (sz - Pcoord[i].coordz) * (sz - Pcoord[i].coordz));

   SUM += 1.0/Dist[i]; 
  }

 for (int i=0; i<8; i++)
  {

   Weight[i] = 1.0 / Dist[i] / SUM ;
 
   Mapsi += Weight[i] * Pindx[i].indxx;
   Mapsj += Weight[i] * Pindx[i].indxy; 
   Mapsk += Weight[i] * Pindx[i].indxz;  
  
  }


}
