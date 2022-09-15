/*********************************************************************
 * wavefield for 3d elastic 1st-order equations
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "constants.h"
#include "fdlib_mem.h"
#include "wav_t.h"

int 
wav_init(gdinfo_t *gdinfo,
               wav_t *V,
               int number_of_levels)
{
  int ierr = 0;

  V->nx   = gdinfo->nx;
  V->ny   = gdinfo->ny;
  V->nz   = gdinfo->nz;
  V->ncmp = 9;
  V->nlevel = number_of_levels;

  V->siz_iy   = V->nx;
  V->siz_iz   = V->nx * V->ny;
  V->siz_icmp = V->nx * V->ny * V->nz;
  V->siz_ilevel = V->siz_icmp * V->ncmp;

  // vars
  // 3 Vi, 6 Tij, 4 rk stages
  V->v5d = (float *) fdlib_mem_calloc_1d_float(V->siz_ilevel * V->nlevel,
                        0.0, "v5d, wf_el3d_1st");
  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                      V->ncmp, 0, "w3d_pos, wf_el3d_1st");
  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(
                      V->ncmp, CONST_MAX_STRLEN, "w3d_name, wf_el3d_1st");
  
  // set value
  for (int icmp=0; icmp < V->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * V->siz_icmp;
  }

  // set values
  int icmp = 0;

  /*
   * 0-3: Vx,Vy,Vz
   * 4-9: Txx,Tyy,Tzz,Txz,Tyz,Txy
   */

  sprintf(cmp_name[icmp],"%s","Vx");
  V->Vx_pos = cmp_pos[icmp];
  V->Vx_seq = 0;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vy");
  V->Vy_pos = cmp_pos[icmp];
  V->Vy_seq = 1;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vz");
  V->Vz_pos = cmp_pos[icmp];
  V->Vz_seq = 2;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txx");
  V->Txx_pos = cmp_pos[icmp];
  V->Txx_seq = 3;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Tyy");
  V->Tyy_pos = cmp_pos[icmp];
  V->Tyy_seq = 4;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Tzz");
  V->Tzz_pos = cmp_pos[icmp];
  V->Tzz_seq = 5;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Tyz");
  V->Tyz_pos = cmp_pos[icmp];
  V->Tyz_seq = 6;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txz");
  V->Txz_pos = cmp_pos[icmp];
  V->Txz_seq = 7;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Txy");
  V->Txy_pos = cmp_pos[icmp];
  V->Txy_seq = 8;
  icmp += 1;

  // set pointer
  V->cmp_pos  = cmp_pos;
  V->cmp_name = cmp_name;

  return ierr;
}

int 
wav_ac_init(gdinfo_t *gdinfo,
               wav_t *V,
               int number_of_levels)
{
  int ierr = 0;

  // Vx,Vy,Vz,P
  V->ncmp = 4;

  V->nx   = gdinfo->nx;
  V->ny   = gdinfo->ny;
  V->nz   = gdinfo->nz;
  V->nlevel = number_of_levels;

  V->siz_iy   = V->nx;
  V->siz_iz   = V->nx * V->ny;
  V->siz_icmp = V->nx * V->ny * V->nz;
  V->siz_ilevel = V->siz_icmp * V->ncmp;

  // vars
  // 3 Vi, 6 Tij, 4 rk stages
  V->v5d = (float *) fdlib_mem_calloc_1d_float(V->siz_ilevel * V->nlevel,
                        0.0, "v5d, wf_ac3d_1st");
  // position of each var
  size_t *cmp_pos = (size_t *) fdlib_mem_calloc_1d_sizet(
                      V->ncmp, 0, "w3d_pos, wf_ac3d_1st");
  // name of each var
  char **cmp_name = (char **) fdlib_mem_malloc_2l_char(
                      V->ncmp, CONST_MAX_STRLEN, "w3d_name, wf_ac3d_1st");
  
  // set value
  for (int icmp=0; icmp < V->ncmp; icmp++)
  {
    cmp_pos[icmp] = icmp * V->siz_icmp;
  }

  // set values
  int icmp = 0;

  /*
   * 0-3: Vx,Vy,Vz
   * 4: P
   */

  sprintf(cmp_name[icmp],"%s","Vx");
  V->Vx_pos = cmp_pos[icmp];
  V->Vx_seq = 0;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vy");
  V->Vy_pos = cmp_pos[icmp];
  V->Vy_seq = 1;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Vz");
  V->Vz_pos = cmp_pos[icmp];
  V->Vz_seq = 2;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","P");
  V->Txx_pos = cmp_pos[icmp];
  V->Txx_seq = 3;
  icmp += 1;

  // set pointer
  V->cmp_pos  = cmp_pos;
  V->cmp_name = cmp_name;

  return ierr;
}

int
wav_check_value(float *restrict w, wav_t *wav)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *ptr = w + icmp * wav->siz_icmp;
    for (size_t iptr=0; iptr < wav->siz_icmp; iptr++)
    {
      if (ptr[iptr] != ptr[iptr])
      {
        fprintf(stderr, "ERROR: NaN occurs at iptr=%d icmp=%d\n", iptr, icmp);
        fflush(stderr);
        exit(-1);
      }
    }
  }

  return ierr;
}

int
wav_zero_edge(gdinfo_t *gdinfo, wav_t *wav,
                                  float *restrict w4d)
{
  int ierr = 0;

  for (int icmp=0; icmp < wav->ncmp; icmp++)
  {
    float *restrict var = w4d + wav->cmp_pos[icmp];

    // z1
    for (int k=0; k < gdinfo->nk1; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j=0; j < gdinfo->ny; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // z2
    for (int k=gdinfo->nk2+1; k < gdinfo->nz; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j=0; j < gdinfo->ny; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // y1
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j=0; j < gdinfo->nj1; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // y2
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j = gdinfo->nj2+1; j < gdinfo->ny; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    }

    // x1
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j = gdinfo->nj1; j <= gdinfo->nj2; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i=0; i < gdinfo->ni1; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    } 

    // x2
    for (int k = gdinfo->nk1; k <= gdinfo->nk2; k++)
    {
      size_t iptr_k = k * gdinfo->siz_iz;
      for (int j = gdinfo->nj1; j <= gdinfo->nj2; j++)
      {
        size_t iptr_j = iptr_k + j * gdinfo->siz_iy;
        for (int i = gdinfo->ni2+1; i < gdinfo->nx; i++)
        {
          size_t iptr = iptr_j + i;
          var[iptr] = 0.0; 
        }
      }
    } 

  } // icmp

  return ierr;
}
int
PG_calcu(float *w_end, float *w_pre, gdinfo_t *gdinfo, float *PG, float *Dis_accu, float dt)
{
  //Dis_accu is displacement accumulation.
  int ni1 = gdinfo->ni1;
  int nj1 = gdinfo->nj1;
  int nk1 = gdinfo->nk1;
  int ni2 = gdinfo->ni2;
  int nj2 = gdinfo->nj2;
  int nk2 = gdinfo->nk2;
  int siz_line = gdinfo->siz_line;
  int siz_slice  = gdinfo->siz_slice;
  int siz_volume = gdinfo->siz_volume;
  // 0-2 Vx1,Vy1,Vz1  it+1 moment  V
  // 0-2 Vx0,Vy0,Vz0  it moment V
  float *Vx1 = w_end + 0*siz_volume;
  float *Vy1 = w_end + 1*siz_volume;
  float *Vz1 = w_end + 2*siz_volume;
  float *Vx0 = w_pre + 0*siz_volume;
  float *Vy0 = w_pre + 1*siz_volume;
  float *Vz0 = w_pre + 2*siz_volume;
  float *PGVx = PG + 0*siz_slice;
  float *PGVy = PG + 1*siz_slice;
  float *PGVz = PG + 2*siz_slice;
  float *PGAx = PG + 3*siz_slice;
  float *PGAy = PG + 4*siz_slice;
  float *PGAz = PG + 5*siz_slice;
  float *PGDx = PG + 6*siz_slice;
  float *PGDy = PG + 7*siz_slice;
  float *PGDz = PG + 8*siz_slice;
  float *D_x = Dis_accu + 0*siz_slice;
  float *D_y = Dis_accu + 1*siz_slice;
  float *D_z = Dis_accu + 2*siz_slice;
  int iptr,iptr1;
  for( int j = nj1; j<=nj2; j++){
    for( int i = ni1; i<=ni2; i++){
      iptr = i + j * siz_line + nk2 * siz_slice;
      iptr1 = i + j * siz_line;
      float Ax, Ay, Az;
      Ax = fabs((Vx1[iptr]-Vx0[iptr])/dt);
      Ay = fabs((Vy1[iptr]-Vy0[iptr])/dt);
      Az = fabs((Vz1[iptr]-Vz0[iptr])/dt);
      D_x[iptr1] += 0.5*(Vx1[iptr]+Vx0[iptr])*dt; 
      D_y[iptr1] += 0.5*(Vy1[iptr]+Vy0[iptr])*dt; 
      D_z[iptr1] += 0.5*(Vz1[iptr]+Vz0[iptr])*dt; 
      if(PGVx[iptr1]<fabs(Vx1[iptr])) PGVx[iptr1]=fabs(Vx1[iptr]);
      if(PGVy[iptr1]<fabs(Vy1[iptr])) PGVy[iptr1]=fabs(Vy1[iptr]);
      if(PGVz[iptr1]<fabs(Vz1[iptr])) PGVz[iptr1]=fabs(Vz1[iptr]);
      if(PGAx[iptr1]<Ax) PGAx[iptr1]=Ax;
      if(PGAy[iptr1]<Ay) PGAy[iptr1]=Ay;
      if(PGAz[iptr1]<Az) PGAz[iptr1]=Az;
      if(PGDx[iptr1]<fabs(D_x[iptr1])) PGDx[iptr1] = fabs(D_x[iptr1]);
      if(PGDy[iptr1]<fabs(D_y[iptr1])) PGDy[iptr1] = fabs(D_y[iptr1]);
      if(PGDz[iptr1]<fabs(D_z[iptr1])) PGDz[iptr1] = fabs(D_z[iptr1]);
      }
    }
return 0;
}
