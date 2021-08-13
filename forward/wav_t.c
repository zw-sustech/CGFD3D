/*********************************************************************
 * wavefield for 3d elastic 1st-order equations
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

  sprintf(cmp_name[icmp],"%s","Txz");
  V->Txz_pos = cmp_pos[icmp];
  V->Txz_seq = 6;
  icmp += 1;

  sprintf(cmp_name[icmp],"%s","Tyz");
  V->Tyz_pos = cmp_pos[icmp];
  V->Tyz_seq = 7;
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

//int
//wav_set_stag_loc(wav_t *wav)
//{
//  wav->Vx_stg_loc = (wav_stag_loc_t **)malloc(wav->nlevel * sizeof(wav_stag_loc_t*));
//  for (int n=0; n < wav->nlevel; n++) {
//    wav->Vx_stg_loc[n] = (wav_stag_loc_t *)malloc(CONST_NDIM * sizeof(wav_stag_loc_t));
//  }
//
//  Vx_stg_loc[];
//}

int
wav_set_cell_loc_stg(wav_t *wav)
{
  // Vx x-half
  wav->Vx_cell_loc[0] = 1;
  wav->Vx_cell_loc[1] = 0;
  wav->Vx_cell_loc[2] = 0;

  // Vy y-half
  wav->Vy_cell_loc[0] = 0;
  wav->Vy_cell_loc[1] = 1;
  wav->Vy_cell_loc[2] = 0;

  // Vz z-half
  wav->Vz_cell_loc[0] = 0;
  wav->Vz_cell_loc[1] = 0;
  wav->Vz_cell_loc[2] = 1;

  // Txx at point
  wav->Txx_cell_loc[0] = 0;
  wav->Txx_cell_loc[1] = 0;
  wav->Txx_cell_loc[2] = 0;

  // Tyy at point
  wav->Tyy_cell_loc[0] = 0;
  wav->Tyy_cell_loc[1] = 0;
  wav->Tyy_cell_loc[2] = 0;

  // Tzz at point
  wav->Tzz_cell_loc[0] = 0;
  wav->Tzz_cell_loc[1] = 0;
  wav->Tzz_cell_loc[2] = 0;
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
