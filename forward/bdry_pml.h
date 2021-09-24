#ifndef BDRY_PML_H
#define BDRY_PML_H

#include "constants.h"
#include "gd_info.h"
#include "gd_t.h"
#include "wav_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
  float *var;
  int nx, ny, nz, ncmp, nlevel;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;
  size_t siz_ilevel;

  size_t *cmp_pos;
  char  **cmp_name;
  size_t *level_pos;

  size_t Vx_pos;
  size_t Vy_pos;
  size_t Vz_pos;
  size_t Txx_pos;
  size_t Tyy_pos;
  size_t Tzz_pos;
  size_t Tyz_pos;
  size_t Txz_pos;
  size_t Txy_pos;

  // for rk scheme
  float *cur;
  float *pre;
  float *rhs;
  float *end;
  float *tmp;

} bdrypml_auxvar_t;

typedef struct
{
  int is_enable; //
  int is_at_sides[CONST_NDIM][2];
  int num_of_layers[CONST_NDIM][2]; //

  int ni1[CONST_NDIM][2];
  int ni2[CONST_NDIM][2];
  int nj1[CONST_NDIM][2];
  int nj2[CONST_NDIM][2];
  int nk1[CONST_NDIM][2];
  int nk2[CONST_NDIM][2];

  float *A[CONST_NDIM][2]; // dim, side, length
  float *B[CONST_NDIM][2]; // dim, side, length
  float *D[CONST_NDIM][2]; // dim, side, length

  // for middile point of staggered grid
  float *Am[CONST_NDIM][2]; // dim, side, length
  float *Bm[CONST_NDIM][2]; // dim, side, length
  float *Dm[CONST_NDIM][2]; // dim, side, length

  bdrypml_auxvar_t auxvar[CONST_NDIM][2];

} bdrypml_t;

/*************************************************
 * function prototype
 *************************************************/

float
bdry_pml_cal_R(float N);

float
bdry_pml_cal_dmax(float L, float Vp, float Rpp);

float
bdry_pml_cal_amax(float fc);

float
bdry_pml_cal_d(float x, float L, float dmax);

float
bdry_pml_cal_a(float x, float L, float amax);

float
bdry_pml_cal_b(float x, float L, float bmax);

void
bdry_pml_set(gdinfo_t *gdinfo,
             gd_t *gd,
             wav_t *wav,
             bdrypml_t *bdrypml,
             int   *neighid, 
             int   in_is_sides[][2],
             int   in_num_layers[][2],
             float in_alpha_max[][2], //
             float in_beta_max[][2], //
             float in_velocity[][2], //
             int verbose);

void
bdry_pml_set_stg(gdinfo_t *gdinfo,
                 gd_t *gd,
                 wav_t *wav,
                 bdrypml_t *bdrypml,
                 int   *neighid, 
                 int   in_is_sides[][2],
                 int   in_num_layers[][2],
                 float in_alpha_max[][2], //
                 float in_beta_max[][2], //
                 float in_velocity[][2], //
                 int verbose);

// alloc auxvar
void
bdry_pml_auxvar_init(int nx, int ny, int nz, 
                     wav_t *wav,
                     bdrypml_auxvar_t *auxvar,
                     const int verbose);

int
bdry_pml_cal_len_dh(gd_t *gd, 
                    int abs_ni1, int abs_ni2,
                    int abs_nj1, int abs_nj2,
                    int abs_nk1, int abs_nk2,
                    int idim,
                    float *avg_L, float *avg_dh);

#endif
