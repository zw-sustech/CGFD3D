#ifndef WF_EL_1ST_H
#define WF_EL_1ST_H

#include "gd_info.h"

/*************************************************
 * structure
 *************************************************/

/*
 * wavefield structure
 */

// wavefield variables elastic 1st eqn: vel + stress
typedef struct {
  float *v5d; // allocated var

  int n1, n2, n3, n4, n5;
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

  // sequential index 0-based
  size_t Vx_seq;
  size_t Vy_seq;
  size_t Vz_seq;
  size_t Txx_seq;
  size_t Tyy_seq;
  size_t Tzz_seq;
  size_t Tyz_seq;
  size_t Txz_seq;
  size_t Txy_seq;

} wav_t;

struct fd_vel_t
{
  float *Vx;
  float *Vy;
  float *Vz;
};

struct fd_stress_t
{
  float *Txx;
  float *Tyy;
  float *Tzz;
};

struct var5d_t
{
  int n1, n2, n3, n4, n5;
  int nx, ny, nz, ncmp, nlevel;
  float *v5d;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;
  size_t siz_ilevel;

  size_t *cmp_pos;
  char  **cmp_name;
};

/*************************************************
 * function prototype
 *************************************************/

int 
wav_init(gdinfo_t *gdinfo,
               wav_t *V,
               int number_of_levels);

int
wav_check_value(float *restrict w, wav_t *wav);

int
wav_zero_edge(gdinfo_t *gdinfo, wav_t *wav,
                                  float *restrict w4d);

int 
wav_ac_init(gdinfo_t *gdinfo,
               wav_t *V,
               int number_of_levels);

#endif
