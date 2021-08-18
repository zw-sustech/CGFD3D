#ifndef SV_EQ1ST_CART_STG_AC_ISO_H
#define SV_EQ1ST_CART_STG_AC_ISO_H

#include "fd_t.h"
#include "gd_info.h"
#include "mympi_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_free.h"
#include "bdry_pml.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

void
sv_eq1st_cart_stg_ac_iso_allstep(
  fdstg_t    *fd,
  gdinfo_t   *gdinfo,
  gd_t       *gdcart,
  md_t       *md,
  src_t      *src,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  wav_t      *wav,
  mympi_t    *mympi,
  iorecv_t   *iorecv,
  ioline_t   *ioline,
  ioslice_t  *ioslice,
  iosnap_t   *iosnap,
  // time
  float dt, int nt_total, float t0,
  char *output_fname_part,
  char *output_dir,
  int qc_check_nan_num_of_step,
  const int output_all, // qc all var
  const int verbose);

void
sv_eq1st_cart_stg_ac_iso_hook(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose);

void
sv_eq1st_cart_stg_ac_iso_moment(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose);

void
sv_eq1st_cart_stg_ac_iso_hook_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  P,
    float *restrict kappa3d, 
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdrypml_t *bdrypml, bdryfree_t *bdryfree,
    const int myid, const int verbose);

void
sv_eq1st_cart_stg_ac_iso_moment_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  P,
    float *restrict slw3d,
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdrypml_t *bdrypml,
    const int myid, const int verbose);

int
sv_eq1st_cart_stg_ac_iso_hook_src(
    float *restrict P,
    float dt, float dx, float dy, float dz,
    src_t *src, // short nation for reference member
    const int myid, const int verbose);

int
sv_eq1st_cart_stg_ac_iso_moment_src(
    float *restrict Vx, float *restrict Vy, float *restrict Vz,
    float *restrict slw3d, float dt, float dx, float dy, float dz,
    src_t *src, // short nation for reference member
    const int myid, const int verbose);

void
sv_eq1st_cart_stg_ac_iso_free_simg(
    float *restrict  P, 
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2, int nz,
    size_t siz_line, size_t siz_slice,
    const int myid, const int verbose);

void
sv_eq1st_cart_stg_ac_iso_free_vimg(
    float *restrict  Vx, float *restrict  Vy, float *restrict  Vz,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2, int nz,
    size_t siz_line, size_t siz_slice,
    const int myid, const int verbose);

#endif
