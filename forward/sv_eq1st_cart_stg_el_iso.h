#ifndef SV_EQ1ST_CART_STG_EL_ISO_H
#define SV_EQ1ST_CART_STG_EL_ISO_H

#include "fd_t.h"
#include "gd_info.h"
#include "mympi_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_t.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

void
sv_eq1st_cart_stg_el_iso_allstep(
  fdstg_t    *fd,
  gdinfo_t   *gdinfo,
  gd_t       *gdcart,
  md_t       *md,
  src_t      *src,
  bdry_t    *bdry,
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
sv_eq1st_cart_stg_el_iso_hook(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdry_t    *bdry,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose);

void
sv_eq1st_cart_stg_el_iso_moment(
  fdstg_t *fd,
  wav_t  *wav,
  gdinfo_t   *gdinfo,
  md_t   *md,
  bdry_t    *bdry,
  src_t *src,
  float dt, float dx, float dy, float dz,
  const int myid, const int verbose);

void
sv_eq1st_cart_stg_el_iso_hook_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict lam3d, float *restrict  mu3d,
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdry_t    *bdry,
    const int myid, const int verbose);

void
sv_eq1st_cart_stg_el_iso_moment_cfspml(
    fdstg_t *fd,
    float *restrict  Vx , float *restrict  Vy , float *restrict  Vz ,
    float *restrict  Txx, float *restrict  Tyy, float *restrict  Tzz,
    float *restrict  Txz, float *restrict  Tyz, float *restrict  Txy,
    float *restrict slw3d,
    float dt, float dx, float dy, float dz,
    int nk2, size_t siz_line, size_t siz_slice,
    bdry_t    *bdry,
    const int myid, const int verbose);

int
sv_eq1st_cart_stg_el_iso_hook_src(
    float *restrict Txx, float *restrict Tyy, float *restrict Tzz,
    float *restrict Txz, float *restrict Tyz, float *restrict Txy,
    float dt, float dx, float dy, float dz,
    src_t *src, // short nation for reference member
    const int myid, const int verbose);

int
sv_eq1st_cart_stg_el_iso_moment_src(
    float *restrict Vx, float *restrict Vy, float *restrict Vz,
    float *restrict slw3d, float dt, float dx, float dy, float dz,
    src_t *src, // short nation for reference member
    const int myid, const int verbose);

int
sv_eq1st_cart_stg_el_iso_dvh2dvz(gdinfo_t   *gdinfo,
                                 md_t       *md,
                                 bdry_t      *bdry,
                                 const int verbose);

void
sv_eq1st_cart_stg_el_iso_free_simg(
    float *restrict  Tzz, float *restrict  Txz, float *restrict  Tyz,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2, int nz,
    size_t siz_line, size_t siz_slice,
    const int myid, const int verbose);

void
sv_eq1st_cart_stg_el_iso_free_vzero(
    float *restrict  Vx, float *restrict  Vy, float *restrict  Vz,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2, int nz,
    size_t siz_line, size_t siz_slice,
    const int myid, const int verbose);

#endif
