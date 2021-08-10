#ifndef MD_EL_ISO_H
#define MD_EL_ISO_H

#include "gd_info.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
  int n1, n2, n3, n4;
  int nx, ny, nz, ncmp;
  float *v4d; // allocated var

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;

  // flag to determine medium type
  int medium_type;

  // rho for all media
  float *rho;

  // for acustic
  float *kappa; // pointer to var

  // for isotropic media
  float *lambda; // pointer to var
  float *mu;

  // for anisotropic media
  float *c11;
  float *c12;
  float *c13;
  float *c14;
  float *c15;
  float *c16;
  float *c22;
  float *c23;
  float *c24;
  float *c25;
  float *c26;
  float *c33;
  float *c34;
  float *c35;
  float *c36;
  float *c44;
  float *c45;
  float *c46;
  float *c55;
  float *c56;
  float *c66;

} md_t;

/*************************************************
 * function prototype
 *************************************************/

int
md_init(gdinfo_t *gdinfo, md_t *md, int media_type);

int
md_import(float *restrict m3d, size_t *restrict m3d_pos, char **restrict m3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

int
md_export(gdinfo_t  *gdinfo,
                 md_t *md,
                 char *fname_coords,
                 char *output_dir);

int
md_gen_test_el_iso(md_t *md);

int
md_gen_test_el_aniso(md_t *md);

int
md_rho_to_slow(float *restrict rho, size_t siz_volume);


#endif
