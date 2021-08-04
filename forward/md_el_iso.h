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

  float *lambda; // pointer to var
  float *mu;
  float *rho;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} mdeliso_t;

/*************************************************
 * function prototype
 *************************************************/

int
md_el_iso_init(gdinfo_t *gdinfo, mdeliso_t *mdeliso);

int
md_el_iso_import(float *restrict m3d, size_t *restrict m3d_pos, char **restrict m3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

int
md_el_iso_export(gdinfo_t  *gdinfo,
                 mdeliso_t *mdeliso,
                 char *fname_coords,
                 char *output_dir);

int
md_el_iso_gen_test(mdeliso_t *mdeliso);

int
md_el_iso_rho_to_slow(float *restrict rho, size_t siz_volume);


#endif
