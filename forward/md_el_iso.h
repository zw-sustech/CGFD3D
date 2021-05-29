#ifndef MD_EL_ISO_H
#define MD_EL_ISO_H

#define MD_EL_ISO_SEQ_RHO    0
#define MD_EL_ISO_SEQ_LAMBDA 1
#define MD_EL_ISO_SEQ_MU     2

void
md_el_iso_init_vars(
    size_t siz_volume,
    int *number_of_vars, 
    float **p_m3d,
    size_t **p_m3d_pos,
    char ***p_m3d_name);

void
md_el_iso_import(float *restrict m3d, size_t *restrict m3d_pos, char **restrict m3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

void
md_el_iso_export(float  *restrict m3d,
                 size_t *restrict m3d_pos,
                 char  **restrict m3d_name,
                 int number_of_vars,
                 int  nx,
                 int  ny,
                 int  nz,
                 char *fname_coords,
                 char *output_dir);

void
md_el_iso_gen_test(
    float *restrict m3d,
    float *restrict x3d,
    float *restrict y3d,
    float *restrict z3d,
    int nx,
    int ny,
    int nz,
    size_t siz_line,
    size_t siz_slice,
    size_t siz_volume);

void
md_el_iso_rho_to_slow(
    float *restrict m3d,
    size_t siz_volume);

#endif
