#ifndef GD_CURV_H
#define GD_CURV_H

#define GD_CURV_SEQ_X3D 0
#define GD_CURV_SEQ_Y3D 1
#define GD_CURV_SEQ_Z3D 2

#define GD_CURV_SEQ_JAC 0
#define GD_CURV_SEQ_XIX 1
#define GD_CURV_SEQ_XIY 2
#define GD_CURV_SEQ_XIZ 3
#define GD_CURV_SEQ_ETX 4
#define GD_CURV_SEQ_ETY 5
#define GD_CURV_SEQ_ETZ 6
#define GD_CURV_SEQ_ZTX 7
#define GD_CURV_SEQ_ZTY 8
#define GD_CURV_SEQ_ZTZ 9

void 
gd_curv_init_c3d(
    size_t siz_volume,
    int *number_of_vars,
    float  **p_c3d,
    size_t **p_c3d_pos,
    char  ***p_c3d_name);

void 
gd_curv_init_g3d(
    size_t siz_volume,
    int *number_of_vars,
    float  **p_g3d,
    size_t **p_g3d_pos,
    char  ***p_g3d_name);

void
gd_curv_metric_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, int *myid2);

void
gd_curv_coord_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, int *myid2);

void
gd_curv_cal_metric(
    float *restrict c3d,
    float *restrict g3d,
    size_t ni1, size_t ni2, size_t nj1, size_t nj2, size_t nk1, size_t nk2,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    size_t fd_len, size_t *restrict fd_indx, float *restrict fd_coef);

void
gd_curv_gen_cart(
    float *restrict c3d,
    size_t siz_voluem,
    size_t nx, float dx, float x0,
    size_t ny, float dy, float y0,
    size_t nz, float dz, float z0);

#endif
