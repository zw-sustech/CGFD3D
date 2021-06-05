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

#define INDEX( i, j, k, ni, nj ) ( ( i ) + ( j ) * ( ni ) + ( k ) * ( ni ) * ( nj ) )

void 
gd_curv_init_c3d(
    size_t siz_volume,
    int *number_of_vars,
    float  **p_c3d,
    size_t **p_c3d_pos,
    char  ***p_c3d_name,
    char  ***p_coord_name);

void 
gd_curv_init_g3d(
    size_t siz_volume,
    int *number_of_vars,
    float  **p_g3d,
    size_t **p_g3d_pos,
    char  ***p_g3d_name);

void
gd_curv_metric_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

void
gd_curv_coord_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

void
gd_curv_cal_metric(
    float *restrict c3d,
    float *restrict g3d,
    int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
    int nx, int ny, int nz,
    size_t siz_line, size_t siz_slice, size_t siz_volume,
    int fd_len, int *restrict fd_indx, float *restrict fd_coef);

void
gd_curv_gen_cart(
    float *restrict c3d,
    size_t siz_voluem,
    int nx, float dx, float x0,
    int ny, float dy, float y0,
    int nz, float dz, float z0);

void
gd_curv_metric_export(float  *restrict g3d,
                      size_t *restrict g3d_pos,
                      char  **restrict g3d_name,
                      int number_of_vars,
                      int  nx,
                      int  ny,
                      int  nz,
                      int  ni1,
                      int  nj1,
                      int  nk1,
                      int  ni,
                      int  nj,
                      int  nk,
                      int  gni1,
                      int  gnj1,
                      int  gnk1,
                      char *fname_coords,
                      char *output_dir);

void
gd_curv_coord_export(float  *restrict c3d,
                     size_t *restrict c3d_pos,
                     char  **restrict c3d_name,
                     int number_of_vars,
                     int  nx,
                     int  ny,
                     int  nz,
                     int  ni1,
                     int  nj1,
                     int  nk1,
                     int  ni,
                     int  nj,
                     int  nk,
                     int  gni1,
                     int  gnj1,
                     int  gnk1,
                     char *fname_coords,
                     char *output_dir);

int gd_grid_z_interp(int xi, int yi, float* z3d, float* zlayer3d, int* NCellPerlay,
                     int* VmapSpacingIsequal, int nLayers, int nx, int ny );
float gd_seval(int ni, float u,
               int n, float x[], float y[],
               float b[], float c[], float d[],
               int *last);
int gd_SPLine( int n, int end1, int end2,
              float slope1, float slope2,
              float x[], float y[],
              float b[], float c[], float d[],
              int *iflag);
void gd_SPL(int n, float *x, float *y, int ni, float *xi, float *yi);
int gd_curv_gen_layer(float *restrict c3d, int nLayers,
                      int* NCellPerlay, int* VmapSpacingIsequal,
                      int nx, int ni, int gni1, int fdx_nghosts, int n_total_grid_x,
                      int ny, int nj, int gnj1, int fdy_nghosts, int n_total_grid_y,
                      int nz, int nk, int gnk1, int fdz_nghosts, int n_total_grid_z);


#endif
