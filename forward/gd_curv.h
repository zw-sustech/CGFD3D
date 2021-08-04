#ifndef GD_CURV_H
#define GD_CURV_H

#include <mpi.h>

#include "constants.h"
#include "gd_info.h"

/*************************************************
 * structure
 *************************************************/

//  default means coordinate
typedef struct {
  int n1, n2, n3, n4;
  int nx, ny, nz, ncmp;
  float *v4d; // allocated var

  //to avoid ref x3d at different funcs
  float *x3d; // pointer to var
  float *y3d;
  float *z3d;

  // min/max
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} gdcurv_t;

//  for metric
typedef struct {
  int n1, n2, n3, n4;
  int nx, ny, nz, ncmp;
  float *v4d; // allocated var

  float *jac; // pointer to var
  float *xi_x;
  float *xi_y;
  float *xi_z;
  float *eta_x;
  float *eta_y;
  float *eta_z;
  float *zeta_x;
  float *zeta_y;
  float *zeta_z;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} gdcurv_metric_t;


/*************************************************
 * function prototype
 *************************************************/

void 
gd_curv_init(gdinfo_t *gdinfo, gdcurv_t *gdcurv);

void 
gd_curv_metric_init(gdinfo_t        *gdinfo,
                    gdcurv_metric_t *metric);
void
gd_curv_metric_cal(gdinfo_t        *gdinfo,
                   gdcurv_t        *gdcurv,
                   gdcurv_metric_t *metric,
                   int fd_len, int *restrict fd_indx, float *restrict fd_coef);
void
gd_curv_metric_exchange(gdinfo_t        *gdinfo,
                        gdcurv_metric_t *metric,
                        int             *neighid,
                        MPI_Comm        topocomm);

void
gd_curv_gen_cart(
  gdcurv_t *gdcurv,
  float dx, float x0,
  float dy, float y0,
  float dz, float z0);

void
gd_curv_metric_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

void
gd_curv_coord_import(float *restrict g3d, size_t *restrict g3d_pos, char **restrict g3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, char *fname_coords);

void
gd_curv_coord_export(
  gdinfo_t *gdinfo,
  gdcurv_t *gdcurv,
  char *fname_coords,
  char *output_dir);

void
gd_curv_metric_export(gdinfo_t        *gdinfo,
                      gdcurv_metric_t *metric,
                      char *fname_coords,
                      char *output_dir);

int
gd_curv_gen_layer(char *in_grid_layer_file,
                      int *grid_layer_resample_factor,
                      int *grid_layer_start,
                      int n_total_grid_x,
                      int n_total_grid_y,
                      int n_total_grid_z,
                      float *restrict x3d,
                      float *restrict y3d,
                      float *restrict z3d,
                      int nx, int ni, int gni1, int fdx_nghosts, 
                      int ny, int nj, int gnj1, int fdy_nghosts, 
                      int nz, int nk, int gnk1, int fdz_nghosts);

int gd_grid_z_interp(float *z3dpart, float *zlayerpart, int *NCellPerlay,
                     int *VmapSpacingIsequal, int nLayers, int nx, int ny);

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

void
gd_curv_set_minmax(gdcurv_t *gdcurv);

#endif
