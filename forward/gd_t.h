#ifndef GD_CURV_H
#define GD_CURV_H

#include <mpi.h>

#include "constants.h"
#include "mympi_t.h"

#define GD_TILE_NX 4
#define GD_TILE_NY 4
#define GD_TILE_NZ 4

/*************************************************
 * structure
 *************************************************/

typedef enum {

  GD_TYPE_CART = 1,
  GD_TYPE_VMAP = 2,
  GD_TYPE_CURV = 3

} gd_type_t;

//  grid coordinate for both cart, vmap and curv
//    to reduce duplicated functions
typedef struct {

  gd_type_t type;

  int ni;
  int nj;
  int nk;
  int nx;
  int ny;
  int nz;
  int ni1;
  int ni2;
  int nj1;
  int nj2;
  int nk1;
  int nk2;

  int glob_ni; // total points without ghost
  int glob_nj;
  int glob_nk;
  int glob_nx; // total points include ghost
  int glob_ny;
  int glob_nz;

  int npoint_ghosts;
  int fdx_nghosts;
  int fdy_nghosts;
  int fdz_nghosts;

  // global index
  int gni1, gnj1, gnk1; // global index, do not accout ghost point
  int gni2, gnj2, gnk2; // global index
  // new naming
  int ni1_to_glob_phys0;
  int nj1_to_glob_phys0;
  int nk1_to_glob_phys0;
  int ni2_to_glob_phys0;
  int nj2_to_glob_phys0;
  int nk2_to_glob_phys0;
  // include ghosts
  int nx1_to_glob_halo0;
  int ny1_to_glob_halo0;
  int nz1_to_glob_halo0;

  int n1, n2, n3, n4;
  int ncmp;
  float *v4d; // allocated var

  //to avoid ref x3d at different funcs
  float *x3d; // pointer to var
  float *y3d;
  float *z3d;

  // for cart grid
  float *x1d;
  float *y1d;
  float *z1d;
  float dx;
  float dy;
  float dz;
  // x0/y0/z0 for grid gen
  float x0_glob;
  float y0_glob;
  float z0_glob;

  // min/max of this thread including ghost points
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;

  // min/max of this thread for points in physical region
  float xmin_phy, xmax_phy;
  float ymin_phy, ymax_phy;
  float zmin_phy, zmax_phy;

  // boundary of each cell for AABB algorithm
  float *cell_xmin;
  float *cell_xmax;
  float *cell_ymin;
  float *cell_ymax;
  float *cell_zmin;
  float *cell_zmax;
  // boundary of tiles by 4x4x4 partition for AABB algorithm
  int   tile_istart[GD_TILE_NX];
  int   tile_iend  [GD_TILE_NX];
  int   tile_jstart[GD_TILE_NY];
  int   tile_jend  [GD_TILE_NY];
  int   tile_kstart[GD_TILE_NZ];
  int   tile_kend  [GD_TILE_NZ];
  float tile_xmin[GD_TILE_NZ][GD_TILE_NY][GD_TILE_NX];
  float tile_xmax[GD_TILE_NZ][GD_TILE_NY][GD_TILE_NX];
  float tile_ymin[GD_TILE_NZ][GD_TILE_NY][GD_TILE_NX];
  float tile_ymax[GD_TILE_NZ][GD_TILE_NY][GD_TILE_NX];
  float tile_zmin[GD_TILE_NZ][GD_TILE_NY][GD_TILE_NX];
  float tile_zmax[GD_TILE_NZ][GD_TILE_NY][GD_TILE_NX];

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t siz_line;
  size_t siz_slice;
  size_t siz_volume; // number of points per var

  size_t *cmp_pos;
  char  **cmp_name;

  // curvilinear coord name,
  char **index_name;
} gd_t;

//  default means coordinate
typedef struct {
  int n1, n2, n3, n4;
  int nx, ny, nz, ncmp;
  float *v4d; // allocated var

  //to avoid ref x3d at different funcs
  float *x3d; // pointer to var
  float *y3d;
  float *z3d;

  // min/max including ghost points
  float xmin, xmax;
  float ymin, ymax;
  float zmin, zmax;

  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t *cmp_pos;
  char  **cmp_name;
} gdcurv_t_nouse;

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
gd_curv_init(gd_t *gdcurv);

void 
gd_curv_metric_init(gd_t        *gdinfo,
                    gdcurv_metric_t *metric);
void
gd_curv_metric_cal(
                   gd_t        *gdcurv,
                   gdcurv_metric_t *metric,
                   int fd_len, int *restrict fd_indx, float *restrict fd_coef);
void
gd_curv_metric_exchange(gd_t        *gdinfo,
                        gdcurv_metric_t *metric,
                        int             *neighid,
                        MPI_Comm        topocomm);

void
gd_curv_gen_cart(
  gd_t *gdcurv,
  float dx, float x0,
  float dy, float y0,
  float dz, float z0);

int
gd_curv_metric_import(gd_t        *gdcurv,
                      gdcurv_metric_t *metric,
                      int is_parallel_netcdf,
                      MPI_Comm comm, 
                      char *fname_coords, char *import_dir);

int
gd_curv_coord_import(gd_t *gdcurv,
                     int is_parallel_netcdf,
                     MPI_Comm comm, 
                     char *fname_coords,
                     char *import_dir);

int
gd_curv_coord_export(
  gd_t *gdcurv,
  int is_parallel_netcdf,
  MPI_Comm comm, 
  char *fname_coords,
  char *output_dir);

void
gd_cart_coord_export(
  gd_t *gdcart,
  char *fname_coords,
  char *output_dir);

int
gd_curv_metric_export(gd_t        *gdinfo,
                      gdcurv_metric_t *metric,
                      int is_parallel_netcdf,
                      MPI_Comm comm, 
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
gd_curv_set_minmax(gd_t *gdcurv);

void 
gd_cart_init_set(gd_t *gdcart,
  float dx, float x0_glob,
  float dy, float y0_glob,
  float dz, float z0_glob);

int
gd_cart_coord_to_glob_indx(
                           gd_t *gdcart,
                           float sx,
                           float sy,
                           float sz,
                           MPI_Comm comm,
                           int myid,
                           int   *ou_si, int *ou_sj, int *ou_sk,
                           float *ou_sx_inc, float *ou_sy_inc, float *ou_sz_inc);

int
gd_curv_coord_to_glob_indx(
                           gd_t *gdcurv,
                           float sx,
                           float sy,
                           float sz,
                           MPI_Comm comm,
                           int myid,
                           int   *ou_si, int *ou_sj, int *ou_sk,
                           float *ou_sx_inc, float *ou_sy_inc, float *ou_sz_inc);


int
gd_curv_coord_to_local_indx(
                        gd_t *gd,
                        float sx, float sy, float sz,
                        int *si, int *sj, int *sk,
                        float *sx_inc, float *sy_inc, float *sz_inc);

int
gd_curv_coord2shift_sample(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    int    nx_sample,
    int    ny_sample,
    int    nz_sample,
    float *si_shift, // interped curv coord
    float *sj_shift,
    float *sk_shift);

  int
gd_curv_coord2index_sample(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_j,
    float *points_k,
    int    nx_sample,
    int    ny_sample,
    int    nz_sample,
    float *si_curv, // interped curv coord
    float *sj_curv,
    float *sk_curv);

  int
gd_curv_coord2index_rdinterp(float sx, float sy, float sz, 
    int num_points,
    float *points_x, // x coord of all points
    float *points_y,
    float *points_z,
    float *points_i, // curv coord of all points
    float *points_j,
    float *points_k,
    float *si_curv, // interped curv coord
    float *sj_curv,
    float *sk_curv);

int
gd_curv_depth_to_axis(
                      gd_t *gdcurv,
                      float sx,
                      float sy,
                      float *sz,
                      MPI_Comm comm,
                      int myid);

float
gd_coord_get_x(gd_t *gd, int i, int j, int k);

float
gd_coord_get_y(gd_t *gd, int i, int j, int k);

float
gd_coord_get_z(gd_t *gd, int i, int j, int k);

int isPointInHexahedron_c(float px,  float py,  float pz,
                        float *vx, float *vy, float *vz);

int point2face(float *hexa1d,float *point, float *p2f);

int face_normal(float (*hexa2d)[3], float *normal_unit);

int
gd_print(gd_t *gd, int verbose);

int
gd_indx_set(gd_t *const gd,
            const mympi_t *const mympi,
            const int number_of_total_grid_points_x,
            const int number_of_total_grid_points_y,
            const int number_of_total_grid_points_z,
                  int abs_num_of_layers[][2],
            const int fdx_nghosts,
            int const fdy_nghosts,
            const int fdz_nghosts,
            const int verbose);

int
gd_lindx_is_inner(int i, int j, int k, gd_t *gdinfo);

int
gd_gindx_is_inner(int gi, int gj, int gk, gd_t *gdinfo);

int
gd_gindx_is_inner_i(int gi, gd_t *gdinfo);

int
gd_gindx_is_inner_j(int gj, gd_t *gdinfo);

int
gd_gindx_is_inner_k(int gk, gd_t *gdinfo);

int
gd_indx_glphy2lcext_i(int gi, gd_t *gdinfo);

int
gd_indx_glphy2lcext_j(int gj, gd_t *gdinfo);

int
gd_indx_glphy2lcext_k(int gk, gd_t *gdinfo);

int
gd_indx_lcext2glphy_i(int i, gd_t *gdinfo);

int
gd_indx_lcext2glphy_j(int j, gd_t *gdinfo);

int
gd_indx_lcext2glphy_k(int k, gd_t *gdinfo);

#endif
