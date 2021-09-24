#ifndef GD_INFO_H
#define GD_INFO_H

#include "mympi_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct {
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

  //int nj1_to_glob_ghost0;
  //int ni1_to_glob_ghost0;
  //int nk1_to_glob_ghost0;

  //int nx1_to_glob_phys0;
  //int nx1_to_glob_ghost0;
  //int ny1_to_glob_phys0;
  //int ny1_to_glob_ghost0;
  //int nz1_to_glob_phys0;
  //int nz1_to_glob_ghost0;

  //int glob_phys_ix1; // gloabl start index along x this thread
  //int glob_phys_ix2; // gloabl end index along x
  //int glob_phys_iy1;
  //int glob_phys_iy2;
  //int glob_phys_iz1;
  //int glob_phys_iz2;

  // size of a single var
  //  the following two naming are same
  size_t siz_iy;
  size_t siz_iz;
  size_t siz_icmp;

  size_t siz_line;
  size_t siz_slice;
  size_t siz_volume; // number of points per var

  // curvilinear coord name,
  char **index_name;
  
  //size_t siz_vars; // volume * num_of_vars, not easy for understand, may named with w3d and aux
} gdinfo_t;

/*************************************************
 * function prototype
 *************************************************/

int
gd_info_set(gdinfo_t *const gdinfo,
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
gd_info_lindx_is_inner(int i, int j, int k, gdinfo_t *gdinfo);

int
gd_info_gindx_is_inner(int gi, int gj, int gk, gdinfo_t *gdinfo);

int
gd_info_gindx_is_inner_i(int gi, gdinfo_t *gdinfo);

int
gd_info_gindx_is_inner_j(int gj, gdinfo_t *gdinfo);

int
gd_info_gindx_is_inner_k(int gk, gdinfo_t *gdinfo);

int
gd_info_ind_glphy2lcext_i(int gi, gdinfo_t *gdinfo);

int
gd_info_ind_glphy2lcext_j(int gj, gdinfo_t *gdinfo);

int
gd_info_ind_glphy2lcext_k(int gk, gdinfo_t *gdinfo);

int
gd_info_ind_lcext2glphy_i(int i, gdinfo_t *gdinfo);

int
gd_info_ind_lcext2glphy_j(int j, gdinfo_t *gdinfo);

int
gd_info_ind_lcext2glphy_k(int k, gdinfo_t *gdinfo);

int
gd_info_print(gdinfo_t *gdinfo);

#endif
