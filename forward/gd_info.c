/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fdlib_mem.h"
#include "constants.h"
#include "gd_info.h"

//
// set grid size
//

int
gd_info_set(gdinfo_t *const gdinfo,
            const mympi_t *const mympi,
            const int number_of_total_grid_points_x,
            const int number_of_total_grid_points_y,
            const int number_of_total_grid_points_z,
            const int *const abs_num_of_layers,
            const int fdx_nghosts,
            int const fdy_nghosts,
            const int fdz_nghosts,
            const int verbose)
{
  int ierr = 0;

  // determine ni
  int nx_et = number_of_total_grid_points_x;

  // double cfspml load
  nx_et += abs_num_of_layers[0] + abs_num_of_layers[1];

  // partition into average plus left at last
  int nx_avg  = nx_et / mympi->nprocx;
  int nx_left = nx_et % mympi->nprocx;

  // should not be less than 2 * fdx_nghosts
  if (nx_avg < 2 * fdx_nghosts) {
    // error
  }

  // should not be less than abs_num_of_layers
  if (nx_avg<abs_num_of_layers[0] || nx_avg<abs_num_of_layers[1]) {
    // error
  }

  // default set to average value
  int ni = nx_avg;
  // subtract nlay for pml node
  if (mympi->neighid[0] == MPI_PROC_NULL) {
    ni -= abs_num_of_layers[0];
  }
  if (mympi->neighid[1] == MPI_PROC_NULL) {
    ni -= abs_num_of_layers[1];
  }
  // first nx_left node add one more point
  if (mympi->topoid[0] < nx_left) {
    ni++;
  }
  // global index
  if (mympi->topoid[0]==0) {
    gdinfo->gni1 = 0;
  } else {
    gdinfo->gni1 = mympi->topoid[0] * nx_avg - abs_num_of_layers[0];
  }
  if (nx_left != 0) {
    gdinfo->gni1 += (mympi->topoid[0] < nx_left)? mympi->topoid[0] : nx_left;
  }

  // determine nj
  int ny_et = number_of_total_grid_points_y;
  // double cfspml load
  ny_et += abs_num_of_layers[2] + abs_num_of_layers[3];
  int ny_avg  = ny_et / mympi->nprocy;
  int ny_left = ny_et % mympi->nprocy;
  if (ny_avg < 2 * fdy_nghosts) {
    // error
  }
  if (ny_avg<abs_num_of_layers[2] || ny_avg<abs_num_of_layers[3]) {
    // error
  }
  int nj = ny_avg;
  if (mympi->neighid[2] == MPI_PROC_NULL) {
    nj -= abs_num_of_layers[2];
  }
  if (mympi->neighid[3] == MPI_PROC_NULL) {
    nj -= abs_num_of_layers[3];
  }
  // not equal divided points given to first ny_left procs
  if (mympi->topoid[1] < ny_left) {
    nj++;
  }
  // global index
  if (mympi->topoid[1]==0) {
    gdinfo->gnj1 = 0;
  } else {
    gdinfo->gnj1 = mympi->topoid[1] * ny_avg - abs_num_of_layers[2];
  }
  if (ny_left != 0) {
    gdinfo->gnj1 += (mympi->topoid[1] < ny_left)? mympi->topoid[1] : ny_left;
  }

  // determine nk
  int nk = number_of_total_grid_points_z;
  gdinfo->gnk1 = 0;
  
  // add ghost points
  int nx = ni + 2 * fdx_nghosts;
  int ny = nj + 2 * fdy_nghosts;
  int nz = nk + 2 * fdz_nghosts;

  gdinfo->ni = ni;
  gdinfo->nj = nj;
  gdinfo->nk = nk;

  gdinfo->nx = nx;
  gdinfo->ny = ny;
  gdinfo->nz = nz;

  gdinfo->ni1 = fdx_nghosts;
  gdinfo->ni2 = gdinfo->ni1 + ni - 1;

  gdinfo->nj1 = fdy_nghosts;
  gdinfo->nj2 = gdinfo->nj1 + nj - 1;

  gdinfo->nk1 = fdz_nghosts;
  gdinfo->nk2 = gdinfo->nk1 + nk - 1;

  // global index end
  gdinfo->gni2 = gdinfo->gni1 + gdinfo->ni - 1;
  gdinfo->gnj2 = gdinfo->gnj1 + gdinfo->nj - 1;
  gdinfo->gnk2 = gdinfo->gnk1 + gdinfo->nk - 1;

  gdinfo->ni1_to_glob_phys0 = gdinfo->gni1;
  gdinfo->ni2_to_glob_phys0 = gdinfo->gni2;
  gdinfo->nj1_to_glob_phys0 = gdinfo->gnj1;
  gdinfo->nj2_to_glob_phys0 = gdinfo->gnj2;
  gdinfo->nk1_to_glob_phys0 = gdinfo->gnk1;
  gdinfo->nk2_to_glob_phys0 = gdinfo->gnk2;
  
  // x dimention varies first
  gdinfo->siz_line   = nx; 
  gdinfo->siz_slice  = nx * ny; 
  gdinfo->siz_volume = nx * ny * nz;

  // new var, will replace above old naming
  gdinfo->siz_iy   = gdinfo->siz_line;
  gdinfo->siz_iz   = gdinfo->siz_slice;
  gdinfo->siz_icmp = gdinfo->siz_volume;

  // set npoint_ghosts according to fdz_nghosts
  gdinfo->npoint_ghosts = fdz_nghosts;

  gdinfo->index_name = fdlib_mem_malloc_2l_char(
                        CONST_NDIM, CONST_MAX_STRLEN, "gdinfo name");

  // grid coord name
  sprintf(gdinfo->index_name[0],"%s","i");
  sprintf(gdinfo->index_name[1],"%s","j");
  sprintf(gdinfo->index_name[2],"%s","k");

  return ierr;
}

/*
 * give a local index ref, check if in this thread
 */

int
gd_info_lindx_is_inner(int i, int j, int k, gdinfo_t *gdinfo)
{
  int is_in = 0;

  if (   i >= gdinfo->ni1 && i <= gdinfo->ni2
      && j >= gdinfo->nj1 && j <= gdinfo->nj2
      && k >= gdinfo->nk1 && k <= gdinfo->nk2)
  {
    is_in = 1;
  }

  return is_in;
}  

/*
 * give a global index ref to phys0, check if in this thread
 */

int
gd_info_gindx_is_inner(int gi, int gj, int gk, gdinfo_t *gdinfo)
{
  int ishere = 0;

  if ( gi >= gdinfo->ni1_to_glob_phys0 && gi <= gdinfo->ni2_to_glob_phys0 &&
       gj >= gdinfo->nj1_to_glob_phys0 && gj <= gdinfo->nj2_to_glob_phys0 &&
       gk >= gdinfo->nk1_to_glob_phys0 && gk <= gdinfo->nk2_to_glob_phys0 )
  {
    ishere = 1;
  }

  return ishere;
}

/*
 * glphyinx, glextind, gp,ge
 * lcphyind, lcextind
 * gl: global
 * lc: local
 * inx: index
 * phy: physical points only, do not count ghost
 * ext: include extended points, with ghots points
 */

int
gd_info_gindx_is_inner_i(int gi, gdinfo_t *gdinfo)
{
  int ishere = 0;

  if ( gi >= gdinfo->ni1_to_glob_phys0 && gi <= gdinfo->ni2_to_glob_phys0)
  {
    ishere = 1;
  }

  return ishere;
}

int
gd_info_gindx_is_inner_j(int gj, gdinfo_t *gdinfo)
{
  int ishere = 0;

  if ( gj >= gdinfo->nj1_to_glob_phys0 && gj <= gdinfo->nj2_to_glob_phys0)
  {
    ishere = 1;
  }

  return ishere;
}

int
gd_info_gindx_is_inner_k(int gk, gdinfo_t *gdinfo)
{
  int ishere = 0;

  if ( gk >= gdinfo->nk1_to_glob_phys0 && gk <= gdinfo->nk2_to_glob_phys0)
  {
    ishere = 1;
  }

  return ishere;
}

/*
 * convert global index to local
 */

int
gd_info_ind_glphy2lcext_i(int gi, gdinfo_t *gdinfo)
{
  return gi - gdinfo->ni1_to_glob_phys0 + gdinfo->npoint_ghosts;
}

int
gd_info_ind_glphy2lcext_j(int gj, gdinfo_t *gdinfo)
{
  return gj - gdinfo->nj1_to_glob_phys0 + gdinfo->npoint_ghosts;
}

int
gd_info_ind_glphy2lcext_k(int gk, gdinfo_t *gdinfo)
{
  return gk - gdinfo->nk1_to_glob_phys0 + gdinfo->npoint_ghosts;
}

/*
 * convert local index to global
 */

int
gd_info_ind_lcext2glphy_i(int i, gdinfo_t *gdinfo)
{
  return i - gdinfo->npoint_ghosts + gdinfo->ni1_to_glob_phys0;
}

int
gd_info_ind_lcext2glphy_j(int j, gdinfo_t *gdinfo)
{
  return j - gdinfo->npoint_ghosts + gdinfo->nj1_to_glob_phys0;
}

int
gd_info_ind_lcext2glphy_k(int k, gdinfo_t *gdinfo)
{
  return k - gdinfo->npoint_ghosts + gdinfo->nk1_to_glob_phys0;
}

/*
 * print for QC
 */

int
gd_info_print(gdinfo_t *gdinfo)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> grid info:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " nx    = %-10d\n", gdinfo->nx);
  fprintf(stdout, " ny    = %-10d\n", gdinfo->ny);
  fprintf(stdout, " nz    = %-10d\n", gdinfo->nz);
  fprintf(stdout, " ni    = %-10d\n", gdinfo->ni);
  fprintf(stdout, " nj    = %-10d\n", gdinfo->nj);
  fprintf(stdout, " nk    = %-10d\n", gdinfo->nk);

  fprintf(stdout, " ni1   = %-10d\n", gdinfo->ni1);
  fprintf(stdout, " ni2   = %-10d\n", gdinfo->ni2);
  fprintf(stdout, " nj1   = %-10d\n", gdinfo->nj1);
  fprintf(stdout, " nj2   = %-10d\n", gdinfo->nj2);
  fprintf(stdout, " nk1   = %-10d\n", gdinfo->nk1);
  fprintf(stdout, " nk2   = %-10d\n", gdinfo->nk2);

  fprintf(stdout, " ni1_to_glob_phys0   = %-10d\n", gdinfo->gni1);
  fprintf(stdout, " ni2_to_glob_phys0   = %-10d\n", gdinfo->gni2);
  fprintf(stdout, " nj1_to_glob_phys0   = %-10d\n", gdinfo->gnj1);
  fprintf(stdout, " nj2_to_glob_phys0   = %-10d\n", gdinfo->gnj2);

  return(0);
}
