#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fd_t.h"
#include "par_t.h"

#define PAR_MAX_STRLEN 1000

struct par_t{

  //-- dirs and file name
  char project_dir  [PAR_MAX_STRLEN];
  char output_dir   [PAR_MAX_STRLEN];
  char grid_dir     [PAR_MAX_STRLEN];
  char metric_dir   [PAR_MAX_STRLEN];
  char media_dir    [PAR_MAX_STRLEN];
  char source_dir   [PAR_MAX_STRLEN];
  char station_dir  [PAR_MAX_STRLEN];
  char log_file_name[PAR_MAX_STRLEN];

  // MPI
  int number_of_mpiprocs_x;
  int number_of_mpiprocs_y;

  // time step
  int   number_of_time_steps;
  float size_of_time_step ;
  int   time_start_index;
  int   time_end_index;
  float time_start;
  float time_end  ;

  // for each block
  char grid_name[PAR_MAX_STRLEN];
  
  // grid
  int  number_of_total_grid_points_x;
  int  number_of_total_grid_points_y;
  int  number_of_total_grid_points_z;

  int  coord_by_import;
  int  coord_by_cartesian;
  int  coord_by_vmap;
  int  coord_by_refine;
  char coord_in_fname[PAR_MAX_STRLEN];

  float cartesian_grid_x0;
  float cartesian_grid_y0;
  float cartesian_grid_z0;
  float cartesian_grid_dx;
  float cartesian_grid_dy;
  float cartesian_grid_dz;
  //int   cartesian_grid_nx;
  //int   cartesian_grid_ny;
  //int   cartesian_grid_nz;

  int   metric_by_import;

  // medium
  int medium_by_import;
  char medium_in_fname[PAR_MAX_STRLEN];

  // boundary, FD_NDIM_2
  char **boundary_type_name;
  
  // abs
  int abs_num_of_layers[FD_NDIM_2];
  float cfspml_alpha_max[FD_NDIM_2];
  float cfspml_beta_max[FD_NDIM_2];
  float cfspml_velocity[FD_NDIM_2];

  // source
  char source_in_fname[PAR_MAX_STRLEN];

  // output
  int number_of_snapshot;
  char **snapshot_name;
  int *snapshot_index_start;
  int *snapshot_index_count;
  int *snapshot_index_stride;
  int *snapshot_time_start;
  int *snapshot_time_count;
  int *snapshot_time_stride;

  // misc
  int check_nan_every_nummber_of_steps;
};

void
par_mpi_get(char *par_fname, int myid, MPI_Comm comm, struct par_t *par, int verbose);

void
par_read_from_file(char *par_fname, int myid, MPI_Comm comm, struct par_t *par, int verbose);

void 
par_read_from_str(const char *str, struct par_t *par);

void
par_print(struct par_t *par);

#endif
