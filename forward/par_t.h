#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fd_t.h"
#include "par_t.h"

#define PAR_MAX_STRLEN 1000
#define PAR_TYPE_STRLEN 50

struct par_t{

  //-- dirs and file name
  //char project_dir  [PAR_MAX_STRLEN];
  char output_dir   [PAR_MAX_STRLEN];
  char grid_dir     [PAR_MAX_STRLEN];
  //char metric_dir   [PAR_MAX_STRLEN];
  char media_dir    [PAR_MAX_STRLEN];
  //char source_dir   [PAR_MAX_STRLEN];
  //char station_dir  [PAR_MAX_STRLEN];
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
  //char grid_name[PAR_MAX_STRLEN];
  
  // grid size
  int  number_of_total_grid_points_x;
  int  number_of_total_grid_points_y;
  int  number_of_total_grid_points_z;

  // boundary, FD_NDIM_2
  char **boundary_type_name;
  
  // abs
  int   abs_num_of_layers[FD_NDIM_2];
  float cfspml_alpha_max[FD_NDIM_2];
  float cfspml_beta_max[FD_NDIM_2];
  float cfspml_velocity[FD_NDIM_2];

  // grid
  char input_grid_type[PAR_TYPE_STRLEN];
  char input_grid_file[PAR_MAX_STRLEN];

  int  coord_by_import;
  int  coord_by_cartesian;
  int  coord_by_vmap;
  int  coord_by_refine;

  float cartesian_grid_origin[FD_NDIM];
  float cartesian_grid_stepsize[FD_NDIM];

  // metric
  char input_metric_type[PAR_TYPE_STRLEN];
  int  metric_by_import;

  // medium
  char input_medium_type[PAR_TYPE_STRLEN];
  int medium_by_import;
  //char input_medium_file[PAR_MAX_STRLEN];

  // source
  char input_source_file[PAR_MAX_STRLEN];

  // output
  // receiver
  char input_receiver_file[PAR_MAX_STRLEN];
  // line
  int number_of_receiver_line;
  int *receiver_line_index_start;
  int *receiver_line_index_incre;
  int *receiver_line_count;
  //int *receiver_line_time_interval;
  char **receiver_line_name;
  // slice
  int number_of_slice_x;
  int number_of_slice_y;
  int number_of_slice_z;
  int *slice_x_index;
  int *slice_y_index;
  int *slice_z_index;
  // snapshot
  int number_of_snapshot;
  char **snapshot_name;
  int *snapshot_index_start;
  int *snapshot_index_count;
  int *snapshot_index_incre;
  int *snapshot_time_start;
  //int *snapshot_time_count; // should output to end 
  int *snapshot_time_incre;
  int *snapshot_save_velocity;
  int *snapshot_save_stress;
  int *snapshot_save_strain;

  // misc
  int check_nan_every_nummber_of_steps;
  int output_all;
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
