#ifndef PAR_T_H
#define PAR_T_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PAR_MAX_STRLEN 1000

struct par_t{
  // dir and file name
  char * in_grid_file;
  char * in_media_file;
  
  // grid
  int  coord_by_import;
  int  coord_by_cartesian;
  int  coord_by_vmap;
  int  coord_by_refine;
  char in_coord_fname[PAR_MAX_STRLEN];

  float cartesian_grid_x0;
  float cartesian_grid_y0;
  float cartesian_grid_z0;
  float cartesian_grid_dx;
  float cartesian_grid_dy;
  float cartesian_grid_dz;
  int   cartesian_grid_nx;
  int   cartesian_grid_ny;
  int   cartesian_grid_nz;

  float grid_x0;
  float grid_dx;
  int   grid_nx;
  float grid_y0;
  float grid_dy;
  int   grid_ny;
  float grid_z0;
  float grid_dz;
  int   grid_nz;

  int   metric_by_import;

  // medium
  int medium_by_import;
  
  // PML
  float cfspml_alpha_max[FD_NDIM][2];
  float cfspml_beta_max[FD_NDIM][2];
  float cfspml_velocity[FD_NDIM][2];
  
  // size
  size_t number_of_points;
  size_t 
};

