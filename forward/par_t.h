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
  
  // PML
  //
  
  // grid
  int  coord_by_import;
  int  coord_by_cartesian;
  int  coord_by_vmap;
  int  coord_by_refine;
  char in_coord_fname[PAR_MAX_STRLEN];

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
  
  // size
  size_t number_of_points;
  size_t 
};

