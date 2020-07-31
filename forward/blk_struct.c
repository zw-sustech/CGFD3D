
/*********************************************************************
 This is the main program  for multi-block forward modeling: 

 general block structure for curv, vmap, cart etc structured-grid fd scheme

 blk_struct.c

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "blk_struct.h"

//
// set grid size
//
int blk_struct_set_grid_size(struct blk_struct *blk_w, size_t nx, size_t ny, size_t nz)
{
  int ierr=0;
  
  blk_w->nx = nx;
  blk_w->ny = ny;
  blk_w->nz = nz;
  
  // x dimention varies first
  blk_w->siz_line   = nx; 
  blk_w->siz_slice  = nx * ny; 
  blk_w->siz_volume = nx * ny * nz;
  
  return ierr;
}

//
// set levels for time int scheme
//
int blk_struct_set_wave_level(struct blk_struct *blk_w, int number_of_levels)
{
  int ierr=0;
  
  blk_w->number_of_wave_levels = number_of_levels;
  
  return ierr;
}
