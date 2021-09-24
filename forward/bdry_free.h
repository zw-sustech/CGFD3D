#ifndef BDRY_FREE_H
#define BDRY_FREE_H

#include "constants.h"
#include "gd_info.h"
#include "gd_t.h"

/*************************************************
 * structure
 *************************************************/

typedef struct
{
  int is_enable; //
  int is_at_sides[CONST_NDIM][2];

  int nx;
  int ny;
  int nz;

  // top
  float *matVx2Vz2; // [j,i, dzVi, dxVi]
  float *matVy2Vz2;

  // bottom
  float *matVx2Vz1;
  float *matVy2Vz1;

  // left
  float *matVy2Vx1;
  float *matVz2Vx1;

  // right
  float *matVy2Vx2;
  float *matVz2Vx2;

  // front
  float *matVx2Vy1;
  float *matVz2Vy1;

  // back
  float *matVx2Vy2;
  float *matVz2Vy2;

} bdryfree_t;

/*************************************************
 * function prototype
 *************************************************/

int
bdry_free_set(gdinfo_t        *gdinfo,
              bdryfree_t      *bdryfree,
              int   *neighid, 
              int   in_is_sides[][2],
              const int verbose);

#endif
