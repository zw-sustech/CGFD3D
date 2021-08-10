/*
 * prep funcs for free surface condition
 *   implementation of the condition is inside sv_
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "bdry_free.h"

/*
 * matrix for velocity gradient conversion
 *  only implement z2 (top) right now
 */

int
bdry_free_set(gdinfo_t        *gdinfo,
              bdryfree_t      *bdryfree,
              int   *neighid, 
              int   in_is_sides[][2],
              const int verbose)
{
  int ierr = 0;

  size_t siz_slice  = gdinfo->siz_iz;

  // default disable
  bdryfree->is_enable = 0;

  // check each side
  for (int idim=0; idim<CONST_NDIM; idim++)
  {
    for (int iside=0; iside<2; iside++)
    {
      int ind_1d = iside + idim * 2;

      bdryfree->is_at_sides  [idim][iside] = in_is_sides[idim][iside];

      // reset 0 if not mpi boundary
      if (neighid[ind_1d] != MPI_PROC_NULL)
      {
        bdryfree->is_at_sides  [idim][iside] = 0;
      }

      // enable if any side valid
      if (bdryfree->is_at_sides  [idim][iside] == 1) {
        bdryfree->is_enable = 1;
      }
    } // iside
  } // idim

  // following only implement z2 (top) right now
  float *matVx2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * CONST_NDIM * CONST_NDIM,
                                      0.0,
                                      "bdry_free_set");

  float *matVy2Vz = (float *)fdlib_mem_calloc_1d_float(
                                      siz_slice * CONST_NDIM * CONST_NDIM,
                                      0.0,
                                      "bdry_free_set");

  bdryfree->matVx2Vz2 = matVx2Vz;
  bdryfree->matVy2Vz2 = matVy2Vz;

  return ierr;
}

