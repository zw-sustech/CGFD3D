#ifndef MY_MPI_H
#define MY_MPI_H

#include <mpi.h>

#include "constants.h"

/*******************************************************************************
 * structure
 ******************************************************************************/

typedef struct {
  int       nprocx;
  int       nprocy;

  int       myid;
  MPI_Comm  comm;

  int    topoid[2];
  //int    neighid[4];
  int    neighid[CONST_NDIM_2];
  MPI_Comm    topocomm;

  size_t siz_sbuff;
  size_t siz_rbuff;
  float *sbuff;
  float *rbuff;

  // for col scheme
  MPI_Request r_reqs[4];
  MPI_Request s_reqs[4];

  // for stg
  MPI_Request r_reqs_vel[4];
  MPI_Request s_reqs_vel[4];
  MPI_Request r_reqs_stress[4];
  MPI_Request s_reqs_stress[4];

  // for macdrp
  size_t **pair_siz_sbuff_x1;
  size_t **pair_siz_sbuff_x2;
  size_t **pair_siz_sbuff_y1;
  size_t **pair_siz_sbuff_y2;
  size_t **pair_siz_rbuff_x1;
  size_t **pair_siz_rbuff_x2;
  size_t **pair_siz_rbuff_y1;
  size_t **pair_siz_rbuff_y2;
  MPI_Request ***pair_r_reqs;
  MPI_Request ***pair_s_reqs;

} mympi_t;

/*******************************************************************************
 * function prototype
 ******************************************************************************/

int
mympi_set(mympi_t *mympi,
          int number_of_mpiprocs_x,
          int number_of_mpiprocs_y,
          MPI_Comm comm, 
          const int myid, const int verbose);

void
mympi_print(mympi_t *mympi);

#endif
