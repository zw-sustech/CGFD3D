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

  // for stg
  MPI_Request r_reqs_vel[4];
  MPI_Request s_reqs_vel[4];
  MPI_Request r_reqs_stress[4];
  MPI_Request s_reqs_stress[4];

  // for col scheme
  MPI_Request r_reqs[4];
  MPI_Request s_reqs[4];

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
