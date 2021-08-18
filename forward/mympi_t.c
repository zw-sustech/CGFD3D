/*********************************************************************
 * mpi for this package
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mympi_t.h"

//
// set grid size
//
int
mympi_set(mympi_t *mympi,
          int number_of_mpiprocs_x,
          int number_of_mpiprocs_y,
          MPI_Comm comm, 
          const int myid, const int verbose)
{
  int ierr = 0;

  mympi->nprocx = number_of_mpiprocs_x;
  mympi->nprocy = number_of_mpiprocs_y;

  mympi->myid = myid;
  mympi->comm = comm;

  // mpi topo, only consider 2d topo
  int pdims[2]   = {number_of_mpiprocs_x, number_of_mpiprocs_y};
  int periods[2] = {0,0};

  // create Cartesian topology
  MPI_Cart_create(comm, 2, pdims, periods, 0, &(mympi->topocomm));

  // get my local x,y coordinates
  MPI_Cart_coords(mympi->topocomm, myid, 2, mympi->topoid);

  // neighour
  MPI_Cart_shift(mympi->topocomm, 0, 1, &(mympi->neighid[0]), &(mympi->neighid[1]));
  MPI_Cart_shift(mympi->topocomm, 1, 1, &(mympi->neighid[2]), &(mympi->neighid[3]));
  // set z dir for bdry condition
  mympi->neighid[4] = MPI_PROC_NULL;
  mympi->neighid[5] = MPI_PROC_NULL;

  return ierr;
}

void
mympi_print(mympi_t *mympi)
{    
  fprintf(stdout, "\n-------------------------------------------------------\n");
  fprintf(stdout, "print mympi info:\n");
  fprintf(stdout, "-------------------------------------------------------\n\n");

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> media info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //if (blk->media_type == MEDIA_TYPE_LAYER)
  //{
  //    strcpy(str, "layer");
  //}
  //else if (blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    strcpy(str, "grid");
  //}
  //fprintf(stdout, " media_type = %s\n", str);
  //if(blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    fprintf(stdout, "\n --> the media filename is:\n");
  //    fprintf(stdout, " velp_file  = %s\n", blk->fnm_velp);
  //    fprintf(stdout, " vels_file  = %s\n", blk->fnm_vels);
  //    fprintf(stdout, " rho_file   = %s\n", blk->fnm_rho);
  //}
  //fprintf(stdout, "\n");
  
  return;
}
