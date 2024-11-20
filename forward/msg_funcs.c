/*
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>


/*
 * message to be printed for log
 */

void
msg_on_start(MPI_Comm  comm,
             int myid,
             char *processor_name,
             const int verbose)
{
  char *msg=
"CGFD3D_MPI, last update: 2024/11/20.\n"
"Major features:\n"
"+ Surface force source\n"
"+ elastic wave simulation\n"
"+ ADE CFS-PML\n"
"+ Traction-imaging free surface implementation\n"
"+ singl/finite-fault/area sources\n"
"+ station/line/snapshot output\n";

  fprintf(stdout,"%s\n",msg);
  fprintf(stdout,"myid=%d at host %s\n", myid, processor_name); 
  fflush(stdout);
  return;
}

