#ifndef MSG_FUNCS_H
#define MSG_FUNCS_H

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
             const int verbose);
#endif
