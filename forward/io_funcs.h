#ifndef IO_FUNCS_H
#define IO_FUNCS_H

void
io_snapshot_export(char *fname,
                   float *restrict var,
                   int nx,
                   int ny,
                   int nz,
                   int *snap_indx,
                   int verbose);

#endif
