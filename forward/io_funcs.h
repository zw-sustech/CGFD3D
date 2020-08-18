#ifndef IO_FUNCS_H
#define IO_FUNCS_H

void
io_build_fname(char *out_dir,
               char *prefix,
               char *subfix,
               int *myid2,
               char *ou_fname);

void
io_build_fname_time(char *out_dir,
                    char *prefix,
                    char *subfix,
                    int *myid2,
                    int  it,
                    char *ou_fname);

void
io_snapshot_export(char *fname,
                   float *restrict var,
                   int nx,
                   int ny,
                   int nz,
                   int *snap_indx,
                   int verbose);

void
io_var3d_export_nc(char   *ou_file,
                   float  *restrict v3d,
                   size_t *restrict v3d_pos,
                   char  **restrict v3d_name,
                   int   number_of_vars,
                   char  **restrict coord_name,
                   int  nx,
                   int  ny,
                   int  nz);

#endif
