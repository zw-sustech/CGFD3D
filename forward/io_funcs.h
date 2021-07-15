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
io_snapshot_export_binary(char *fname,
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

int
io_read_locate_station(char *in_filenm, 
                       int   glob_phys_ix1, // gloabl start index along x this thread
                       int   glob_phys_ix2, // gloabl end index along x
                       int   glob_phys_iy1,
                       int   glob_phys_iy2,
                       int   glob_phys_iz1,
                       int   glob_phys_iz2,
                       int   ni1,
                       int   nj1,
                       int   nk1,
                       size_t siz_line,
                       size_t siz_slice,
                       float *x3d, float *y3d, float *z3d,
                       struct fd_sta_all_t *sta_info);

#endif
