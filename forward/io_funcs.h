#ifndef IO_FUNCS_H
#define IO_FUNCS_H

#include "constants.h"
#include "gd_t.h"
#include "io_funcs.h"
#include "md_t.h"
#include "wav_t.h"

/*************************************************
 * structure
 *************************************************/

// for stations output

// single station
typedef struct
{
  float x;
  float y;
  float z;
  float di;
  float dj;
  float dk;
  int   i;
  int   j;
  int   k;
  int   indx1d[CONST_2_NDIM];
  float *seismo;
  float *vi;
  float *tij;
  int   seq_id; // seqential id, 0-based
  char  name[CONST_MAX_NMLEN];
} iorecv_one_t;

typedef struct
{
  int  file_type_sac;
  int  total_number; // total nr from input
  int  nr_here; // in this thread
  int  max_nt;

  int  save_velocity;
  int  save_stress;
  int  save_strain;

  // for block output
  int nt_per_out;  // time block size
  int nt_this_out;  //time block size of this output
  int it_to_this;   // cur it relative to it0 of cur output
  int it0_to_start ;   // cur it relative to begin

  // for netcdf
  int ncid;
  int varid_vi[CONST_NDIM]; // vel var
  int varid_tij[CONST_NDIM_2]; // stress var
  int varid_eij[CONST_NDIM_2]; // strain var

  iorecv_one_t *recvone;

} iorecv_t;

// line output
typedef struct
{
  int    line_seq; // line number, for name from input file
  int    total_number; // total nr of this line
  int    nr_here; // in this thread

  int     *recv_seq; // recv seq in this line
  size_t  *recv_iptr;
  float  *recv_x; // for sac output
  float  *recv_y; // for sac output
  float  *recv_z; // for sac output

  float  *vi;
  float  *tij;

  char   name[CONST_MAX_NMLEN];

  // for netcdf
  int ncid;
  int varid_vi[CONST_NDIM]; // vel var
  int varid_tij[CONST_NDIM_2]; // stress var
  int varid_eij[CONST_NDIM_2]; // strain var

} ioline_one_t;

typedef struct
{
  int     num_of_lines_total; 
  int     num_of_lines_here; 
  int     max_nt;

  int  save_velocity;
  int  save_stress;
  int  save_strain;

  // for block output
  int nt_per_out;  // time block size
  int nt_this_out;  //time block size of this output
  int it_to_this;   // cur it relative to it0 of cur output
  int it0_to_start ;   // cur it relative to begin

  ioline_one_t *lineone;

  //int    *line_nr; // number of receivers, for name from input file
  //int    *line_seq; // line number, for name from input file
  ////int    **recv_ir;
  ////int    **recv_jr;
  ////int    **recv_kr;
  //int    **recv_seq; // recv seq in this line
  //int    **recv_iptr;
  //float  **recv_x; // for sac output
  //float  **recv_y; // for sac output
  //float  **recv_z; // for sac output
  //float  **recv_seismo;
  //char   **line_name;
} ioline_t;

// snapshot output
typedef struct
{
  // for esti size of working space var
  size_t siz_max_wrk;

  int num_snap_total; // num of snap form .json
  int num_of_snap; // num of snap in this thread

  int *in_this_proc; // if output in this proc

  int *i1;
  int *j1;
  int *k1;
  int *ni;
  int *nj;
  int *nk;
  int *di;
  int *dj;
  int *dk;
  int *it1;
  int *dit;
  int *out_vel;
  int *out_stress;
  int *out_strain;
  int *out_coord;

  int *i1_to_glob;
  int *j1_to_glob;
  int *k1_to_glob;

  // to switch parallel netcdf
  int *ni_total;
  int *nj_total;
  int *nk_total;

  char **fname_main;
} iosnap_t;

// for nc output

typedef struct
{
  int num_of_snap;
  int num_snap_total; // num of snap form .json
  int is_parallel;

  int *in_this_proc; // if output in this proc

  int *ncid;
  int *timeid;
  int *varid_V;  // [num_of_snap*CONST_NDIM];
  int *varid_T;  // [num_of_snap*CONST_NDIM_2];
  int *varid_E;  // [num_of_snap*CONST_NDIM_2];
  int *cur_it ;  // [num_of_snap];

  int *start_i;
  int *start_j;
  int *start_k;

  char **fname;
}
iosnap_nc_t;

/*************************************************
 * function prototype
 *************************************************/

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

int
io_get_nextline(FILE *fp, char *str, int length);

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

void
io_snapshot_export_binary(char *fname,
                   float *restrict var,
                   int nx,
                   int ny,
                   int nz,
                   int *snap_indx,
                   int verbose);


void
io_snapshot_locate(gd_t *gdinfo,
                   iosnap_t *iosnap,
                    int  number_of_snapshot,
                    char **snapshot_name,
                    int *snapshot_index_start,
                    int *snapshot_index_count,
                    int *snapshot_index_incre,
                    int *snapshot_time_start,
                    int *snapshot_time_incre,
                    int *snapshot_save_velocity,
                    int *snapshot_save_stress,
                    int *snapshot_save_strain,
                    int *snapshot_save_coord,
                    char *output_dir);

int
io_snap_nc_create(iosnap_t *iosnap, 
                  iosnap_nc_t *iosnap_nc,
                  gd_t        *gdinfo,
                  char *output_fname_part,
                  float *buff,
                  int is_parallel_netcdf,
                  MPI_Comm comm, 
                  int *topoid);

int
io_snap_nc_put(iosnap_t *iosnap,
               iosnap_nc_t *iosnap_nc,
               gd_t    *gdinfo,
               md_t    *md,
               wav_t   *wav,
               float *restrict w4d,
               float *restrict buff,
               int   nt_total,
               int   it,
               float time,
               int is_run_out_vel,     // for stg, out vel and stress at sep call
               int is_run_out_stress,  // 
               int is_incr_cur_it);     // for stg, should output cur_it once

int
io_snap_pack_buff(float *restrict var,
                  size_t siz_line,
                  size_t siz_slice,
                  int starti,
                  int counti,
                  int increi,
                  int startj,
                  int countj,
                  int increj,
                  int startk,
                  int countk,
                  int increk,
                  float *restrict buff);

int
io_snap_stress_to_strain_eliso(float *restrict lam3d,
                               float *restrict mu3d,
                               float *restrict Txx,
                               float *restrict Tyy,
                               float *restrict Tzz,
                               float *restrict Tyz,
                               float *restrict Txz,
                               float *restrict Txy,
                               float *restrict Exx,
                               float *restrict Eyy,
                               float *restrict Ezz,
                               float *restrict Eyz,
                               float *restrict Exz,
                               float *restrict Exy,
                               size_t siz_line,
                               size_t siz_slice,
                               int starti,
                               int counti,
                               int increi,
                               int startj,
                               int countj,
                               int increj,
                               int startk,
                               int countk,
                               int increk);

int
io_snap_nc_close(iosnap_nc_t *iosnap_nc);

int
iosnap_print(iosnap_t *iosnap);


int
io_recv_read_locate(gd_t     *gd,
                    iorecv_t *iorecv,
                    int       nt_total,
                    int       nt_per_out,
                    int       file_type_sac,
                    int       save_velocity,
                    int       save_stress,
                    int       save_strain,
                    char     *in_filenm,
                    MPI_Comm  comm,
                    int       myid,
                    int       verbose);

int
io_recv_nc_create(iorecv_t *iorecv,
                  float stept,
                  int is_parallel_netcdf,
                  MPI_Comm comm, 
                  int myid,
                  char *fname_mpi,
                  char *output_dir);

int
io_recv_keep(iorecv_t *iorecv, float *restrict w4d,
             int it, int siz_icmp);

int
io_recv_nc_put(iorecv_t *iorecv,
               md_t     *md,
               int it,
               int is_parallel_netcdf);

int
io_recv_nc_close(iorecv_t *iorecv, int is_parallel_netcdf);

int
io_recv_output_sac(iorecv_t *iorecv,
                   md_t     *md,
                   float dt,
                   char **cmp_name,
                   char *evtnm,
                   char *output_dir,
                   char *err_message);

int
iorecv_print(iorecv_t *iorecv);

int
io_line_locate(gd_t     *gd,
               ioline_t *ioline,
               int       nt_total,
               int       nt_per_out,
               int       save_velocity,
               int       save_stress,
               int       save_strain,
               int    number_of_receiver_line,
               int   *receiver_line_index_start,
               int   *receiver_line_index_incre,
               int   *receiver_line_count,
               char **receiver_line_name);

int
io_line_nc_create(ioline_t *ioline,
                  float stept,
                  int is_parallel_netcdf,
                  MPI_Comm comm, 
                  int myid,
                  char *fname_mpi,
                  char *output_dir);

int
io_line_keep(ioline_t *ioline, float *restrict w4d,
             int it, int siz_icmp);

int
io_line_nc_put(ioline_t *ioline,
               md_t     *md,
               int it,
               int is_parallel_netcdf);

int
io_line_nc_close(ioline_t *ioline, int is_parallel_netcdf);

int
PG_slice_output(float *PG,  gd_t *gdinfo, 
                float *buff,
                int is_parallel_netcdf,
                MPI_Comm comm, 
      char *output_dir, char *frame_coords, int* topoid);

#endif
