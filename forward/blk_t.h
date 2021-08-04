#ifndef BLK_T_H
#define BLK_T_H

#include "constants.h"
#include "fd_t.h"
#include "gd_info.h"
#include "mympi_t.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "src_funcs.h"
#include "bdry_free.h"
#include "bdry_pml.h"
#include "io_funcs.h"

/*******************************************************************************
 * structure
 ******************************************************************************/

typedef struct
{
  // name for output file name
  char name[CONST_MAX_STRLEN];

  // fd
  fd_t *fd;

  // mpi
  mympi_t *mympi;

  // grid index info
  gdinfo_t *gdinfo;
  
  // coordnate: x3d, y3d, z3d
  gdcurv_t *gdcurv;

  // grid metrics: jac, xi_x, etc
  gdcurv_metric_t *gdcurv_metric;

  // media: rho, lambda, mu etc
  mdeliso_t *mdeliso;

  // wavefield:
  wfel1st_t *wfel1st;
  
  // source term
  src_t *src;
  
  // free surface
  bdryfree_t *bdryfree;
  
  // pml
  bdrypml_t *bdrypml;
  // exp
  //bdryexp_t *bdryexp;
  
  // io
  iorecv_t  *iorecv;
  ioline_t  *ioline;
  iosnap_t  *iosnap;
  ioslice_t *ioslice;

  // fname and dir
  char output_fname_part[CONST_MAX_STRLEN];
  // wavefield output
  char output_dir[CONST_MAX_STRLEN];
  // seperate grid output to save grid for repeat simulation
  char grid_export_dir[CONST_MAX_STRLEN];
  // seperate medium output to save medium for repeat simulation
  char media_export_dir[CONST_MAX_STRLEN];

  // exchange between blocks
  // point-to-point values
  int num_of_conn_points;
  int *conn_this_indx;
  int *conn_out_blk;
  int *conn_out_indx;
  // interp point
  int num_of_interp_points;
  int *interp_this_indx;
  int *interp_out_blk;
  int *interp_out_indx;
  float *interp_out_dxyz;

  // mem usage
  size_t number_of_float;
  size_t number_of_btye;
} blk_t;

/*******************************************************************************
 * function prototype
 ******************************************************************************/

int
blk_init(blk_t *blk,
         const int myid, const int verbose);

// set str
int
blk_set_output(blk_t *blk,
               mympi_t *mympi,
               char *output_dir,
               char *grid_export_dir,
               char *media_export_dir,
               const int verbose);

void
blk_mesg_init(mympi_t *mympi,
                int ni,
                int nj,
                int nk,
                int fdx_nghosts,
                int fdy_nghosts,
                int num_of_vars);

void
blk_pack_mesg(float *restrict w_cur,float *restrict sbuff,
                 int num_of_vars, gdinfo_t *gdinfo,
                 int   fdx_nghosts, int   fdy_nghosts);

void
blk_unpack_mesg(float *restrict rbuff,float *restrict w_cur,
                 int num_of_vars, gdinfo_t *gdinfo,
                 int   fdx_nghosts, int   fdy_nghosts);

int
blk_print(blk_t *blk);

#endif
