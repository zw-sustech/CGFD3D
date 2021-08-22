#ifndef SV_EQ1ST_CURV_COL_H
#define SV_EQ1ST_CURV_COL_H

#include "fd_t.h"
#include "gd_info.h"
#include "mympi_t.h"
#include "gd_t.h"
#include "md_t.h"
#include "wav_t.h"
#include "src_t.h"
#include "bdry_free.h"
#include "bdry_pml.h"
#include "io_funcs.h"

/*************************************************
 * function prototype
 *************************************************/

void
sv_eq1st_curv_col_allstep(
  fd_t            *fd,
  gdinfo_t        *gdinfo,
  gdcurv_metric_t *metric,
  md_t      *md,
  src_t      *src,
  bdryfree_t *bdryfree,
  bdrypml_t  *bdrypml,
  wav_t  *wav,
  mympi_t    *mympi,
  iorecv_t   *iorecv,
  ioline_t   *ioline,
  ioslice_t  *ioslice,
  iosnap_t   *iosnap,
  // time
  float dt, int nt_total, float t0,
  char *output_fname_part,
  char *output_dir,
  int qc_check_nan_num_of_step,
  const int output_all, // qc all var
  const int verbose);

int
sv_eq1st_curv_graves_Qs(float *w, int ncmp, float dt, gdinfo_t *gdinfo, md_t *md);

#endif
