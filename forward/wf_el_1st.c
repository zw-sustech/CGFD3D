/*********************************************************************
 * wavefield for 3d elastic 1st-order equations
 **********************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "netcdf.h"

#include "fdlib_mem.h"
#include "fd_t.h"
#include "wf_el_1st.h"

void 
wf_el_1st_init_vars(
    size_t siz_volume,
    int number_of_levels,
    int *number_of_vars,
    float  **p_w3d,
    size_t **p_w3d_pos,
    char  ***p_w3d_name)
{
  const int num_wave_vars = WF_EL_1ST_NVAR; // 9
  /*
   * 0-3: Vx,Vy,Vz
   * 4-9: Txx,Tyy,Tzz,Txz,Tyz,Txy
   */

  // vars
  // 3 Vi, 6 Tij, 4 rk stages
  float *w3d = (float *) fdlib_mem_calloc_1d_float(siz_volume * num_wave_vars * number_of_levels,
                                                   0.0,
                                                   "w3d, wf_el3d_1st");

  // position of each var
  size_t *w3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet(num_wave_vars,
                                                         0,
                                                         "w3d_pos, wf_el3d_1st");

  // name of each var
  char **w3d_name = (char **) fdlib_mem_malloc_2l_char(num_wave_vars,
                                                       FD_MAX_STRLEN,
                                                       "w3d_name, wf_el3d_1st");

  // set values
  int ivar = WF_EL_1ST_SEQ_Vx;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Vx");

  ivar = WF_EL_1ST_SEQ_Vy;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Vy");

  ivar = WF_EL_1ST_SEQ_Vz;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Vz");

  ivar = WF_EL_1ST_SEQ_Txx;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Txx");

  ivar = WF_EL_1ST_SEQ_Tyy;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Tyy");

  ivar = WF_EL_1ST_SEQ_Tzz;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Tzz");

  ivar = WF_EL_1ST_SEQ_Txz;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Txz");

  ivar = WF_EL_1ST_SEQ_Tyz;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Tyz");

  ivar = WF_EL_1ST_SEQ_Txy;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Txy");

  // set return values
  *number_of_vars = num_wave_vars;
  *p_w3d = w3d;
  *p_w3d_pos = w3d_pos;
  *p_w3d_name = w3d_name;
}

void
wf_el_1st_check_value(float *restrict w, size_t siz_volume)
{
  for (int ivar=0; ivar<WF_EL_1ST_NVAR; ivar++)
  {
    float *ptr = w + ivar * siz_volume;
    for (size_t iptr=0; iptr<siz_volume; iptr++)
    {
      if (ptr[iptr] != ptr[iptr])
      {
        fprintf(stderr, "ERROR: NaN occurs at iptr=%d ivar=%d\n", iptr, ivar);
        fflush(stderr);
        exit(-1);
      }
    }
  }
}
