/*********************************************************************
 * wavefield for 3d elastic 1st-order equations
 **********************************************************************/

#include <stdio.h>
#include <string.h>

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
                                                   "wf_el3d_1st");

  // position of each var
  size_t *w3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet(num_wave_vars,
                                                         0,
                                                         "wf_el3d_1st");

  // name of each var
  char **w3d_name = (char **) fdlib_mem_malloc_2l_char(num_wave_vars,
                                                       FD_MAX_STRLEN,
                                                       "wf_el3d_1st");

  // set values
  int ivar = WF_EL3D_1ST_SEQ_VX;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Vx");

  ivar = WF_EL3D_1ST_SEQ_VY;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Vy");

  ivar = WF_EL3D_1ST_SEQ_VZ;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Vz");

  ivar = WF_EL3D_1ST_SEQ_TXX;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Txx");

  ivar = WF_EL3D_1ST_SEQ_TYY;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Tyy");

  ivar = WF_EL3D_1ST_SEQ_TZZ;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Tzz");

  ivar = WF_EL3D_1ST_SEQ_TXZ;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Txz");

  ivar = WF_EL3D_1ST_SEQ_TYZ;
  w3d_pos[ivar] = ivar * siz_volume;
  strcpy(w3d_name[ivar],"Tyz");

  ivar = WF_EL3D_1ST_SEQ_TXY;
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
  for (int ivar=0; i<WF_EL_1ST_NVAR; i++)
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
