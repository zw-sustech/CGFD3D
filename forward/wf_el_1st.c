
/*********************************************************************
 * wavefield for 3d elastic 1st-order equations
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "netcdf.h"
#include "wf_el3d_1st.h"

// not finished
int wf_el_1st_init_vars(size_t siz_volume, int number_of_levels, int *number_of_vars,
      float **p_w3d, size_t **p_w3d_pos, char ***p_w3d_name)
{
    int ivar;

    *number_of_vars = 9;
    /*
     * 0-3: Vx,Vy,Vz
     * 4-9: Txx,Tyy,Tzz,Txz,Tyz,Txy
     */

    // vars
    // 3 Vi, 6 Tij, 4 rk stages
    *p_w3d = (float *) fdlib_mem_calloc_1d_float( 
                 siz_volume * (*number_of_vars) * number_of_levels, 0.0, "wf_el3d_1st");

    // position of each var
    *p_w3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet( 
                 *number_of_vars, 0, "wf_el3d_1st");

    // name of each var
    *p_w3d_name = (char **) fdlib_mem_malloc_2d_char( 
                 *number_of_vars, FDSYS_MAX_STR_LEN, "wf_el3d_1st");

    // set values
    ivar = WF_EL3D_1ST_SEQ_VX;
    (*p_w3d_pos)[ivar] = ivar * siz_volume;
    strcpy((*p_w3d_name)[ivar],"Vx");

    ivar = WF_EL3D_1ST_SEQ_VY;
    (*p_w3d_pos)[ivar] = ivar * siz_volume;
    strcpy((*p_w3d_name)[ivar],"Vy");

    ivar = WF_EL3D_1ST_SEQ_VZ;
    (*p_w3d_pos)[ivar] = ivar * siz_volume;
    strcpy((*p_w3d_name)[ivar],"Vz");

    // need to set other vars

    return 0;
}

// set cfs pml vars
// not finished
int wf_el_1st_init_cfspml_vars(size_t ni, size_t nj, size_t nk, int number_of_levels, int number_of_vars,
      int *restrict boundary_itype,
      int *restrict abs_number_of_layers, size_t *restrict abs_vars_dimpos, float **p_abs_vars)
{
    int ivar;
    size_t siz_aux_volume;

    *number_of_vars = 9;
    /*
     * 0-3: Vx,Vy,Vz
     * 4-9: Txx,Tyy,Tzz,Txz,Tyz,Txy
     */

    size_pml_vars = 0;

    // x1 x2
    for (i=0; i<2)
    {
      abs_vars_dimpos[i] = size_pml_vars;

      if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML) {
        // contain all vars at each side, include rk scheme 4 levels vars
        size_pml_vars += abs_number_of_layers[i] * nj * nk * number_of_vars * number_of_levels;
      }
    }

    // y1 y2
    for (i=3; i<4)
    {
      abs_vars_dimpos[i] = size_pml_vars;

      if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML) {
        size_pml_vars += abs_number_of_layers[i] * ni * nk * number_of_vars * number_of_levels;
      }
    }

    // z1 z2
    for (i=3; i<4)
    {
      abs_vars_dimpos[i] = size_pml_vars;

      if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML) {
        size_pml_vars += abs_number_of_layers[i] * ni * nj * number_of_vars * number_of_levels;
      }
    }

    // vars
    *p_abs_vars = (float *) fdlib_mem_calloc_1d_float( 
                 siz_pml_vars, 0.0, "wf_el3d_1st_cfspml");

    // set
    *abs_vars_size = size_pml_vars;

    // need to set other vars

    return 0;
}
