
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
int wf_el_1st_init_cfspml_vars(size_t nx, size_t ny, size_t nz, int number_of_levels, int number_of_vars,
      int *restrict abs_numbers, float **p_abs_vars)
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
