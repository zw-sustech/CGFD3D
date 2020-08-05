
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

    *number_of_vars = WF_EL_1ST_NVAR; // 9
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
int wf_el_1st_init_cfspml_vars(int number_of_levels, int number_of_vars,
      int *restrict boundary_itype,
      size_t **restrict abs_blk_indx,
      size_t *abs_blk_vars_level_size,
      size_t *restrict abs_blk_vars_siz_volume,
      size_t *restrict abs_blk_vars_blkpos,
      float *p_abs_blk_vars)
{
    size_t i;

    size_t siz_level;

    siz_level = 0;

    for (i=0; i < FD_NDIM_2; i++)
    {
      abs_blk_vars_blkpos[i] = siz_level;

      // init to 0
      abs_blk_vars_siz_volume[i] = 0;

      // set if pml
      if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML)
      {
        abs_blk_vars_siz_volume[i] =   (abs_blk_indx[i][1] - abs_blk_indx[i][0] + 1)
                                     * (abs_blk_indx[i][3] - abs_blk_indx[i][2] + 1)
                                     * (abs_blk_indx[i][5] - abs_blk_indx[i][4] + 1);

        // add to total size
        siz_level += abs_blk_vars_siz_volume[i] * number_of_vars;
      }
    }

    *abs_blk_vars_level_size = siz_level;

    // vars
    // contain all vars at each side, include rk scheme 4 levels vars
    *p_abs_blk_vars = (float *) fdlib_mem_calloc_1d_float( 
                 siz_level * number_of_levels,
                 0.0, "wf_el_1st_init_cfspml_vars");

    return 0;
}

/*
int wf_el_1st_init_cfspml_vars(int number_of_levels, int number_of_vars,
      int *restrict boundary_itype,
      size_t **restrict abs_blk_indx,
      size_t *abs_blk_vars_level_size,
      size_t *restrict abs_blk_vars_siz_volume,
      size_t *restrict abs_blk_vars_blkpos,
      float *p_abs_blk_vars)
{
    size_t i;

    size_t siz_volum;

    ***p_abs_blk_vars = (float **) fdlib_mem_malloc_1d_float( 
               FD_NDIM_2, "wf_el_1st_init_cfspml_vars");

    for (i=0; i < FD_NDIM_2; i++)
    {
      // init to 0
      abs_blk_vars_siz_volume[i] = 0;

      // set if pml
      if (boundary_itype[i] == FD_BOUNDARY_TYPE_CFSPML)
      {
        abs_blk_vars_siz_volume[i] =   (abs_blk_indx[i][1] - abs_blk_indx[i][0] + 1)
                                     * (abs_blk_indx[i][3] - abs_blk_indx[i][2] + 1)
                                     * (abs_blk_indx[i][5] - abs_blk_indx[i][4] + 1);

        // vars
        // contain all vars at each side, include rk scheme 4 levels vars
        (***p_abs_blk_vars)[i] = (float *) fdlib_mem_calloc_1d_float( 
                 abs_blk_vars_siz_volume[i] * number_of_vars * number_of_levels,
                 0.0, "wf_el_1st_init_cfspml_vars");
      }
    }

    return 0;
}
*/
