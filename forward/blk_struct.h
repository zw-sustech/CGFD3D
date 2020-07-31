#ifndef BLK_STRUCT_H
#define BLK_STRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct blk_struct{

    // size
    size_t siz_line;
    size_t siz_slice;
    size_t siz_volume; // number of points per var
    size_t siz_vars; // volume * num_of_vars


    //
    // grid index
    //
    int ni1, ni2, nk1, nk2, ni, nk;
    int nx1, nx2, nz1, nz2, nx, nz;

    //
    // grid metrics
    //  x3d, y3d, z3d, jac, xi_x, etc
    float  *g3d; // grid 3d vars
    int     number_of_grid_vars;
    size_t *g3d_pos;
    char  **g3d_name;

    //
    // media
    //  rho, lambda, mu etc
    float  *m3d; // media 3d vars
    int     number_of_medium_vars;
    size_t *m3d_pos;
    char  **m3d_name;

    // Wavefield
    //
    float  *w3d; // wavefield 3d vars
    //size_t pos_vx, pos_vz, pos_txx, pos_tyy
    int     number_of_wave_vars;
    int     number_of_wave_levels;
    size_t *w3d_pos;
    char  **w3d_name;

    // auxiliary vars
    //
    float  *aux; // auxiliary vars
    int     number_of_aux_vars;
    int     number_of_aux_level;
    int     aux_cur_var_indx;  // the index of current var
    size_t *aux_cur_p; // current position of aux
    size_t *aux_pos;
    char  **aux_name;

    // boundary 
    int *boundary_itype;

    // abs
    int     abs_numbers[2*3];

    size_t  abs_coefs_size;
    size_t  abs_coefs_dimpos[2*3];
    float   *abs_coefs;

    size_t  abs_vars_size;
    float   *abs_vars;
    size_t  abs_vars_dimpos[2*3];

    // free surface
    float *matVx2Vz, *matVy2Vz;

    // mem usage
    size_t number_of_float;
    size_t number_of_btye;
};

#endif
