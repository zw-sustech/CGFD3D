#ifndef BLK_STRUCT_H
#define BLK_STRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct blk_t{

    // size of a single var
    size_t siz_line;
    size_t siz_slice;
    size_t siz_volume; // number of points per var

    //size_t siz_vars; // volume * num_of_vars, not easy for understand, may named with w3d and aux

    //
    // grid index
    //
    int ni1, ni2, nj1, nj2, nk1, nk2, ni, nj, nk;
    int nx1, nx2, ny1, ny2, nz1, nz2, nx, ny, nz;

    //
    // grid metrics
    //  x3d, y3d, z3d, jac, xi_x, etc
    float  *g3d; // grid 3d vars
    int     g3d_num_of_vars;
    size_t *g3d_pos; // for each var
    char  **g3d_name;

    //
    // media
    //  rho, lambda, mu etc
    float  *m3d; // media 3d vars
    int     m3d_num_of_vars;
    size_t *m3d_pos;
    char  **m3d_name;

    // Wavefield
    //
    float  *w3d; // wavefield 3d vars
    //size_t pos_vx, pos_vz, pos_txx, pos_tyy
    int     w3d_num_of_vars;
    int     w3d_num_of_levels;;
    size_t  w3d_size_per_level;
    size_t *w3d_pos; // for vars in a single level
    char  **w3d_name;

    // auxiliary vars
    //
    float  *a3d; // auxiliary vars
    int     a3d_num_of_vars;
    int     a3d_num_of_levels;
    size_t  a3d_size_per_level;
    size_t *a3d_pos;
    char  **a3d_name;

    // boundary 
    int *boundary_itype[ FD_NDIM_2 ];

    // free surface
    float *matVx2Vz, *matVy2Vz;

    // source term
    int num_of_force;
    int num_of_moment;

    //
    // abs
    //
    int     abs_num_of_layers[ FD_NDIM_2 ];

    //  _blk means all dims, eliminate _dimpos var
    size_t  abs_blk_indx[FD_NDIM_2][FD_NDIM_2];

    // may add abs_numb_of_vars later for mpml

    //size_t  abs_blk_coefs_size[FD_NDIM_2]; // size of all coef for one block, seems no use
    //size_t  abs_coefs_dimpos[2*3];

    // abs_blk_coefs[idim] + abs_num_of_layers[idim] for next coef
    float   **abs_blk_coefs; // [dim][A,B,D etc]

    //size_t  abs_blk_vars_size_per_level[FD_NDIM_2]; // 

    size_t  abs_blk_vars_level_size; // size of vars per level
    size_t  abs_blk_vars_siz_volume[FD_NDIM_2]; // size of single var in abs_blk_vars for each block
    size_t  abs_blk_vars_blkpos[FD_NDIM_2]; // start pos of each blk in first level; other level similar
    float   *abs_blk_vars; //  order: vars_one_block -> all block -> one level -> 4 levels
    //size_t  abs_vars_dimpos[2*3];

    //
    // connection to other blk or mpi thread
    //
    size_t size_of_buff;
    int *inn_bdry_blk_id; // blk id
    size_t *inn_bdry_blk_indx_pair; // this blk indx to bdry_blk indx

    int out_mpi_num_of_neig; //
    int    *out_mpi_neig_ids; //
    size_t *out_mpi_neig_buff_size; //
    int *out_mpi_neigh_blk_ids; // mpi id and blk id
    size_t *out_bdry_blk


    // mem usage
    size_t number_of_float;
    size_t number_of_btye;
};

#endif
