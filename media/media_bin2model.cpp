/******************************************************************************
 *
 * This function is used to discretize a given grid model (binary version) 
 *  to a calculation grid.
 *
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    03/2022: Created by Luqian Jiang 
 *
 *******************************************************************************/
#include <iostream>
#include <vector>
#include <string.h>
#include <math.h>
#include "media_geometry3d.hpp"
#include "media_bin2model.hpp"
#include "media_utility.hpp"
#include "media_read_file.hpp"


/* =================== for C call======================== */
int media_bin2model_el_iso(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    float xmin, float xmax,
    float ymin, float ymax,
    int grid_type,
    int  *bin_order,    // eg, [2, 0, 1]=[z, x, y] 0:x, 1:y, 2:z
    int  *bin_size,     // [ndim1, ndim2, ndim3],
    float  *bin_spacing,  // [dh1, dh2, dh3],
    float  *bin_origin,   // [h0_1, h0_2, h0_3],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs  )
{
  
    // get the dim order and set error reporting
    int dimx = -1, dimy = -1, dimz = -1;

    for (int i = 0; i < 3; i++) {
      if (bin_order[i] == 0) dimx = i;
      else if (bin_order[i] == 1) dimy = i;
      else if (bin_order[i] == 2) dimz = i;
      else {
        fprintf(stderr, "Error: wrong number in bin_order: bin_order[%d] = %d! \n"\
                        "       %d is a wrong order, it must be 0, 1, or 2. \n",
                i, bin_order[i], bin_order[i] );
        fflush(stderr);
        exit(1);
      }
    }

    if (dimx == -1) {
      fprintf(stderr, "Error: missing the order of x in bin_order.\n");
      fprintf(stderr, "       bin_order=%d,%d,%d\n", bin_order[0],bin_order[1],bin_order[2]);
      fflush(stderr);
      exit(1);
    }
    if (dimy == -1) {
      fprintf(stderr, "Error: missing the order of y in bin_order.\n");
      fprintf(stderr, "       bin_order=%d,%d,%d\n", bin_order[0],bin_order[1],bin_order[2]);
      fflush(stderr);
      exit(1);
    }
    if (dimz == -1) {
      fprintf(stderr, "Error: missing the order of z in bin_order.\n");
      fprintf(stderr, "       bin_order=%d,%d,%d\n", bin_order[0],bin_order[1],bin_order[2]);
      fflush(stderr);
      exit(1);
    }
    
    // get the read range index
    // x and y should be increasing
    float  bin_x0 = bin_origin[dimx]; 
    float  bin_y0 = bin_origin[dimy]; 
    float  bin_z0 = bin_origin[dimz]; 
    float  bin_dx = bin_spacing[dimx]; 
    float  bin_dy = bin_spacing[dimy];
    float  bin_dz = bin_spacing[dimz];

    int ix_start = (xmin < bin_x0 ? 0:( xmin - bin_x0 ) / bin_dx); 
    if (ix_start >= bin_size[dimx] || xmax < bin_x0) {
      fprintf(stderr, "Error: The given model does not coincide with the calculation area\n"\
                      "       in x-direction. Please check the model settings!");
      fflush(stderr);
      exit(1);
    }
    int ix_end = ceil( ( xmax - bin_x0 ) / bin_dx );
    if (ix_end >= bin_size[dimx])
      ix_end = bin_size[dimx]-1;

    int iy_start = (ymin < bin_y0 ? 0: (ymin - bin_y0)/ bin_dy); 
    if (iy_start >= bin_size[dimy] || ymax < bin_y0) {
      fprintf(stderr, "Error: The given model does not coincide with the calculation area\n"\
                      "       in y-direction. Please check the model settings!");
      fflush(stderr);
      exit(1);
    }
    int iy_end = ceil( ( ymax - bin_y0 ) / bin_dy ); 
    if (iy_end >= bin_size[dimy])
      iy_end = bin_size[dimy]-1;

    // judge whether the calculation area does coincide with the model in z direction
    float zmin = FLT_MAX, zmax = -FLT_MAX;
    size_t siz_volume =  nx * ny * nz;
    size_t siz_z3d = (grid_type == GRID_CART ? nz:siz_volume);

    for (size_t i = 0; i < siz_z3d; i++) {
      zmin = std::min(zmin, z3d[i]);
      zmax = std::max(zmax, z3d[i]);
    }
    
    float bin_zn = bin_z0 + bin_dz * (bin_size[dimz]-1);
    if (zmin >= std::max(bin_zn, bin_z0) || zmax <= std::min(bin_z0, bin_zn) ) {
      fprintf(stderr, "Error: The given model does not coincide with the calculation area\n"\
                      "       in z-direction. Please check the model settings!");
      fflush(stderr);
      exit(1);
    }

    // the coordinate vectors of the given model (after cutting by [min max] )
    size_t bin_nx = ix_end-ix_start+1;
    size_t bin_ny = iy_end-iy_start+1;
    size_t bin_nz = bin_size[dimz];
    std::vector<float> xvec(bin_nx);
    std::vector<float> yvec(bin_ny);
    std::vector<float> zvec(bin_nz);

    for (size_t ix = 0; ix < bin_nx; ix++) {
      xvec[ix] = bin_x0 + (ix_start+ix)*bin_dx;
    }

    for (size_t iy = 0; iy < bin_ny; iy++) {
      yvec[iy] = bin_y0 + (iy_start+iy)*bin_dy;
    }

    for (size_t iz = 0; iz < bin_nz; iz++) {
      zvec[iz] = bin_z0 + iz*bin_dz;
    }

    int bin_start[3];
    int bin_end[3];
    
    bin_start[dimx] = ix_start; bin_end[dimx] = ix_end;
    bin_start[dimy] = iy_start; bin_end[dimy] = iy_end;
    bin_start[dimz] = 0; bin_end[dimz] = bin_nz-1;

    //- Read bin file 
    size_t bin_line = bin_nx;
    size_t bin_slice = bin_nx * bin_ny;
    size_t bin_volume = bin_slice * bin_nz;
    float *bin_rho = new float[bin_volume];
    float *bin_vp  = new float[bin_volume];
    float *bin_vs  = new float[bin_volume];

    fprintf(stdout, "- reading model file: %s, \n", bin_file_rho);
    read_bin_file(bin_file_rho, bin_rho, dimx, dimy, dimz, bin_start, bin_end, bin_size, bin_line, bin_slice);

    fprintf(stdout, "                      %s, \n", bin_file_vp);
    read_bin_file(bin_file_vp,  bin_vp,  dimx, dimy, dimz, bin_start, bin_end, bin_size, bin_line, bin_slice);

    fprintf(stdout, "                      %s, \n", bin_file_vs);
    read_bin_file(bin_file_vs,  bin_vs,  dimx, dimy, dimz, bin_start, bin_end, bin_size, bin_line, bin_slice);

    // media parameterization
    parameterization_bin_el_iso_loc(rho3d, lam3d, mu3d, x3d, y3d, z3d, nx, ny, nz, grid_type, 
        xvec, yvec, zvec, bin_rho, bin_vp, bin_vs);

    delete [] bin_rho;
    delete [] bin_vp;
    delete [] bin_vs;
    return 0;
}

void parameterization_bin_el_iso_loc(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    int grid_type,
    std::vector<float> &xvec, 
    std::vector<float> &yvec, 
    std::vector<float> &zvec,
    float *bin_rho,
    float *bin_vp,
    float *bin_vs ) 
{

  size_t siz_line   = nx;
  size_t siz_slice  = nx * ny;
  size_t siz_volume = nx * ny * nz;

  float slow_k = 1.0/(nz-1); // for print progress
  std::cout << "- discrete model from the binary file\n\n";
  for (size_t k = 0; k < nz; ++k) {
    printProgress(slow_k*k);
    for (size_t j = 0; j < ny; ++j) {
      for (size_t i = 0; i < nx; ++i) {
        size_t indx =  i + j * siz_line + k * siz_slice;
        size_t indx_z = indx;          // for vmap and curv: z
        size_t indx_x = i, indx_y = j; // for vmap and cart: x, y 
        if (grid_type == GRID_CART) {
          indx_z = k;                // for cart: z
        } else if (grid_type == GRID_CURV) {
          indx_x = indx;
          indx_y = indx;             // for curv: x, y
        }
        float rho, vp, vs;
        rho = TrilinearInterpolation(xvec, yvec, zvec, bin_rho, x3d[indx_x], y3d[indx_y], z3d[indx_z]);    
        vp  = TrilinearInterpolation(xvec, yvec, zvec, bin_vp , x3d[indx_x], y3d[indx_y], z3d[indx_z]);    
        vs  = TrilinearInterpolation(xvec, yvec, zvec, bin_vs , x3d[indx_x], y3d[indx_y], z3d[indx_z]); 

        mu3d[indx]  = vs*vs*rho;
        lam3d[indx] = vp*vp*rho-2.0*mu3d[indx];
        rho3d[indx] = rho;
      }
    }
  }

}



