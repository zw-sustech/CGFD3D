#ifndef _MEDIA_BIN2MODEL_
#define _MEDIA_BIN2MODEL_

#include "media_geometry3d.hpp"
#include "media_read_file.hpp"

#ifdef __cplusplus
extern "C" {
#endif
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
    int *bin_order,    // eg, [2, 0, 1]=[z, x, y] 0:x, 1:y, 2:z
    int *bin_size,     // [ndim1, ndim2, ndim3],
    float  *bin_spacing,  // [dh1, dh2, dh3],
    float  *bin_origin,   // [h0_1, h0_2, h0_3],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs  );

int media_bin2model_vis_iso(
    float *Qp3d,
    float *Qs3d,
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
    const char *bin_file_Qp,
    const char *bin_file_Qs  );

#ifdef __cplusplus
}
#endif

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
    float *bin_vs ); 

void parameterization_bin_vis_iso_loc(
    float *Qp3d,
    float *Qs3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    int grid_type,
    std::vector<float> &xvec, 
    std::vector<float> &yvec, 
    std::vector<float> &zvec,
    float *bin_Qp,
    float *bin_Qs );

#endif  // MEDIA_BIN2MODEL 
