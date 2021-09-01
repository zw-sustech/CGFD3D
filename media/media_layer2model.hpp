
#ifndef __MEDIA_LAYER2MODEL__
#define __MEDIA_LAYER2MODEL__

#include "media_geometry3d.hpp"
#include "media_utility.hpp"

int AssignLayerMediaPara2Point_onecmp(
    int ix, int iy, int iz, 
    Point3 A,  
    inter_t interfaces,
    float &var);

int AssignLayerMediaPara2Point_el_iso(
    int ix, int iy, int iz, 
    Point3 A, 
    inter_t interfaces,
    float &vp, float &vs, float &rho);

// number used for divided the mesh
int MediaNumberAtPoint(
    Point3 A, 
    inter_t interfaces);

void parametrization_oncmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, // 1: xyz1d; 2: xy1d, z3d; 3: xyz3d.
    inter_t interfaces,
    float *var3d);

void parametrization_el_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, // 1: xyz1d; 2: xy1d, z3d; 3: xyz3d.
    inter_t interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d);

void parametrization_el_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d);

//void md_el_iso_rho_to_slow(float *rho, size_t siz_volume);

#ifdef __cplusplus
extern "C" {
#endif
int media_layer2model_curv_el_iso(
        float *lam3d,
        float *mu3d,
        float *rho3d,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        int nx,
        int ny,
        int nz,
        const char *in_rho_file,
        const char *in_vp_file,
        const char *in_vs_file,
        const char *equivalent_medium_method); 
int media_layer2model_curv_onecmp(
        float *var3d,
        const float *x3d, // nx*ny*nz
        const float *y3d,
        const float *z3d,
        int nx,
        int ny, 
        int nz,
        const char *in_var_file,
        const char *average_method);
#ifdef __cplusplus
}
#endif
 

#endif /* __PRE_LAYER2MODEL__ */
