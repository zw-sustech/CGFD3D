
#ifndef __MEDIA_LAYER2MODEL__
#define __MEDIA_LAYER2MODEL__

#include "media_geometry3d.hpp"
#include "media_utility.hpp"


// number used for divided the mesh
int AssignMediaPara2Point(
    int ix, int iy, int iz, 
    Point3 A, 
    int NI, 
    Interfaces *interfaces,
    float &vp, 
    float &vs,
    float &rho);

int MediaNumberAtPoint(
    Point3 A, 
    int NI, 
    Interfaces *interface);

void isotropic_loc(
    //grid info
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NI,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d );

void isotropic_har(
    //grid info
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NI,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d);

//void md_el_iso_rho_to_slow(float *rho, size_t siz_volume);

#ifdef __cplusplus
extern "C" {
#endif
void media_el_iso_layer2model(
    float *lam3d,
    float *mu3d,
    float *rho3d,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    const char *interfaces_file,
    const char *equivalent_medium_method); 
#ifdef __cplusplus
}
#endif
 

#endif /* __PRE_LAYER2MODEL__ */
