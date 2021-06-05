
#ifndef __PRE_LAYER2MODEL__
#define __PRE_LAYER2MODEL__

#include "media_geometry3d.hpp"
//using namespace std;
//#define MEDIA_LOC 0
//#define MEDIA_HAR 1
//#define MEDIA_TTI 2

// number used for divided the mesh
#define NG 8

// Interface information from the interface file
struct Interfaces{
    // all the interfaces are given by the same interface_file mesh, 
    //  so maybe we can store the interface mesh info out of the struct  
    size_t  NX = 0;
    size_t  NY = 0;
    float   DX = 0.0;
    float   DY = 0.0;
    float MINX = 0.0;
    float MINY = 0.0;

    float *depth    = nullptr;
    float *vp       = nullptr;
    float *vp_grad  = nullptr;
    float *vs       = nullptr;
    float *vs_grad  = nullptr;
    float *rho      = nullptr;
    float *rho_grad = nullptr;
};

int NumOfValues(std::vector<int> v, int NI);

int AssignMediaPara2Point(
    int ix, int iy, int iz, 
    Point3 A, 
    int NI, 
    Interfaces *interfaces,
    float &vp, 
    float &vs,
    float &rho);

void GenerateHalfGrid(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    float *halfGridx,
    float *halfGridy,
    float *halfGridz);


Point3 *MeshSubdivide(Mesh3 M);

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
//#include <iostream>
// medium parameterization
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
    size_t siz_line,
    size_t siz_slice,
    size_t siz_volume, 
    const char *interfaces_file,
    const char *equivalent_medium_method); 
#ifdef __cplusplus
}
#endif
 

#endif /* __PRE_LAYER2MODEL__ */
