#ifndef __MEDIA_UTILITY__
#define __MEDIA_UTILITY__
#include "media_geometry3d.hpp"

// number used for divided the mesh
#define NG 8

// Interface information from the interface file
struct inter_t{
    // all the interfaces are given by the same interface_file mesh.
    size_t  NI = 0; // number of interfaces
    size_t  NX = 0;
    size_t  NY = 0;
    float   DX = FLT_MAX;
    float   DY = FLT_MAX;
    float MINX = FLT_MAX;
    float MINY = FLT_MAX;

    // ni*slice

    float *elevation= nullptr;

    // for aco + el + anisotropy (vpv, rho)
    float *vp       = nullptr;
    float *vp_grad  = nullptr;
    float *vp_pow   = nullptr;
    float *rho      = nullptr;
    float *rho_grad = nullptr;
    float *rho_pow  = nullptr;
    // el + amosotropy_thomsen (vsv)
    float *vs       = nullptr;
    float *vs_grad  = nullptr;
    float *vs_pow   = nullptr;

    // for anisotropy (thomsen)
    float *epsilon      = nullptr; 
    float *epsilon_grad = nullptr; 
    float *epsilon_pow  = nullptr; 
    float *delta        = nullptr;
    float *delta_grad   = nullptr;
    float *delta_pow    = nullptr;
    float *gamma        = nullptr;
    float *gamma_grad   = nullptr;
    float *gamma_pow    = nullptr;
    float *azimuth      = nullptr;
    float *azimuth_grad = nullptr;
    float *azimuth_pow  = nullptr;
    float *dip          = nullptr;
    float *dip_grad     = nullptr;
    float *dip_pow      = nullptr;

    // for one component
    float *var      = nullptr;
    float *var_grad = nullptr;
    float *var_pow  = nullptr;

    // for anisotropy c_ij, call one component
};

bool isEqual(float a, float b);

int NumOfValues(std::vector<int> v, int NI);

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

int findNearestNeighborIndex(
    float value, std::vector<float> &x);

int findNearestGreaterIndex(
    float value, std::vector<float> &x);

float BilinearInterpolation(
    std::vector<float> &x, 
    std::vector<float> &y, 
    float *v,
    float xq,
    float yq );

#endif
