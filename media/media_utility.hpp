#ifndef __MEDIA_UTILITY__
#define __MEDIA_UTILITY__
#include "media_geometry3d.hpp"


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

    float *altitude = nullptr;
    float *vp       = nullptr;
    float *vp_grad  = nullptr;
    float *vs       = nullptr;
    float *vs_grad  = nullptr;
    float *rho      = nullptr;
    float *rho_grad = nullptr;
};

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