#ifndef __MEDIA_GRID2MODEL__
#define __MEDIA_GRID2MODEL__

#include "media_geometry3d.hpp"
#include "media_read_interface_file.hpp"

int AssignGridMediaPara2Point(
    int ix, int iy, int iz, 
    Point3 A, 
    int NI, 
    Interfaces *interfaces,
    float &vp, 
    float &vs,
    float &rho);

int LayerNumberAtPoint(
    Point3 A, 
    int NL,
    std::vector<int> NGz, 
    Interfaces *interfaces);

void iso_grid_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NL, 
    std::vector <int> NGz,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d);

void iso_grid_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NL, 
    std::vector <int> NGz,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d);

#ifdef __cplusplus
extern "C" {
#endif
void media_el_iso_grid2model(
    float *lam3d,
    float *mu3d,
    float *rho3d,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    const char *grid_file,
    const char *equivalent_medium_method); 
#ifdef __cplusplus
}
#endif

#endif // MEIDA_GRID2MODEL
