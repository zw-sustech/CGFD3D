#ifndef _MEDIA_GRID2MODEL_
#define _MEDIA_GRID2MODEL_

#include "media_geometry3d.hpp"
#include "media_read_file.hpp"


#ifdef __cplusplus
extern "C" {
#endif
// --- 0. one component
int media_grid2model_onecmp(
    float *var3d,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 1. acoustic isotropic
int media_grid2model_ac_iso(
    float *rho3d,
    float *kappa3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 2. elastic isotropic
int media_grid2model_el_iso(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

//--- 3. elastic vti
int media_grid2model_el_vti(
    float *rho,
    float *c11, 
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

int media_grid2model_el_aniso(
    float *rho,
    float *c11, float *c12, float *c13,
    float *c14, float *c15, float *c16,
    float *c22, float *c23, float *c24, 
    float *c25, float *c26, float *c33,
    float *c34, float *c35, float *c36,
    float *c44, float *c45, float *c46,
    float *c55, float *c56, float *c66,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method);

#ifdef __cplusplus
}
#endif

int AssignGridMediaPara2Point(
    size_t ix, size_t iy, size_t iz, 
    Point3 A, 
    inter_t &interfaces,
    int media_type,
    std::vector<float> &var);

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for grid2model
void CalPointValue_grid(int media_type, 
                   inter_t &interfaces,
                   size_t slice, 
                   std::vector<float> &xvec,  /* interface mesh */
                   std::vector<float> &yvec,
                   int NI, 
                   Point3 &A,
                   std::vector<float> &elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var);


//- 0. assign the parameter directly (use the local values): one component 
void parametrization_grid_onecmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *var3d);

//- 1. assign the parameter directly (use the local values): isotropic, acoustic 
void parametrization_grid_ac_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *kappa);

//- 2. assign the parameter directly (use the local values): elastic isotropic 
void parametrization_grid_el_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *rho3d, 
    float *lam3d,
    float *mu3d);

//- 3. Assign the parameter directly (use the local values): elastic vti
void parametrization_grid_el_vti_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho);

//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_grid_el_aniso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    float *rho);

//======================= averaging/equivalent medium method =========================

/* 
 * For half-grid point, mark the interface number.  
 */
int LayerNumberAtPoint(
    Point3 A, 
    int NL,
    std::vector<int> &NGz, 
    inter_t &interfaces);

void MarkLayerNumber(
    int grid_type, 
    float *Hx, float *Hy, float *Hz,
    size_t nx, size_t ny, size_t nz,
    int NL, std::vector<int> &NGz,
    int *MaterNum,
    inter_t &interfaces);

//- 0.1 assign the parameter by volume harmonic averaging
//- one component
int parametrization_grid_onecmp_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *var3d);

//- 0.2 assign the parameter by volume arithmetic averaging
//- one component
int parametrization_grid_onecmp_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *var3d);

//- 1.1 assign the parameter by volume arithmetic and harmonic averaging method
//- acoustic isotropic 
int parametrization_grid_ac_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *kappa);

//- 1.2 assign the parameter by volume arithmetic averaging method
//- acoustic isotropic 
int parametrization_grid_ac_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *kappa);

//- 2.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic isotropic
// Moczo et al., 2002 
int parametrization_grid_el_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *lam3d,
    float *mu3d );

//- 2.1 assign the parameter by volume arithmetic averaging method
//- elastic isotropic
int parametrization_grid_el_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *lam3d,
    float *mu3d);

//- 3.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic vti
int parametrization_grid_el_vti_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho);

//- 3.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic vti
int parametrization_grid_el_vti_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho);

//- 4.1 assign the parameter by volume arithmetic averaging method
//- elastic tti/anisotropic
int parametrization_grid_el_aniso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
     float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    float *rho);

//- 4.2 assign the parameter by volume arithmetic averaging method
//- elastic tti
int parametrization_grid_el_aniso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
     float *c11,
    float *c12,
    float *c13,
    float *c14,
    float *c15,
    float *c16,
    float *c22,
    float *c23,
    float *c24,
    float *c25,
    float *c26,
    float *c33,
    float *c34,
    float *c35,
    float *c36,
    float *c44,
    float *c45,
    float *c46,
    float *c55,
    float *c56,
    float *c66,
    float *rho);

#endif // MEIDA_GRID2MODEL
