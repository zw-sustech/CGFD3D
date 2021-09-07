#ifndef __MEDIA_LAYER2MODEL__
#define __MEDIA_LAYER2MODEL__

#include "media_geometry3d.hpp"
#include "media_utility.hpp"


void PrintIsPointOutOfInterfaceRange(Point3 A, 
    int ix, int iy, int iz, 
    float MINX, float MAXX, float MINY, float MAXY);

int AssignLayerMediaPara2Point(
    int ix, int iy, int iz,         /* To print error messages */ 
    Point3 A,  
    inter_t interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var);


/* Calculate the value of the point for different media type (to avoid multiple geometric calculations) */
void CalPointValue(int media_type, 
                   inter_t interfaces,
                   size_t slice, 
                   std::vector<float> xvec,  /* interface mesh */
                   std::vector<float> yvec,
                   Point3 A,
                   std::vector<float> elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var);

/* 
 * For half-grid point, marked the materials number.  
 */
int MediaNumberAtPoint(
        Point3 A,  
        inter_t interfaces);

/* 0. assign the parameter directly (use the local values): one component */
void parametrization_oncmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t interfaces,
    float *var3d);

/* 1. assign the parameter directly (use the local values): isotropic */
void parametrization_el_iso_loc(
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

/*  2. Assign the parameter directly (use the local values): vti-prem */
/* reference: Dziewonski and Anderson, 1981, Preliminary reference Earth model */
void parametrization_el_vti_prem_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *c11_3d,
    float *c33_3d,
    float *c44_3d,
    float *c66_3d,
    float *c13_3d,
    float *rho_3d);

/* 3. Assign the parameter directly (use the local values): vti-thomsen */
/* reference: Thomsen, 1986, Weak elastic anisotropy */
void parametrization_el_vti_thomsen_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *c11_3d,
    float *c33_3d,
    float *c44_3d,
    float *c66_3d,
    float *c13_3d,
    float *rho_3d);

/* 4. Assign the parameter directly (use the local values): vti-cij */
void parametrization_el_vti_cij_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *c11_3d,
    float *c33_3d,
    float *c44_3d,
    float *c66_3d,
    float *c13_3d,
    float *rho_3d);

/* 5. Assign the parameter directly (use the local values): tti-thomsen */
/* reference: Thomsen, 1986, Weak elastic anisotropy */
void parametrization_el_tti_thomsen_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *c11_3d,
    float *c12_3d,
    float *c13_3d,
    float *c14_3d,
    float *c15_3d,
    float *c16_3d,
    float *c22_3d,
    float *c23_3d,
    float *c24_3d,
    float *c25_3d,
    float *c26_3d,
    float *c33_3d,
    float *c34_3d,
    float *c35_3d,
    float *c36_3d,
    float *c44_3d,
    float *c45_3d,
    float *c46_3d,
    float *c55_3d,
    float *c56_3d,
    float *c66_3d,
    float *rho_3d);

/* 6. Assign the parameter directly (use the local values): tti-bond */
void parametrization_el_tti_bond_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *c11_3d,
    float *c12_3d,
    float *c13_3d,
    float *c14_3d,
    float *c15_3d,
    float *c16_3d,
    float *c22_3d,
    float *c23_3d,
    float *c24_3d,
    float *c25_3d,
    float *c26_3d,
    float *c33_3d,
    float *c34_3d,
    float *c35_3d,
    float *c36_3d,
    float *c44_3d,
    float *c45_3d,
    float *c46_3d,
    float *c55_3d,
    float *c56_3d,
    float *c66_3d,
    float *rho_3d);

/* 7. Assign the parameter directly (use the local values): tti-cij */
void parametrization_el_tti_cij_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *c11_3d,
    float *c12_3d,
    float *c13_3d,
    float *c14_3d,
    float *c15_3d,
    float *c16_3d,
    float *c22_3d,
    float *c23_3d,
    float *c24_3d,
    float *c25_3d,
    float *c26_3d,
    float *c33_3d,
    float *c34_3d,
    float *c35_3d,
    float *c36_3d,
    float *c44_3d,
    float *c45_3d,
    float *c46_3d,
    float *c55_3d,
    float *c56_3d,
    float *c66_3d,
    float *rho_3d);

/* 7. assign the parameter directly (use the local values): isotropic, acoustic */
void parametrization_ac_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t interfaces,
    float *kappa,
    float *rho3d);

/* har: elastic isotropic, every grid type */
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

/* har: elastic vti, every grid type */
/* har: elastic tti, every grid type */
/* har: acoustic isotropic, every grid type */


#ifdef __cplusplus
extern "C" {
#endif
//---- 0. one component
int media_layer2model_onecmp(float *var3d,
                             const float *x3d, 
                             const float *y3d, 
                             const float *z3d, 
                             int nx,
                             int ny, 
                             int nz,
                             int grid_type,
                             const char *in_var_file,
                             const char *average_method);

//---- 1. elastic isotropic
int media_layer2model_el_iso(
        float *lam3d,
        float *mu3d,
        float *rho3d,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        int nx,
        int ny,
        int nz,
        int grid_type,
        const char *in_rho_file,
        const char *in_vp_file,
        const char *in_vs_file,
        const char *equivalent_medium_method);


//--- 2. elastic vti_prem
int media_layer2model_el_vti_prem(
        float *rho,
        float *c11,
        float *c33,
        float *c44,
        float *c66,
        float *c13,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        int nx,
        int ny,
        int nz,
        int grid_type,
        const char *in_rho_file,
        const char *in_vph_file,
        const char *in_vpv_file,
        const char *in_vsh_file,
        const char *in_vsv_file,
        const char *in_eta_file,
        const char *equivalent_medium_method);

//--- 3. elastic vti_thomsen
int media_layer2model_vti_thomsen(
        float *rho,
        float *c11,
        float *c33,
        float *c44,
        float *c66,
        float *c13,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        int nx,
        int ny,
        int nz,
        int grid_type, 
        const char *in_rho_file,
        const char *in_vpv_file,
        const char *in_vsv_file,
        const char *in_epsilon_file,
        const char *in_delta_file,
        const char *in_gamma_file,
        const char *equivalent_medium_method);

//--- 4. elastic vti_cij
int media_layer2model_vti_cij(
        float *rho,
        float *c11,
        float *c33,
        float *c44,
        float *c66,
        float *c13,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        int nx,
        int ny,
        int nz,
        int grid_type, 
        const char *in_rho_file,
        const char *in_c11_file,
        const char *in_c33_file,
        const char *in_c44_file,
        const char *in_c66_file,
        const char *in_c13_file,
        const char *equivalent_medium_method);


//--- 5. elastic tti_thomsen
int media_layer2model_tti_thomsen(
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
        int grid_type,
        const char *in_rho_file,
        const char *in_vpv_file,
        const char *in_vsv_file,
        const char *in_epsilon_file,
        const char *in_delta_file,
        const char *in_gamma_file,
        const char *in_azimuth_file,
        const char *in_dip_file,
        const char *equivalent_medium_method);


//--- 6. elastic tti_bond
int media_layer2model_tti_bond(
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
        int grid_type, 
        const char *in_rho_file,
        const char *in_c11_file,
        const char *in_c33_file,
        const char *in_c44_file,
        const char *in_c66_file,
        const char *in_c13_file,
        const char *in_azimuth_file,
        const char *in_dip_file,
        const char *equivalent_medium_method);

//--- 7. elastic tti_cij
int media_layer2model_tti_cij(
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
        int grid_type,
        char *in_rho_file,
        char *in_c11_file,
        char *in_c12_file,
        char *in_c13_file,
        char *in_c14_file,
        char *in_c15_file,
        char *in_c16_file,
        char *in_c22_file,
        char *in_c23_file,
        char *in_c24_file,
        char *in_c25_file,
        char *in_c26_file,
        char *in_c33_file,
        char *in_c34_file,
        char *in_c35_file,
        char *in_c36_file,
        char *in_c44_file,
        char *in_c45_file,
        char *in_c46_file,
        char *in_c55_file,
        char *in_c56_file,
        char *in_c66_file,
        const char *equivalent_medium_method);

//---- 8. acoustic isotropic
// grid type 1
int media_layer2model_ac_iso(
        float *rho3d,
        float *kappa3d,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        int nx,
        int ny,
        int nz,
        int grid_type,
        const char *in_rho_file,
        const char *in_vp_file,
        const char *equivalent_medium_method);

#ifdef __cplusplus
}
#endif
 

#endif /* __PRE_LAYER2MODEL__ */
