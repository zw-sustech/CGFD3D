// for C code call
#define MEDIA_USE_CART 1
#define MEDIA_USE_VMAP 2
#define MEDIA_USE_CURV 3

/*--------------------------- layer2model --------------------- */
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

//---- 7. acoustic isotropic
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

/*-------------- grid2model -------------*/
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
