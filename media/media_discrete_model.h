#ifndef MEDIA_DISCRETE_MODEL_H
#define MEDIA_DISCRETE_MODEL_H

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
                             size_t nx,
                             size_t ny, 
                             size_t nz,
                             int grid_type, 
                             const char *in_var_file,
                             const char *average_method);

//---- 1. elastic isotropic
int media_layer2model_ac_iso(
        float *rho3d,
        float *kappa3d,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file,
        const char *equivalent_medium_method);

//----  2. elastic isotropic
int media_layer2model_el_iso(
        float *lam3d,
        float *mu3d,
        float *rho3d,
        const float *x3d, 
        const float *y3d, 
        const float *z3d, 
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file,
        const char *equivalent_medium_method);

//--- 3. elastic vti
int media_layer2model_el_vti(
        float *rho,
        float *c11,
        float *c33,
        float *c55,
        float *c66,
        float *c13,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file, 
        const char *equivalent_medium_method);

//--- 4. elastic anisotropic/TTI
int media_layer2model_el_aniso(
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
        size_t nx,
        size_t ny,
        size_t nz,
        int grid_type, 
        const char *in_3lay_file,
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

#endif
