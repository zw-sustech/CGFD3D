#ifndef MEDIA_DISCRETE_MODEL_H
#define MEDIA_DISCRETE_MODEL_H

// for C code call
#define MEDIA_USE_CART 1
#define MEDIA_USE_VMAP 2
#define MEDIA_USE_CURV 3

/*------------ bin2model --------------*/
int media_bin2model_el_iso(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    float xmin, float xmax,
    float ymin, float ymax,
    int grid_type,
    int  *bin_order,    // eg, [2, 0, 1]=[z, x, y] 0:x, 1:y, 2:z
    int  *bin_size,     // [ndim1, ndim2, ndim3],
    float  *bin_spacing,  // [dh1, dh2, dh3],
    float  *bin_origin,   // [h0_1, h0_2, h0_3],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs  );
//---- 0. viscoelastic for Qp,Qs
int media_bin2model_vis_iso(
    float *Qp3d,
    float *Qs3d,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    float xmin, float xmax,
    float ymin, float ymax,
    int grid_type,
    int  *bin_order,    // eg, [2, 0, 1]=[z, x, y] 0:x, 1:y, 2:z
    int  *bin_size,     // [ndim1, ndim2, ndim3],
    float  *bin_spacing,  // [dh1, dh2, dh3],
    float  *bin_origin,   // [h0_1, h0_2, h0_3],
    const char *bin_file_Qp,
    const char *bin_file_Qs  );

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
                             const char *average_method,
                             int myid);

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
        const char *equivalent_medium_method,
        int myid);

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
        const char *equivalent_medium_method,
        int myid);

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
        const char *equivalent_medium_method,
        int myid);

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
        const char *equivalent_medium_method,
        int myid); 


/*-------------- grid2model -------------*/

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
    const char *equivalent_medium_method,
    int myid);

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
    const char *equivalent_medium_method, 
    int myid);

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
    const char *equivalent_medium_method,
    int myid);

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
    const char *equivalent_medium_method, 
    int myid);

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
    const char *equivalent_medium_method,
    int myid); 

//================== bin2model ===================
int media_bin2model_el_iso(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    size_t nx, size_t ny, size_t nz,
    float xmin, float xmax,
    float ymin, float ymax,
    int grid_type,
    int *bin_order,    // eg, [2, 0, 1]=[z, x, y] 0:x, 1:y, 2:z
    int *bin_size,     // [ndim1, ndim2, ndim3],
    float  *bin_spacing,  // [dh1, dh2, dh3],
    float  *bin_origin,   // [h0_1, h0_2, h0_3],
    const char *bin_file_rho,
    const char *bin_file_vp,
    const char *bin_file_vs );

#endif
