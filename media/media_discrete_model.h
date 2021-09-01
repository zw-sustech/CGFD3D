// for C code call
int media_layer2model_curv_el_iso(
        float *lam3d,
        float *mu3d,
        float *rho3d,
        const float *x3d,
        const float *y3d,
        const float *z3d,
        int nx,
        int ny,
        int nz,
        const char *in_rho_file,
        const char *in_vp_file,
        const char *in_vs_file,
        const char *equivalent_medium_method);

int media_layer2model_curv_onecmp(
        float *var3d,
        const float *x3d, // nx*ny*nz
        const float *y3d,
        const float *z3d,
        int nx,
        int ny, 
        int nz,
        const char *in_var_file,
        const char *average_method); 

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
