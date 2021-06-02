// for C code call
void pre_el_iso_layer2model(
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
