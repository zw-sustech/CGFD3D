#ifndef __MEDIA_UTILITY__
#define __MEDIA_UTILITY__
#include "media_geometry3d.hpp"

// for media type
#define ONE_COMPONENT       0 /* var */
#define ELASTIC_ISOTROPIC   1 /* vp, vs, rho*/
#define ELASTIC_VTI_PREM    2 /* vph, vpv, vsh, vsv, rho, eta */
#define ELASTIC_VTI_THOMSEN 3 /* vpv, vsv, epsilon, delta, gamma, rho */
#define ELASTIC_VTI_CIJ     4 /* c11 c33 c44 c66 c13 rho */
#define ELASTIC_TTI_THOMSEN 5 /* vp0, vs0, epsilon, delta, gamma, rho, azimuth, dip */
#define ELASTIC_TTI_BOND    6 /* c11 c33 c44 c66 c13 rho azimuth dip */
#define ELASTIC_TTI_CIJ     7 /* c11 c12 c13 c14 c15 c16 c22 c23 c24 c25 c26 c33 c34 c35 c36 c44 c45 c46 c55 c56 c66 rho*/
#define ACOUSTIC_ISOTROPIC  8 /* vp, rho */

/* 
 * number used for divided the mesh, 
 * for equivalent medium parameterization
 */
#define NG 8

/* Interface information from the interface file (different media type) */
struct inter_t{
    /* all the interfaces are given by the same interface_file mesh. */
    size_t  NI = 0; // number of interfaces
    size_t  NX = 0;
    size_t  NY = 0;
    float   DX = FLT_MAX;
    float   DY = FLT_MAX;
    float MINX = FLT_MAX;
    float MINY = FLT_MAX;

    // ni*slice
    float *elevation= nullptr;

    // for acoustic + elastic
    float *vp       = nullptr;
    float *rho      = nullptr;
    float *vp_grad  = nullptr;
    float *rho_grad = nullptr;
    float *vp_pow   = nullptr;
    float *rho_pow  = nullptr;
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

    // for vti_prem
    float *vph      = nullptr;
    float *vpv      = nullptr;
    float *vsh      = nullptr;
    float *vsv      = nullptr;
    float *eta      = nullptr;
    float *vph_grad = nullptr;
    float *vpv_grad = nullptr;
    float *vsh_grad = nullptr;
    float *vsv_grad = nullptr;
    float *eta_grad = nullptr;
    float *vph_pow  = nullptr;
    float *vpv_pow  = nullptr;
    float *vsh_pow  = nullptr;
    float *vsv_pow  = nullptr;
    float *eta_pow  = nullptr;

    // for one component
    float *var      = nullptr;
    float *var_grad = nullptr;
    float *var_pow  = nullptr;

    // for anisotropy c_ij, call one component
    float *c11 = nullptr;
    float *c12 = nullptr;
    float *c13 = nullptr;
    float *c14 = nullptr;
    float *c15 = nullptr;
    float *c16 = nullptr;
    float *c22 = nullptr;
    float *c23 = nullptr;
    float *c24 = nullptr;
    float *c25 = nullptr;
    float *c26 = nullptr;
    float *c33 = nullptr;
    float *c34 = nullptr;
    float *c35 = nullptr;
    float *c36 = nullptr;
    float *c44 = nullptr;
    float *c45 = nullptr;
    float *c46 = nullptr;
    float *c55 = nullptr;
    float *c56 = nullptr;
    float *c66 = nullptr;
    float *c11_grad = nullptr;
    float *c12_grad = nullptr;
    float *c13_grad = nullptr;
    float *c14_grad = nullptr;
    float *c15_grad = nullptr;
    float *c16_grad = nullptr;
    float *c22_grad = nullptr;
    float *c23_grad = nullptr;
    float *c24_grad = nullptr;
    float *c25_grad = nullptr;
    float *c26_grad = nullptr;
    float *c33_grad = nullptr;
    float *c34_grad = nullptr;
    float *c35_grad = nullptr;
    float *c36_grad = nullptr;
    float *c44_grad = nullptr;
    float *c45_grad = nullptr;
    float *c46_grad = nullptr;
    float *c55_grad = nullptr;
    float *c56_grad = nullptr;
    float *c66_grad = nullptr;
    float *c11_pow = nullptr;
    float *c12_pow = nullptr;
    float *c13_pow = nullptr;
    float *c14_pow = nullptr;
    float *c15_pow = nullptr;
    float *c16_pow = nullptr;
    float *c22_pow = nullptr;
    float *c23_pow = nullptr;
    float *c24_pow = nullptr;
    float *c25_pow = nullptr;
    float *c26_pow = nullptr;
    float *c33_pow = nullptr;
    float *c34_pow = nullptr;
    float *c35_pow = nullptr;
    float *c36_pow = nullptr;
    float *c44_pow = nullptr;
    float *c45_pow = nullptr;
    float *c46_pow = nullptr;
    float *c55_pow = nullptr;
    float *c56_pow = nullptr;
    float *c66_pow = nullptr;
/*
    ~inter_t() {
        // ni*slice
        if (elevation != nullptr) delete [] elevation;
    
        // for acoustic + elastic
        if (vp       != nullptr) delete [] vp      ;
        if (rho      != nullptr) delete [] rho     ;
        if (vp_grad  != nullptr) delete [] vp_grad ;
        if (rho_grad != nullptr) delete [] rho_grad;
        if (vp_pow   != nullptr) delete [] vp_pow  ;
        if (rho_pow  != nullptr) delete [] rho_pow ;
        if (vs       != nullptr) delete [] vs      ;
        if (vs_grad  != nullptr) delete [] vs_grad ;
        if (vs_pow   != nullptr) delete [] vs_pow  ;
    
        // for anisotropy (thomsen)
        if (epsilon      != nullptr) delete [] epsilon     ;
        if (epsilon_grad != nullptr) delete [] epsilon_grad;
        if (epsilon_pow  != nullptr) delete [] epsilon_pow ;
        if (delta        != nullptr) delete [] delta       ;
        if (delta_grad   != nullptr) delete [] delta_grad  ;
        if (delta_pow    != nullptr) delete [] delta_pow   ;
        if (gamma        != nullptr) delete [] gamma       ;
        if (gamma_grad   != nullptr) delete [] gamma_grad  ;
        if (gamma_pow    != nullptr) delete [] gamma_pow   ;
        if (azimuth      != nullptr) delete [] azimuth     ;
        if (azimuth_grad != nullptr) delete [] azimuth_grad;
        if (azimuth_pow  != nullptr) delete [] azimuth_pow ;
        if (dip          != nullptr) delete [] dip         ;
        if (dip_grad     != nullptr) delete [] dip_grad    ;
        if (dip_pow      != nullptr) delete [] dip_pow     ;
    
        // for vti_prem
        if (vph      != nullptr) delete [] vph     ;
        if (vpv      != nullptr) delete [] vpv     ;
        if (vsh      != nullptr) delete [] vsh     ;
        if (vsv      != nullptr) delete [] vsv     ;
        if (eta      != nullptr) delete [] eta     ;
        if (vph_grad != nullptr) delete [] vph_grad;
        if (vpv_grad != nullptr) delete [] vpv_grad;
        if (vsh_grad != nullptr) delete [] vsh_grad;
        if (vsv_grad != nullptr) delete [] vsv_grad;
        if (eta_grad != nullptr) delete [] eta_grad;
        if (vph_pow  != nullptr) delete [] vph_pow ;
        if (vpv_pow  != nullptr) delete [] vpv_pow ;
        if (vsh_pow  != nullptr) delete [] vsh_pow ;
        if (vsv_pow  != nullptr) delete [] vsv_pow ;
        if (eta_pow  != nullptr) delete [] eta_pow ;
    
        // for one component
        if (var      != nullptr) delete [] var     ; 
        if (var_grad != nullptr) delete [] var_grad;
        if (var_pow  != nullptr) delete [] var_pow ;
    
        // for anisotropy c_ij, call one component
        if (c11 != nullptr) delete [] c11;
        if (c12 != nullptr) delete [] c12;
        if (c13 != nullptr) delete [] c13;
        if (c14 != nullptr) delete [] c14;
        if (c15 != nullptr) delete [] c15;
        if (c16 != nullptr) delete [] c16;
        if (c22 != nullptr) delete [] c22;
        if (c23 != nullptr) delete [] c23;
        if (c24 != nullptr) delete [] c24;
        if (c25 != nullptr) delete [] c25;
        if (c26 != nullptr) delete [] c26;
        if (c33 != nullptr) delete [] c33;
        if (c34 != nullptr) delete [] c34;
        if (c35 != nullptr) delete [] c35;
        if (c36 != nullptr) delete [] c36;
        if (c44 != nullptr) delete [] c44;
        if (c45 != nullptr) delete [] c45;
        if (c46 != nullptr) delete [] c46;
        if (c55 != nullptr) delete [] c55;
        if (c56 != nullptr) delete [] c56;
        if (c66 != nullptr) delete [] c66;
        if (c11_gard != nullptr) delete [] c11_gard;
        if (c12_gard != nullptr) delete [] c12_gard;
        if (c13_gard != nullptr) delete [] c13_gard;
        if (c14_gard != nullptr) delete [] c14_gard;
        if (c15_gard != nullptr) delete [] c15_gard;
        if (c16_gard != nullptr) delete [] c16_gard;
        if (c22_gard != nullptr) delete [] c22_gard;
        if (c23_gard != nullptr) delete [] c23_gard;
        if (c24_gard != nullptr) delete [] c24_gard;
        if (c25_gard != nullptr) delete [] c25_gard;
        if (c26_gard != nullptr) delete [] c26_gard;
        if (c33_gard != nullptr) delete [] c33_gard;
        if (c34_gard != nullptr) delete [] c34_gard;
        if (c35_gard != nullptr) delete [] c35_gard;
        if (c36_gard != nullptr) delete [] c36_gard;
        if (c44_gard != nullptr) delete [] c44_gard;
        if (c45_gard != nullptr) delete [] c45_gard;
        if (c46_gard != nullptr) delete [] c46_gard;
        if (c55_gard != nullptr) delete [] c55_gard;
        if (c56_gard != nullptr) delete [] c56_gard;
        if (c66_gard != nullptr) delete [] c66_gard;
        if (c11_pow != nullptr) delete [] c11_pow;
        if (c12_pow != nullptr) delete [] c12_pow;
        if (c13_pow != nullptr) delete [] c13_pow;
        if (c14_pow != nullptr) delete [] c14_pow;
        if (c15_pow != nullptr) delete [] c15_pow;
        if (c16_pow != nullptr) delete [] c16_pow;
        if (c22_pow != nullptr) delete [] c22_pow;
        if (c23_pow != nullptr) delete [] c23_pow;
        if (c24_pow != nullptr) delete [] c24_pow;
        if (c25_pow != nullptr) delete [] c25_pow;
        if (c26_pow != nullptr) delete [] c26_pow;
        if (c33_pow != nullptr) delete [] c33_pow;
        if (c34_pow != nullptr) delete [] c34_pow;
        if (c35_pow != nullptr) delete [] c35_pow;
        if (c36_pow != nullptr) delete [] c36_pow;
        if (c44_pow != nullptr) delete [] c44_pow;
        if (c45_pow != nullptr) delete [] c45_pow;
        if (c46_pow != nullptr) delete [] c46_pow;
        if (c55_pow != nullptr) delete [] c55_pow;
        if (c56_pow != nullptr) delete [] c56_pow;
        if (c66_pow != nullptr) delete [] c66_pow;        
    }
*/
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

/* Reference: Bond, 1943, The Mathematics of the Physical Properties of Crystals */
/* theta: dip, phi: azimuth*/
void BondTransform(float c11, float c12, float c13, float c14, float c15, float c16, 
                   float c22, float c23, float c24, float c25, float c26, float c33,
                   float c34, float c35, float c36, float c44, float c45, float c46,
                   float c55, float c56, float c66, float theta, float phi, 
                   float &c11_tti, float &c12_tti, float &c13_tti, float &c14_tti, float &c15_tti, float &c16_tti, 
                   float &c22_tti, float &c23_tti, float &c24_tti, float &c25_tti, float &c26_tti, float &c33_tti,
                   float &c34_tti, float &c35_tti, float &c36_tti, float &c44_tti, float &c45_tti, float &c46_tti,
                   float &c55_tti, float &c56_tti, float &c66_tti); 


/*------ find point and interpolation -------*/
int findLastGreaterEqualIndex(
    float value, 
    std::vector<float> &x);

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

/*---- matrix: just for Bond transform ----*/
template <typename T>
class Matrix;

template <typename T>
std::ostream &operator<<(std::ostream& out, const Matrix<T> &mat);

template <typename T>
class Matrix{
private:
    int row;
    int col;
    T* p;
public:
    Matrix(int r, int c);
    Matrix(int r, int c, T* initval);
    Matrix(const Matrix<T> &mat); // copy
    ~Matrix();
    Matrix<T> operator*(const Matrix &mat);
    Matrix<T> operator=(const Matrix &mat);
    Matrix<T> transpose();
    T &operator()(int i, int j)const;
    friend std::ostream &operator<< <T>(std::ostream& out, const Matrix<T> &mat);
};
/*-----------------------------------------*/

#endif
