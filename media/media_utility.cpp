#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include <string.h>
//#include <Eigen>
#include "media_utility.hpp"

bool isEqual(float a, float b) {
    return abs(a-b) < FLT_EPSILON;
}

void GenerateHalfGrid(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    float *halfGridx,
    float *halfGridy,
    float *halfGridz) 
{
    /* 
     * half grid point is inside the grid, 
     * so the number of grid point in one direction is np-1.
     */
    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;

    for (size_t k = 0; k < nz-1; k++) {
        for (size_t j = 0; j < ny-1; j++) {
            for (size_t i = 0; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice;

                /* clockwise: conducive to debug */
                halfGridx[indx] = (Gridx[indx]       + Gridx[indx+1]                  + 
                    Gridx[indx+siz_line+1]           + Gridx[indx+siz_line]           + 
                    Gridx[indx+siz_slice]            + Gridx[indx+siz_slice+1]        + 
                    Gridx[indx+siz_slice+siz_line+1] + Gridx[indx+siz_line+siz_slice]  )/8.0;

                halfGridy[indx] = (Gridy[indx]       + Gridy[indx+1]                  + 
                    Gridy[indx+siz_line+1]           + Gridy[indx+siz_line]           + 
                    Gridy[indx+siz_slice]            + Gridy[indx+siz_slice+1]        + 
                    Gridy[indx+siz_slice+siz_line+1] + Gridy[indx+siz_line+siz_slice]  )/8.0;

                halfGridz[indx] = (Gridz[indx]       + Gridz[indx+1]                  + 
                    Gridz[indx+siz_line+1]           + Gridz[indx+siz_line]           + 
                    Gridz[indx+siz_slice]            + Gridz[indx+siz_slice+1]        + 
                    Gridz[indx+siz_slice+siz_line+1] + Gridz[indx+siz_line+siz_slice]  )/8.0;
                
            }
        }
    }

}


/* 
 * How many different values in the vector, 
 *  NI is the upper limitation of the v[i].
 */
int NumOfValues(std::vector<int> v, int NI) 
{
    std::vector<int> tab(NI, 0);
    int num = 0;
    for (size_t i = 0; i < v.size(); i++) {
        if (tab[v[i]] == 0)
            num++;
        tab[v[i]]++;
    }
    return num;
}

// one mesh is subdivided into ng segments in one direction
Point3 *MeshSubdivide(Mesh3 M) {

    int ng = NG;
    int siz_line = ng+1;
    int siz_slice = (ng+1) * siz_line;
    int siz_volume = (ng+1) * siz_slice;

    Point3 *gd = new Point3[siz_volume];

    for (int k = 0; k <= ng; k++) {
        for (int j = 0; j <= ng; j++) {
            for (int i = 0; i <= ng; i++) {
                int indx = i + j * siz_line + k * siz_slice;
                // 1st 0/1: y, 2nd 0/1: z
                Point3 Lx00 = M.v[0] + (M.v[1]-M.v[0])/ng*i;
                Point3 Lx10 = M.v[3] + (M.v[2]-M.v[3])/ng*i;
                Point3 Lx01 = M.v[4] + (M.v[5]-M.v[4])/ng*i;
                Point3 Lx11 = M.v[7] + (M.v[6]-M.v[7])/ng*i;
                // 0: upper plane, 1: Lower plane (for z-axis)
                Point3 Lxy0 = Lx00 + (Lx10-Lx00)/ng*j;
                Point3 Lxy1 = Lx01 + (Lx11-Lx01)/ng*j;
                // interpolation in the z-axis direction.
                gd[indx] = Lxy0 + (Lxy1-Lxy0)/ng*k;           
            }
        }
    }
        
    return gd;
}

/* 
 * Reference: Bond, 1943, The Mathematics of the Physical Properties of Crystals 
 * Zhu and Dorman, 2000, Two-dimensional, three-component wave propagation in a transversely 
 *  isotropic medium with arbitrary-orientation–finite-element modeling
 * 
 * theta: dip, phi: azimuth
 */
void BondTransform(float c11, float c12, float c13, float c14, float c15, float c16, 
                   float c22, float c23, float c24, float c25, float c26, float c33,
                   float c34, float c35, float c36, float c44, float c45, float c46,
                   float c55, float c56, float c66, float theta, float phi, 
                   float &c11_tti, float &c12_tti, float &c13_tti, float &c14_tti, float &c15_tti, float &c16_tti, 
                   float &c22_tti, float &c23_tti, float &c24_tti, float &c25_tti, float &c26_tti, float &c33_tti,
                   float &c34_tti, float &c35_tti, float &c36_tti, float &c44_tti, float &c45_tti, float &c46_tti,
                   float &c55_tti, float &c56_tti, float &c66_tti) 
{
    float a11 =  cos(theta)*sin(phi), a12 = cos(phi), a13 =  sin(theta)*sin(phi);
    float a21 = -cos(theta)*cos(phi), a22 = sin(phi), a23 = -sin(theta)*cos(phi);
    float a31 = -sin(theta),          a32 = 0,        a33 =  cos(theta);   

    float R_tmp[36] = {
        a11*a11, a12*a12, a13*a13, 2*a12*a13,       2*a13*a11,       2*a11*a12,
        a21*a21, a22*a22, a23*a23, 2*a22*a23,       2*a23*a21,       2*a21*a22,
        a31*a31, a32*a32, a33*a33, 2*a32*a33,       2*a33*a31,       2*a31*a32,
        a21*a31, a22*a32, a23*a33, a22*a33+a23*a32, a21*a33+a23*a31, a22*a31+a21*a32,
        a31*a11, a32*a12, a33*a13, a12*a33+a13*a32, a13*a31+a11*a33, a11*a31+a12*a31,
        a11*a21, a12*a22, a13*a23, a12*a23+a13*a22, a13*a21+a11*a23, a11*a22+a12*a21};
    float C0_tmp[36] = {
        c11, c12, c13, c14, c15, c16,
        c12, c22, c23, c24, c25, c26,
        c13, c23, c33, c34, c35, c36,
        c14, c24, c34, c44, c45, c46,
        c15, c25, c35, c45, c55, c56,
        c16, c26, c36, c46, c56, c66};
    
    Matrix<float> R(6,6,R_tmp);
    Matrix<float> C0(6,6,C0_tmp);
    Matrix<float> C(6,6);
    // C = R * C_vti * R^T
    C = R * C0 *R.transpose();

    c11_tti = C(1,1); c12_tti = C(1,2); c13_tti = C(1,3); c14_tti = C(1,4); c15_tti = C(1,5); c16_tti = C(1,6);
    c22_tti = C(2,2); c23_tti = C(2,3); c24_tti = C(2,4); c25_tti = C(2,5); c26_tti = C(2,6);
    c33_tti = C(3,3); c34_tti = C(3,4); c35_tti = C(3,5); c36_tti = C(3,6);
    c44_tti = C(4,4); c45_tti = C(4,5); c46_tti = C(4,6);
    c55_tti = C(5,5); c56_tti = C(5,6); c66_tti = C(6,6);
}

/* 
 * Find the last index of the x-vector greater than or equal to value, (x[i] >= value)
 *  Just for the ordered array x is from largest to smallest.
 *  (used to find the nearest elevation, from top to bottom)
 *  If value > x[all], return -1.
 */
int findLastGreaterEqualIndex(
    float value, 
    std::vector<float> &x)
{

    if (value > x[0]) return -1;

    float dist, newDist;
    int indx = -1, n = x.size();
    dist = FLT_MAX;

    for (int i = 0; i < n; i++) {
        if (x[i] == -FLT_MAX) continue;

        newDist = x[i]-value;
        
        if (newDist < 0) break; 
        if (newDist >= 0 && newDist <= dist) {
            dist = newDist;
            indx = i;
        }
    }

    return indx;
}


/* 
 * Find the index of the nearest value and no less than the value (value <= x[i]), 
 *    Just for the ordered array x is from largest to smallest.
 */
int findNearestGreaterIndex(
    float value, 
    std::vector<float> &x)
{
    float dist, newDist;
    int idx;

    idx = -1;    // expolation: INF(idx = -1) or
    dist = FLT_MAX;
    for(size_t i = 0; i < x.size(); i++) {
        if (x[i] == -FLT_MAX) continue;

        /* 
         * If the value is equal to the value in x, return the first index of the value,
         * If the value is not in the array x, return the last index which x[index] > value.
         */ 
        newDist = x[i] - value;   
        if (newDist >= 0 && newDist <= dist) {
            if (newDist <= 1e-6) 
                return i;
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

/* 
 * Find the index of the nearest value (x[i] <= value), 
 *  the x vector can be arbitrary
 */
int findNearestNeighborIndex(
    float value, 
    std::vector<float> &x)
{
    float dist, newDist;
    int idx = -1, nx = x.size();
    // expolation: INF(idx = -1) or

    dist = FLT_MAX;
    for(int i = 0; i < nx; i++) {
        if (x[i] == -FLT_MAX) continue;
        newDist = value - x[i];
        if (newDist >= 0 && newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

/* 
 * BilinearInterpolate to the (xq, vq) point 
 * It's just for the problum: (x,y) is a grid and the x, y vector are increment 
 * If x, y are arbitrary, we need to re-detrimine points for interpolation
 */
float BilinearInterpolation(
    std::vector<float> &x, 
    std::vector<float> &y, 
    float *v,
    float xq,
    float yq )
{    

    float vq;
    int NX = x.size();
    int NY = y.size();

    int ix = -1; //findNearestNeighborIndex(xq, x);
    int iy = -1; //findNearestNeighborIndex(yq, y);
    ix = (xq-x[0])/(x[1]-x[0]);  // just for the given uniform interface mash.
    iy = (yq-y[0])/(y[1]-y[0]);

    if (ix >= 0 && iy >= 0) {
        int indx = iy*NX + ix;
        if (ix == NX-1 || iy == NY-1) {
            vq = v[indx];
        }
        else {
            float area = (x[ix+1] - x[ix]) * (y[iy+1] - y[iy]); 

            vq =  v[indx]      * (x[ix+1] - xq) * (y[iy+1] - yq) 
                + v[indx+1]    * (xq -   x[ix]) * (y[iy+1] - yq)
                + v[indx+NX]   * (x[ix+1] - xq) * (yq - y[iy]  )
                + v[indx+1+NX] * (xq -   x[ix]) * (yq - y[iy]  );

            vq /= area;
        }
    }
    else
        vq = FLT_MAX;

    return vq;
}


/*---------- for matrix: just used in bond transform -----------*/
// construct: init
template<typename T>
Matrix<T>::Matrix(int r, int c):row(r), col(c){   
    p = new T[row*col];
}

// construct: init
template<typename T>
Matrix<T>::Matrix(int r, int c, T* initval){
    row = r;
    col = c;
    int sz = row*col;
    p = new T[sz];
    for (int i = 0; i < sz; i++) 
        p[i] = initval[i];
}

// construct: copy
template<typename T>
Matrix<T>::Matrix(const Matrix<T> &mat){
    row = mat.row;
    col = mat.col;
    int sz = row*col;
    p = new T[sz];
    for (int i = 0; i < sz; i++) {
        p[i] = mat.p[i];
    }
}

// destruct
template<typename T>
Matrix<T>::~Matrix(){
    if (p != nullptr)
        delete[] p;
}

// get value
template<typename T>
T& Matrix<T>::operator()(int i, int j) const {
    if (i < 0 || i > row-1 || i < 0 || j > col -1)
        throw "Out of bound!";
    return p[i*col+j];
}

// operator: euqal
template<typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &mat){
    if (row != mat.row || col != mat.col) {
        fprintf(stderr,"The row and column do not match!");
        fflush(stderr);
        exit(1);
    } 
    for (int i = 0; i < row*col; i++)
        p[i] = mat.p[i];
    return *this;
}

// multiply
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &mat){
    Matrix<T> ans(row, mat.col);
    if (col != mat.row) {
        fprintf(stderr,"Can not multiply!");
        fflush(stderr);
        exit(1);
    } else {
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < mat.col; j++) {
                ans.p[i*mat.col + j] = 0;
                for (int k = 0; k < col; k++) {
                    ans.p[i*mat.col + j] += p[i*col + k]*mat.p[k*mat.col + j];
                }
            }
        }
    }
    return ans;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() {
    Matrix<T> transposed(col, row);
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            transposed(j,i) = (*this)(i,j);
        }
    }
    return transposed;
}

// print
template<typename T>
std::ostream& operator<<(std::ostream &out, const Matrix<T> &mat) {
    for (int i = 0; i < mat.row; i++) {
        for (int j = 0; j < mat.col; j++) {
            out << mat(i,j);
            if (j == mat.col-1) 
                out << std::endl;
            else
                out << " ";
        }
    }
    return out;
}
/*--------------------------------------------*/