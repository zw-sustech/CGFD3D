/***************************************************************************
 *
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    06/2021: Created by Luqian Jiang 
 *
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version. 
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the IMPLIED WARRANTY OF
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 * 
 * You should have received a copy of the GNU General Public License along 
 * with this program. If not, see <https://www.gnu.org/licenses/>.
 *
 ***************************************************************************/


#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include <string.h>
//#include <Eigen>
#include "media_utility.hpp"

// for report error: int -> string
std::map<int, std::string> create_md2str_map() 
{
    std::map<int, std::string> m;
    m[ONE_COMPONENT] = "one_componet";
    m[ACOUSTIC_ISOTROPIC] = "acoustic_isotropic";
    m[ELASTIC_ISOTROPIC] = "elastic_isotropic"; 
    m[ELASTIC_VTI_PREM] = "elastic_vti_prem";
    m[ELASTIC_VTI_THOMSEN] = "elastic_vti_thomsen"; 
    m[ELASTIC_VTI_CIJ] = "elastic_vti_cij";
    m[ELASTIC_TTI_THOMSEN] = "elastic_tti_thomsen"; 
    m[ELASTIC_TTI_BOND] = "elastic_tti_bond";
    m[ELASTIC_ANISO_CIJ] = "elastic_aniso_cij";
    return m;
}

bool isEqual(float a, float b) {
    return fabs(a-b) < FLT_EPSILON;
}

void printProgress(float slowk) {

    if (slowk > 1) slowk = 1;
    int p = slowk * 50;

    std::cout << "\33[1A";
    std::cout << "  [" + std::string(p, '=') + ">" + std::string(50-p, ' ') << "]" << std::endl;

    fflush(stdout);
}

void PrintIsPointOutOfInterfaceRange(Point3 A, 
    int ix, int iy, int iz, 
    float MINX, float MAXX, float MINY, float MAXY) 
{
    if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
        fprintf(stderr,"Error: Grid(%d, %d, %d) = (%f, %f, %f) is out of the INTERFACES MESH (x in [%f %f], y in [%f %f]), "\
                       "please check the interfaces file!\n", ix, iy, iz, A.x, A.y, A.z, MINX, MAXX, MINY, MAXY);
        fflush(stderr);
        exit(1);
    }

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

    float dist = FLT_MAX, newDist;
    int indx = -1, n = x.size();

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
 * Find the first index of the x-vector greater than or equal to value, (x[i] >= value)
 *  Just for the ordered array x is from largest to smallest.
 *  Once there is x[i] = value, return.
 *  If value > x[all], return -1.
 */
int findFirstGreaterEqualIndex(
    float value, 
    std::vector<float> &x)
{
    float dist = FLT_MAX, newDist;
    int idx = -1, n = x.size();

    for(size_t i = 0; i < n; i++) {
        if (x[i] == -FLT_MAX) continue;
         
        newDist = x[i] - value;   
        if (newDist >= 0 && newDist <= dist) {
            // If the value is equal to the value in x, return the first index of the value,
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

Point3 PlaneBilinearIntersectSegment(
  float x0, float y0, float x1, float y1,
  float z0, float z1, float z2, float z3,
  const Point3 &P0, const Point3 &P1, bool &isInter)
{
  isInter = true;
  Vector3 d = P1 - P0;
  float area = (x1-x0)*(y1-y0);
  float t = -1.0;

  if (isEqual(d.x, 0.0) && isEqual(d.y, 0.0) && isEqual(d.z, 0.0)) {
    isInter = false;
    return P0;
  } else if (isEqual(d.x, 0.0) && isEqual(d.y, 0.0) && (!isEqual(d.z, 0.0))) {
    t = P0.x*P0.y*z0 - P0.y*x1*z0 - P0.x*y1*z0 + x1*y1*z0 - 
        P0.x*P0.y*z1 + P0.y*x0*z1 + P0.x*y1*z1 - x0*y1*z1 - 
        P0.x*P0.y*z2 + P0.y*x1*z2 + P0.x*y0*z2 - x1*y0*z2 + 
        P0.x*P0.y*z3 - P0.y*x0*z3 - P0.x*y0*z3 + x0*y0*z3 -
        area * P0.z;
    t = t*1.0/(area*d.z);
  } else if (isEqual(d.x, 0.0) && (!isEqual(d.y, 0.0)) ) {
    t = P0.x*P0.y*z0 - P0.y*x1*z0 - P0.x*y1*z0 + x1*y1*z0 - 
        P0.x*P0.y*z1 + P0.y*x0*z1 + P0.x*y1*z1 - x0*y1*z1 - 
        P0.x*P0.y*z2 + P0.y*x1*z2 + P0.x*y0*z2 - x1*y0*z2 + 
        P0.x*P0.y*z3 - P0.y*x0*z3 - P0.x*y0*z3 + x0*y0*z3 - area*P0.z;

    t /= (area*d.z - d.y*P0.x*z0 + d.y*x1*z0 + d.y*P0.x*z1 - d.y*x0*z1 + 
                     d.y*P0.x*z2 - d.y*x1*z2 - d.y*P0.x*z3 + d.y*x0*z3 );
  } else if (isEqual(d.y, 0.0) && (!isEqual(d.x, 0.0)) ){  
    t = P0.x*P0.y*z0 - P0.y*x1*z0 - P0.x*y1*z0 + x1*y1*z0 - 
        P0.x*P0.y*z1 + P0.y*x0*z1 + P0.x*y1*z1 - x0*y1*z1 - 
        P0.x*P0.y*z2 + P0.y*x1*z2 + P0.x*y0*z2 - x1*y0*z2 + 
        P0.x*P0.y*z3 - P0.y*x0*z3 - P0.x*y0*z3 + x0*y0*z3 - area*P0.z;

    t /= (area*d.z - d.x*P0.y*z0 + d.x*y1*z0 + d.x*P0.y*z1 - d.x*y1*z1 
                   + d.x*P0.y*z2 - d.x*y0*z2 - d.x*P0.y*z3 + d.x*y0*z3);
  } else {
    float t0 = d.y*P0.x*z0 + d.x*P0.y*z0 - d.y*x1*z0 - d.x*y1*z0 - 
               d.y*P0.x*z1 - d.x*P0.y*z1 + d.y*x0*z1 + d.x*y1*z1 - 
               d.y*P0.x*z2 - d.x*P0.y*z2 + d.y*x1*z2 + d.x*y0*z2 +
               d.y*P0.x*z3 + d.x*P0.y*z3 - d.y*x0*z3 - d.x*y0*z3;
    t0 /= area;
    t0 -= d.z;
    float t_sqrt1 =-d.y*P0.x*z0 - d.x*P0.y*z0 + d.y*x1*z0 + d.x*y1*z0 + 
                    d.y*P0.x*z1 + d.x*P0.y*z1 - d.y*x0*z1 - d.x*y1*z1 + 
                    d.y*P0.x*z2 + d.x*P0.y*z2 - d.y*x1*z2 - d.x*y0*z2 - 
                    d.y*P0.x*z3 - d.x*P0.y*z3 + d.y*x0*z3 + d.x*y0*z3;
    t_sqrt1 /= area;
    t_sqrt1 += d.z;
    t_sqrt1 = t_sqrt1 * t_sqrt1;
    float t_sqrt2 = -d.x*d.y*z0 + d.x*d.y*z1 + d.x*d.y*z2 - d.x*d.y*z3;
    t_sqrt2 *= (4.0/area);
    float t_sqrt3 =(P0.x*P0.y*z0 + P0.y*x1*z0 + P0.x*y1*z0 - x1*y1*z0 + 
                    P0.x*P0.y*z1 - P0.y*x0*z1 - P0.x*y1*z1 + x0*y1*z1 + 
                    P0.x*P0.y*z2 - P0.y*x1*z2 - P0.x*y0*z2 + x1*y0*z2 - 
                    P0.x*P0.y*z3 + P0.y*x0*z3 + P0.x*y0*z3 - x0*y0*z3)/area;
    t_sqrt3 = P0.z - t_sqrt3;
    float t_sqrt  = t_sqrt1 - t_sqrt2 * t_sqrt3;
    if (t_sqrt < 0) {
      fprintf(stderr, "Error: t_sqrt <0 \n");
      exit(1);
    }
    t_sqrt = sqrtf(t_sqrt);
    float t_deno = (-d.x*d.y*z0 + d.x*d.y*z1 + d.x*d.y*z2 - d.x*d.y*z3) * 2.0;
    t_deno /= area;
  
    float t1 = (t0 + t_sqrt)/t_deno;
    float t2 = (t0 - t_sqrt)/t_deno;
    if (t1 <= 1.0 && t1 >= 0.0) {
      t = t1;
    } else if (t2 <= 1.0 && t2 >= 0.0) {
      t = t2;
    } else {
      isInter = false;
    }
  }
  Point3 interP = P0 + (P1-P0)*t;
  if (t <= 1.0 && t >= 0.0) {
    isInter = true;
    return interP;
  } else {
    isInter = false;
  }

  return P0;
}

/* 
 * TrilinearInterpolate to the (xq, yq, zq) point 
 * It's just for the problum: (x,y) is a grid and the x, y vector are increment 
 * Just for bin2model, 
 */
float TrilinearInterpolation(
    std::vector<float> &x, 
    std::vector<float> &y, 
    std::vector<float> &z, 
    float *v,
    float xq,
    float yq, 
    float zq)
{    
    float vq;
    int nx = x.size();
    int ny = y.size();
    int nz = z.size();
    float dx = x[1]-x[0];
    float dy = y[1]-y[0];
    float dz = z[1]-z[0];

//    int ix0 = xq <= x[0]? 0:floor((xq-x[0])/dx); 
//    int iy0 = yq <= y[0]? 0:floor((yq-y[0])/dy);
//    int ix2 = xq >= x[nx-1]? nx-1:ceil((xq-x[0])/dx); 
//    int iy2 = yq >= y[ny-1]? ny-1:ceil((yq-y[0])/dy);
    int ix0 = floor((xq-x[0])/dx); 
    int iy0 = floor((yq-y[0])/dy);
    int ix2 = ceil((xq-x[0])/dx); 
    int iy2 = ceil((yq-y[0])/dy);
    // out of bound, or every bound
    ix0 = ix0 < 0 ? 0:ix0;
    ix2 = ix2 < 0 ? 0:ix2;
    ix0 = ix0 > nx-1 ? nx-1:ix0;
    ix2 = ix2 > nx-1 ? nx-1:ix2;
    iy0 = iy0 < 0 ? 0:iy0;
    iy2 = iy2 < 0 ? 0:iy2;
    iy0 = iy0 > ny-1 ? ny-1:iy0;
    iy2 = iy2 > ny-1 ? ny-1:iy2;

    // z can be increasing or descreasing
    int iz0 = floor((zq-z[0])/dz);
    int iz2 = ceil((zq-z[0])/dz);
    // judge out of bound by 
    iz0 = iz0 < 0 ? 0:iz0;
    iz0 = iz0 > nz-1?nz-1:iz0;
    iz2 = iz2 > nz-1?nz-1:iz2;
    iz2 = iz2 < 0 ? 0:iz2;

    size_t siz_line = nx;
    size_t siz_slice = nx*ny;

    // at the point
    if (ix0 == ix2 && iy0 == iy2 && iz0 == iz2) {
      return v[ix0 + iy0*siz_line + iz0*siz_slice];
    }

    // interp in every direction, can reduce judgement
    float v000 = v[ix0 + iy0*siz_line + iz0*siz_slice];
    float v002 = v[ix0 + iy0*siz_line + iz2*siz_slice];
    float v200 = v[ix2 + iy0*siz_line + iz0*siz_slice];
    float v202 = v[ix2 + iy0*siz_line + iz2*siz_slice];
    float v020 = v[ix0 + iy2*siz_line + iz0*siz_slice];
    float v022 = v[ix0 + iy2*siz_line + iz2*siz_slice];
    float v220 = v[ix2 + iy2*siz_line + iz0*siz_slice];
    float v222 = v[ix2 + iy2*siz_line + iz2*siz_slice];
    //- xdir
    float v100, v102, v122, v120;
    // xq = x[i]
    if (ix0 == ix2) { 
        v100 = v000;
        v102 = v002;
        v122 = v022;
        v120 = v020;
    } else {
        float w = (xq-x[ix0])/(x[ix2]-x[ix0]);
        v100 = v000 + w * (v200-v000);
        v102 = v002 + w * (v202-v002);
        v122 = v022 + w * (v222-v022);
        v120 = v020 + w * (v220-v020);   
    }
    // y-dir
    float v110, v112;
    if (iy0 == iy2) {
        v110 = v100;
        v112 = v102;
    } else {
        float w = (yq - y[iy0])/(y[iy2]-y[iy0]);
        v110 = v100 + w * (v120-v100);
        v112 = v102 + w * (v122-v102); 
    }
      
    // z-dir
    float v111;
    if (iz0 == iz2) {
        v111 = v110;
    } else {
        v111 = v110 + (zq-z[iz0])/(z[iz2]-z[iz0]) * (v112-v110); 
    }

    return v111;
}

//============== for matrix: just used in bond transform ===============
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
    if (i < 0 || i > row-1 || i < 0 || j > col-1) {
        fprintf(stderr,"i=%d, j=%d out of bound!\n", i, j);
        fflush(stderr);
        exit(1);
    }
    return p[i*col+j];
}

// operator: euqal
template<typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T> &mat){
    if (row != mat.row || col != mat.col) {
        fprintf(stderr,"The row and column do not match!\n");
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

// divided
template<typename T>
Matrix<T> Matrix<T>::operator/(T f){
  Matrix<T> ans(row, col);
  
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      ans.p[i*col + j] = p[i*col+j]/f;
    }
  }
  
  return ans;
}

// addition
template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &mat){
    Matrix<T> ans(row, col);
    if (col != mat.col || row != mat.row) {
        fprintf(stderr,"ERROR: Matrix + Two matric should be the same size!");
        fflush(stderr);
        exit(1);
    } else {
      for (int i = 0; i < row; i++) {
        for (int j = 0; j < mat.col; j++) {
          ans.p[i*col + j] = p[i*col + j] + mat.p[i*col + j];
        }
      }
    }
    return ans;
}

// minus
template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &mat){
    Matrix<T> ans(row, col);
    if (col != mat.col || row != mat.row) {
        fprintf(stderr,"ERROR: Matrix - Two matric should be the same size!");
        fflush(stderr);
        exit(1);
    } else {
      for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
          ans.p[i*col + j] = p[i*col + j] - mat.p[i*col + j];
        }
      }
    }
    return ans;
}

// transpose
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


// inverse of the 3x3 matrix
template<typename T>
Matrix<T> Matrix<T>::inverse3x3() {
  if (row != 3 || col != 3) {
      fprintf(stderr,"The row and column can just be 3!\n");
      fflush(stderr);
      exit(1);
  } 
  Matrix<T> inversed(3, 3, p);
  float inv[3][3];
  inv[0][0] = p[1*col+1]*p[2*col+2] - p[2*col+1]*p[1*col+2];
  inv[0][1] = p[2*col+1]*p[0*col+2] - p[0*col+1]*p[2*col+2];
  inv[0][2] = p[0*col+1]*p[1*col+2] - p[0*col+2]*p[1*col+1];
  inv[1][0] = p[1*col+2]*p[2*col+0] - p[1*col+0]*p[2*col+2];
  inv[1][1] = p[0*col+0]*p[2*col+2] - p[2*col+0]*p[0*col+2];
  inv[1][2] = p[1*col+0]*p[0*col+2] - p[0*col+0]*p[1*col+2];
  inv[2][0] = p[1*col+0]*p[2*col+1] - p[1*col+1]*p[2*col+0];
  inv[2][1] = p[2*col+0]*p[0*col+1] - p[0*col+0]*p[2*col+1];
  inv[2][2] = p[0*col+0]*p[1*col+1] - p[0*col+1]*p[1*col+0];

  float det = inv[0][0] * p[0*col] 
            + inv[0][1] * p[1*col] 
            + inv[0][2] * p[2*col];
  det = 1.0/det;          
        
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      inversed.p[i*col+j] = inv[i][j] * det;
    }
  }

  return inversed;
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

// re-statement for the compiler determines the type
template class Matrix<float>;
template class Matrix<int>;


//====================================================================


//========== for equivalent medium parameterization method ===========

/* 
 * How many different values in the vector, 
 *  NI is the upper limitation of the v[i].
 */
int NumOfValues(std::vector<int> v, int NI) 
{
    std::vector<int> tab(NI, 0);
    int num = 0;
    for (size_t i = 0; i < v.size(); i++) {
        // If higher than interface, do not apply equivalent medium
        if (v[i] == -1)
            return 0;

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
                // 1st 0/1: x-axis, 2nd 0/1: z-axis
                Point3 Ly00 = M.v[0] + (M.v[1]-M.v[0])/ng*j;
                Point3 Ly10 = M.v[3] + (M.v[2]-M.v[3])/ng*j;
                Point3 Ly01 = M.v[4] + (M.v[5]-M.v[4])/ng*j;
                Point3 Ly11 = M.v[7] + (M.v[6]-M.v[7])/ng*j;
                // 0: upper plane, 1: Lower plane (for z-axis)
                Point3 Lxy0 = Ly00 + (Ly10-Ly00)/ng*i;
                Point3 Lxy1 = Ly01 + (Ly11-Ly01)/ng*i;
                // interpolation in the z-axis direction.
                gd[indx] = Lxy0 + (Lxy1-Lxy0)/ng*k;           
            }
        }
    }    
    return gd;
}

void GenerateHalfGrid(
    size_t nx, 
    size_t ny, 
    size_t nz,
    int grid_type, 
    const float *Gx,   // gridx
    const float *Gy,   // gridy
    const float *Gz,   // gridz
    float **Hx,        // half gridx    
    float **Hy,        // half gridy
    float **Hz)        // half gridz
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    /* 
     * half grid point is inside the grid, 
     * so the number of grid point in one direction is np-1.
     */
    if (grid_type == GRID_CART) {
        
        *Hx = new float [nx]; 
        *Hy = new float [ny]; 
        *Hz = new float [nz]; 

        for (size_t i = 0; i < nx-1; i++)
            (*Hx)[i] = (Gx[i+1] + Gx[i])/2.0;
        (*Hx)[nx-1] = Gx[nx-1];

        for (size_t i = 0; i < ny-1; i++)
            (*Hy)[i] = (Gy[i+1] + Gy[i])/2.0;
        (*Hy)[ny-1] = Gy[ny-1];

        for (size_t i = 0; i < nz-1; i++)
            (*Hz)[i] = (Gz[i+1] + Gz[i])/2.0;
        (*Hz)[nz-1] = Gz[nz-1];

    } else if (grid_type == GRID_VMAP) {

        *Hx = new float [nx]; 
        *Hy = new float [ny]; 
        *Hz = new float [siz_volume]; 

        for (size_t i = 0; i < nx-1; i++)
            (*Hx)[i] = (Gx[i+1] + Gx[i])/2.0;
        (*Hx)[nx-1] = Gx[nx-1];

        for (size_t i = 0; i < ny-1; i++)
            (*Hy)[i] = (Gy[i+1] + Gy[i])/2.0;
        (*Hy)[ny-1] = Gy[ny-1];

        for (size_t i = 0; i < siz_volume; ++i)
            (*Hz)[i] = Gz[i];

        for (size_t k = 0; k < nz-1; k++) {
            for (size_t j = 0; j < ny-1; j++) {
                for (size_t i = 0; i < nx-1; i++) {
                    size_t indx =  i + j * siz_line + k * siz_slice;
                    (*Hz)[indx] = (Gz[indx]           + Gz[indx+1]                  + 
                        Gz[indx+siz_line+1]           + Gz[indx+siz_line]           + 
                        Gz[indx+siz_slice]            + Gz[indx+siz_slice+1]        + 
                        Gz[indx+siz_slice+siz_line+1] + Gz[indx+siz_line+siz_slice]  )/8.0;
     
                }
            }
        }
    } else {

        *Hx = new float [siz_volume]; 
        *Hy = new float [siz_volume]; 
        *Hz = new float [siz_volume]; 

        for (size_t i = 0; i < siz_volume; i++) {
            (*Hx)[i] = Gx[i];
            (*Hy)[i] = Gy[i];
            (*Hz)[i] = Gz[i];
        }

        for (size_t k = 0; k < nz-1; k++) {
            for (size_t j = 0; j < ny-1; j++) {
                for (size_t i = 0; i < nx-1; i++) {
                    size_t indx =  i + j * siz_line + k * siz_slice;
    
                    /* clockwise: conducive to debug */
                    (*Hx)[indx] = (Gx[indx]           + Gx[indx+1]                  + 
                        Gx[indx+siz_line+1]           + Gx[indx+siz_line]           + 
                        Gx[indx+siz_slice]            + Gx[indx+siz_slice+1]        + 
                        Gx[indx+siz_slice+siz_line+1] + Gx[indx+siz_line+siz_slice]  )/8.0;
    
                    (*Hy)[indx] = (Gy[indx]           + Gy[indx+1]                  + 
                        Gy[indx+siz_line+1]           + Gy[indx+siz_line]           + 
                        Gy[indx+siz_slice]            + Gy[indx+siz_slice+1]        + 
                        Gy[indx+siz_slice+siz_line+1] + Gy[indx+siz_line+siz_slice]  )/8.0;
    
                    (*Hz)[indx] = (Gz[indx]           + Gz[indx+1]                  + 
                        Gz[indx+siz_line+1]           + Gz[indx+siz_line]           + 
                        Gz[indx+siz_slice]            + Gz[indx+siz_slice+1]        + 
                        Gz[indx+siz_slice+siz_line+1] + Gz[indx+siz_line+siz_slice] )/8.0;
                    
                }
            }
        }
    }

}

Mesh3 GenerateHalfMesh(int grid_type,
                size_t ix, size_t iy, size_t iz, size_t indx, 
                size_t siz_line, size_t siz_slice, size_t siz_volume,
                float *Hx, float *Hy, float *Hz) 
{
    /* 
     * The H(ix, iy, iz) = Trilinear of G(ix, iy, iz) -> G(ix+1, iy+1, iz+1),
     *  If we want to get the averaging of the point (ix, iy, iz), 
     *  Mesh H(ix-1, iy-1, iz-1) -> H(ix, iy, iz) is need.
     */
    Mesh3 M( Point3(Hx[ix-1], Hy[iy-1], Hz[iz-1]),
             Point3(Hx[ix-1], Hy[iy  ], Hz[iz-1]),
             Point3(Hx[ix  ], Hy[iy  ], Hz[iz-1]),
             Point3(Hx[ix  ], Hy[iy-1], Hz[iz-1]),
             Point3(Hx[ix-1], Hy[iy-1], Hz[iz  ]),
             Point3(Hx[ix-1], Hy[iy  ], Hz[iz  ]),
             Point3(Hx[ix  ], Hy[iy  ], Hz[iz  ]),
             Point3(Hx[ix  ], Hy[iy-1], Hz[iz  ]) );
    
    if (grid_type == GRID_VMAP) {
        size_t indx0 = indx-1-siz_line-siz_slice;
        size_t indx1 = indx-1         -siz_slice;
        size_t indx2 = indx           -siz_slice;
        size_t indx3 = indx  -siz_line-siz_slice;
        size_t indx4 = indx-1-siz_line;
        size_t indx5 = indx-1         ;
        size_t indx6 = indx           ;
        size_t indx7 = indx  -siz_line;

        Mesh3 M( Point3(Hx[ix-1], Hy[iy-1], Hz[indx0]),
                 Point3(Hx[ix-1], Hy[iy  ], Hz[indx1]),
                 Point3(Hx[ix  ], Hy[iy  ], Hz[indx2]),
                 Point3(Hx[ix  ], Hy[iy-1], Hz[indx3]),
                 Point3(Hx[ix-1], Hy[iy-1], Hz[indx4]),
                 Point3(Hx[ix-1], Hy[iy  ], Hz[indx5]),
                 Point3(Hx[ix  ], Hy[iy  ], Hz[indx6]),
                 Point3(Hx[ix  ], Hy[iy-1], Hz[indx7]) );
        return M;

    } else if (grid_type == GRID_CURV) {
        size_t indx0 = indx-1-siz_line-siz_slice;
        size_t indx1 = indx-1         -siz_slice;
        size_t indx2 = indx           -siz_slice;
        size_t indx3 = indx  -siz_line-siz_slice;
        size_t indx4 = indx-1-siz_line;
        size_t indx5 = indx-1         ;
        size_t indx6 = indx           ;
        size_t indx7 = indx  -siz_line;
        Mesh3 M( Point3( Hx[indx0], Hy[indx0], Hz[indx0] ),
                 Point3( Hx[indx1], Hy[indx1], Hz[indx1] ),
                 Point3( Hx[indx2], Hy[indx2], Hz[indx2] ),
                 Point3( Hx[indx3], Hy[indx3], Hz[indx3] ),
                 Point3( Hx[indx4], Hy[indx4], Hz[indx4] ),
                 Point3( Hx[indx5], Hy[indx5], Hz[indx5] ),
                 Point3( Hx[indx6], Hy[indx6], Hz[indx6] ),
                 Point3( Hx[indx7], Hy[indx7], Hz[indx7] ));
        return M;
    }

    return M;
}

//======================for vti and tti ====================================
void para2vti(
    std::vector<float> var, // input var
    int media_type, // return cij
    float &c11_3d,
    float &c33_3d,
    float &c55_3d,
    float &c66_3d,
    float &c13_3d,
    float &rho_3d) 
{
    /* Calculate cij */
    if (media_type == ELASTIC_VTI_PREM) {
    //-r eference: Dziewonski and Anderson, 1981, Preliminary reference Earth model
    /* rho, vph, vpv, vsh, vsv, rho */
        float vph = var[1], vpv = var[2];
        float vsh = var[3], vsv = var[4];
        float rho = var[0], eta = var[5];
    
        c11_3d = vph*vph*rho;
        c33_3d = vpv*vpv*rho;
        c66_3d = vsh*vsh*rho;
        c55_3d = vsv*vsv*rho;
        rho_3d = rho; 
        c13_3d = eta*(c11_3d - 2.0*c55_3d);

    } else if (media_type == ELASTIC_VTI_THOMSEN) {
    //- reference: Thomsen, 1986, Weak elastic anisotropy 
    /* rho, vp0 (alpha), vs0 (beta), epsilon, delta, gamma */
        float rho = var[0];
        float vp0 = var[1], vs0 = var[2];
        float epsil = var[3], delta = var[4];
        float gamma = var[5];
        float c33 = vp0*vp0*rho, c55 = vs0*vs0*rho;
        rho_3d = rho; 
        c33_3d = c33;
        c55_3d = c55;                
        c66_3d = 2.0*gamma * c55 + c55;
        c11_3d = 2.0*epsil * c33 + c33;
        c13_3d = sqrt( 2.0*delta*c33*(c33-c55) + (c33-c55)*(c33-c55) ) - c55;

    } else if (media_type == ELASTIC_VTI_CIJ) {
    /* rho, c11 c33 c55 c66 c13*/
        c11_3d = var[1];
        c33_3d = var[2];
        c55_3d = var[3];
        c66_3d = var[4];
        c13_3d = var[5];
        rho_3d = var[0];

    } else {
        fprintf(stderr,"Error: Unknow VTI type, for code check, please contact Luqian Jiang!\n");
        fflush(stderr);
        exit(1);
    }
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
    theta *= (PI/180.0);
    phi   *= (PI/180.0);
    float a11 =  cos(theta)*sin(phi), a12 = cos(phi), a13 =  sin(theta)*sin(phi);
    float a21 = -cos(theta)*cos(phi), a22 = sin(phi), a23 = -sin(theta)*cos(phi);
    float a31 = -sin(theta),          a32 = 0,        a33 =  cos(theta);   

    float R_tmp[36] = {
        a11*a11, a12*a12, a13*a13,     2*a12*a13,     2*a13*a11,     2*a11*a12,
        a21*a21, a22*a22, a23*a23,     2*a22*a23,     2*a23*a21,     2*a21*a22,
        a31*a31, a32*a32, a33*a33,     2*a32*a33,     2*a33*a31,     2*a31*a32,
        a21*a31, a22*a32, a23*a33, a22*a33+a23*a32, a21*a33+a23*a31, a22*a31+a21*a32,
        a31*a11, a32*a12, a33*a13, a12*a33+a13*a32, a13*a31+a11*a33, a11*a32+a12*a31,
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

    c11_tti = C(0,0); c12_tti = C(0,1); c13_tti = C(0,2); c14_tti = C(0,3); c15_tti = C(0,4); c16_tti = C(0,5);
    c22_tti = C(1,1); c23_tti = C(1,2); c24_tti = C(1,3); c25_tti = C(1,4); c26_tti = C(1,5);
    c33_tti = C(2,2); c34_tti = C(2,3); c35_tti = C(2,4); c36_tti = C(2,5);
    c44_tti = C(3,3); c45_tti = C(3,4); c46_tti = C(3,5);
    c55_tti = C(4,4); c56_tti = C(4,5); c66_tti = C(5,5);
}

/* 
 * Bond transform according to a new coordinate system
 */
void BondTransform_CoordSystem(float c11, float c12, float c13, float c14, float c15, float c16, 
                               float c22, float c23, float c24, float c25, float c26, float c33,
                               float c34, float c35, float c36, float c44, float c45, float c46,
                               float c55, float c56, float c66, 
                               const Vector3 &a, const Vector3 &b, const Vector3 &c, 
                               float &c11_tti, float &c12_tti, float &c13_tti, float &c14_tti, float &c15_tti, float &c16_tti, 
                               float &c22_tti, float &c23_tti, float &c24_tti, float &c25_tti, float &c26_tti, float &c33_tti,
                               float &c34_tti, float &c35_tti, float &c36_tti, float &c44_tti, float &c45_tti, float &c46_tti,
                               float &c55_tti, float &c56_tti, float &c66_tti) 
{
    Vector3 xaxis(1.0, 0.0, 0.0), yaxis(0.0, 1.0, 0.0), zaxis(0.0, 0.0, 1.0);

    float xa = CosinDir(xaxis, a),  xb = CosinDir(xaxis, b), xc = CosinDir(xaxis, c);
    float ya = CosinDir(yaxis, a),  yb = CosinDir(yaxis, b), yc = CosinDir(yaxis, c);
    float za = CosinDir(zaxis, a),  zb = CosinDir(zaxis, b), zc = CosinDir(zaxis, c);

    float R_tmp[36] = {
        xa*xa, xb*xb, xc*xc,     2*xb*xc,     2*xc*xa,     2*xa*xb,
        ya*ya, yb*yb, yc*yc,     2*yb*yc,     2*yc*ya,     2*ya*yb,
        za*za, zb*zb, zc*zc,     2*zb*zc,     2*zc*za,     2*za*zb,
        ya*za, yb*zb, yc*zc, yb*zc+yc*zb, ya*zc+yc*za, yb*za+ya*zb,
        za*xa, zb*xb, zc*xc, xb*zc+xc*zb, xc*za+xa*zc, xa*zb+xb*za,
        xa*ya, xb*yb, xc*yc, xb*yc+xc*yb, xc*ya+xa*yc, xa*yb+xb*ya};

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
    C = R * C0 * R.transpose();

    c11_tti = C(0,0); c12_tti = C(0,1); c13_tti = C(0,2); c14_tti = C(0,3); c15_tti = C(0,4); c16_tti = C(0,5);
    c22_tti = C(1,1); c23_tti = C(1,2); c24_tti = C(1,3); c25_tti = C(1,4); c26_tti = C(1,5);
    c33_tti = C(2,2); c34_tti = C(2,3); c35_tti = C(2,4); c36_tti = C(2,5);
    c44_tti = C(3,3); c45_tti = C(3,4); c46_tti = C(3,5);
    c55_tti = C(4,4); c56_tti = C(4,5); c66_tti = C(5,5);
}


void para2tti(std::vector<float> var, // input var
             int media_type, // return cij
             float &c11_3d,
             float &c12_3d,
             float &c13_3d,
             float &c14_3d,
             float &c15_3d,
             float &c16_3d,
             float &c22_3d,
             float &c23_3d,
             float &c24_3d,
             float &c25_3d,
             float &c26_3d,
             float &c33_3d,
             float &c34_3d,
             float &c35_3d,
             float &c36_3d,
             float &c44_3d,
             float &c45_3d,
             float &c46_3d,
             float &c55_3d,
             float &c56_3d,
             float &c66_3d,
             float &rho_3d)

{
    if (media_type == ELASTIC_TTI_THOMSEN) {
    /* reference: Thomsen, 1986, Weak elastic anisotropy */    
        /* rho, vp0, vs0, epsilon, delta, gamma, azimuth, dip */
        /* Calculate cij */
        float rho = var[0];
        float vp0 = var[1], vs0 = var[2];
        float epsil = var[3], delta = var[4];
        float gamma = var[5];
        float azimuth = var[6];
        float dip = var[7];
        float c33 = vp0*vp0*rho;
        float c44 = vs0*vs0*rho;
        float c66 = 2.0*gamma  *c44 + c44;
        float c11 = 2.0*epsil*c33 + c33;
        float c13 = sqrt( 2.0*delta*c33*(c33-c44) + (c33-c44)*(c33-c44) ) - c44;
        rho_3d = rho; 
    
        BondTransform(c11, c11-2*c66, c13, 0, 0, 0,  // cij
                      c11, c13, 0, 0, 0,             // c2j: 23-26
                      c33, 0, 0, 0,                  // c3j
                      c44, 0, 0,                     // c4j
                      c44, 0, c66, var[7], var[6],
                      c11_3d, c12_3d, c13_3d, c14_3d, c15_3d, c16_3d,
                      c22_3d, c23_3d, c24_3d, c25_3d, c26_3d, c33_3d,
                      c34_3d, c35_3d, c36_3d, c44_3d, c45_3d, c46_3d,
                      c55_3d, c56_3d, c66_3d);

    } else if (media_type == ELASTIC_TTI_BOND) {

        /* rho c11 c33 c55 c66 c13 azimuth dip */
        float c11 = var[1];
        float c33 = var[2];
        float c44 = var[3];
        float c66 = var[4];
        float c13 = var[5];
        rho_3d = var[0];
    
        BondTransform(c11, c11-2*c66, c13, 0, 0, 0,  // cij
                      c11, c13, 0, 0, 0,             // c2j: 23-26
                      c33, 0, 0, 0,                  // c3j
                      c44, 0, 0,                     // c4j
                      c44, 0, c66, var[7], var[6],
                      c11_3d, c12_3d, c13_3d, c14_3d, c15_3d, c16_3d,
                      c22_3d, c23_3d, c24_3d, c25_3d, c26_3d, c33_3d,
                      c34_3d, c35_3d, c36_3d, c44_3d, c45_3d, c46_3d,
                      c55_3d, c56_3d, c66_3d);

    } else if (media_type == ELASTIC_ANISO_CIJ) {
    /* rho c11 c12 c13 c14 c15 c16 c22 c23 c24 c25 c26 c33 c34 c35 c36 c44 c45 c46 c55 c56 c66*/
    
        /* Calculate cij */
        rho_3d = var[0];
        c11_3d = var[1];
        c12_3d = var[2];
        c13_3d = var[3];
        c14_3d = var[4];
        c15_3d = var[5];
        c16_3d = var[6];
        c22_3d = var[7];
        c23_3d = var[8];
        c24_3d = var[9];
        c25_3d = var[10];
        c26_3d = var[11];
        c33_3d = var[12];
        c34_3d = var[13];
        c35_3d = var[14];
        c36_3d = var[15];
        c44_3d = var[16];
        c45_3d = var[17];
        c46_3d = var[18];
        c55_3d = var[19];
        c56_3d = var[20];
        c66_3d = var[21];

    } else {
        fprintf(stderr,"Error: Unknow TTI type, for code check, please contact Luqian Jiang!\n");
        fflush(stderr);
        exit(1);
    }
}


//================== edge table
// edge table for 
int edgeTable[256]={
    0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
    0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
    0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
    0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
    0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
    0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
    0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
    0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
    0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
    0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
    0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
    0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
    0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
    0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
    0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
    0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
    0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
    0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
    0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
    0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
    0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
    0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
    0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
    0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
    0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
    0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
    0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
    0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
    0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
    0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
    0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0   };

//============================= fitting plane ======================================

//bool FitHeightPlane(std::set<Point3> &Ps, Plane &plane)
//{
//  int n = Ps.size();
//  if (n < 3) return false;
//
//  float sum_xx = 0.0, sum_xy = 0.0, sum_x = 0.0;
//  float sum_yy = 0.0, sum_y = 0.0;
//  float sum_xz = 0.0, sum_yz = 0.0, sum_z = 0.0;
//  for (std::set<Point3>::iterator it = Ps.begin(); it != Ps.end(); it++){
//    sum_xx += (it->x * it->x);
//    sum_xy += (it->x * it->y);
//    sum_yy += (it->y * it->y);
//    sum_xz += (it->x * it->z);
//    sum_yz += (it->y * it->z);
//    sum_x  += (it->x);
//    sum_y  += (it->y);
//    sum_z  += (it->z);
//  } 
//  float a[9] = {sum_xx, sum_xy, sum_x, 
//                sum_xy, sum_yy, sum_y, 
//                sum_x, sum_y, (float)n};
//  float b[3] = {sum_xz, sum_yz, sum_z};
//  Matrix<float> A(3,3, a), B(3,1, b);
//
//  Matrix<float> x = A.inverse3x3() * B;
//
//  plane.a = x(0,0);
//  plane.b = x(1,0);
//  plane.c =  -1.0;
//  plane.d = x(2,0);
//
//  return true;
//}
//

bool FitHeightPlane(std::set<Point3> &Ps, Plane &plane)
{
  int n = Ps.size();
  if (n < 3) return false;

  Vector3 mean(0,0,0);
  for (std::set<Point3>::iterator it = Ps.begin(); it != Ps.end(); it++){
    mean = mean + *it;
  }
  mean = mean/(n*1.0);

  float sum_xx = 0.0, sum_xy = 0.0, sum_yy = 0.0;
  float sum_xz = 0.0, sum_yz = 0.0;
  for (std::set<Point3>::iterator it = Ps.begin(); it != Ps.end(); it++){
    Vector3 diff = *it - mean;

    sum_xx += (diff.x * diff.x);
    sum_xy += (diff.x * diff.y);
    sum_yy += (diff.y * diff.y);
    sum_xz += (diff.x * diff.z);
    sum_yz += (diff.y * diff.z);
  }
  float det = sum_xx*sum_yy - sum_xy*sum_xy;
  if (!isEqual(det, 0.0)) {
    float a0 = (sum_yy*sum_xz-sum_xy*sum_yz )/det;
    float a1 = (sum_xx*sum_yz-sum_xy*sum_xz )/det;
    float a2 = mean.z - a0*mean.x -a1*mean.y;
    plane.a = -1.0*a0;
    plane.b = -1.0*a1;
    plane.c = 1.0;
    plane.d = -1.0*a2;
    return true;
  } else {
    plane = Plane(0,0,0,0);
    return false;
  }

  return true;
}

// build new coordinate from plane
bool buildCoord_fromPlane(Plane &plane, Vector3 &u, Vector3 &v, Vector3 &n) 
{
  if (isEqual(plane.a,0.0) && isEqual(plane.b, 0.0) && isEqual(plane.c, 0.0))
    return false;

  n = {plane.a, plane.b, plane.c};
  n = n/sqrtf(n.x*n.x + n.z*n.z + n.y*n.y);

  if (fabs(n.x) >= fabs(n.y) ) {
    float invLength = 1.0/sqrtf(n.x*n.x + n.z*n.z);
    u.x = n.z*invLength;
    u.y = 0;
    u.z = -n.x*invLength;
  } else {
    float invLength = 1.0/sqrt(n.y*n.y + n.z*n.z); 
    u.x = 0.0;
    u.y =  n.z * invLength;
    u.z = -n.y * invLength;
  }
  v = Cross(n,u);
  return true;
}
