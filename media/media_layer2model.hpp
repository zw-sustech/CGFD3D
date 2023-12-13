/***********************************************************************
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
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


#ifndef _MEDIA_LAYER2MODEL_
#define _MEDIA_LAYER2MODEL_

#include "media_geometry3d.hpp"
#include "media_utility.hpp"


#ifdef __cplusplus
extern "C" {
#endif
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

#ifdef __cplusplus
}
#endif



int AssignLayerMediaPara2Point(
    size_t ix, size_t iy, size_t iz,         /* To print error messages */ 
    Point3 A,  
    inter_t &interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var); 

//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for layer2model
void CalPointValue_layer(int media_type, 
                   inter_t &interfaces,
                   size_t slice, 
                   std::vector<float> &xvec,  /* interface mesh */
                   std::vector<float> &yvec,
                   Point3 &A,
                   std::vector<float> &elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var);

//- 0. assign the parameter directly (use the local values): one component 
void parametrization_layer_onecmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *var3d,
    int myid);

//- 1. assign the parameter directly (use the local values): isotropic, acoustic 
void parametrization_layer_ac_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *kappa,
    float *rho3d,
    int myid);

//- 2. assign the parameter directly (use the local values): elastic isotropic 
void parametrization_layer_el_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d,
    int myid);

//- 3. Assign the parameter directly (use the local values): elastic vti
void parametrization_layer_el_vti_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid);

//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_layer_el_aniso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
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
    float *rho,
    int myid);

void MarkInterfaceNumber(
        int grid_type,
        float *Hx, float *Hy, float *Hz,
        size_t nx, size_t ny, size_t nz,
        int *MaterNum, // nx*ny*nz
        inter_t &interfaces);

//- 0.1 Assign the parameter by volume harmonic averaging
//- one component
void parametrization_layer_onecmp_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *var3d,
    int myid);

//- 0.1 Assign the parameter by volume arithmetic averaging
//- one component
void parametrization_layer_onecmp_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *var3d,
    int myid);

//- 1.0 Assign the parameter by volume harmonic averaging (kappa)
//- acoustic isortopic
void parametrization_layer_ac_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *kappa, 
    float *rho3d,
    int myid);

//- 1.1 Assign the parameter by volume arithmetic averaging (kappa)
//- acoustic isortopic
void parametrization_layer_ac_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *kappa, 
    float *rho3d,
    int myid);

//- 2.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic isotropic
// ref: Moczo, 2002, 3D Heterogeneous Staggered-Grid Finite-Difference Modeling of Seismic
//      Motion with Volume Harmonic and Arithmetic Averaging of Elastic Moduli and Densities
void parametrization_layer_el_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d,
    int myid);

//- 2.1 Assign the parameter by volume arithmetic averaging method 
//- elasic isotropic
void parametrization_layer_el_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    inter_t &interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d,
    int myid);

//- 3.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic vti
void parametrization_layer_el_vti_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid);

//- 3.1 Assign the parameter by volume arithmetic averaging method 
//- elasic vti
void parametrization_layer_el_vti_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid);

//- 4.0 Assign the parameter by volume arithmetic and harmonic averaging method 
//- elasic aniso
void parametrization_layer_el_aniso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
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
    float *rho,
    int myid);

//- 4.1 Assign the parameter by volume arithmetic averaging method 
//- elasic tti
void parametrization_layer_el_aniso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
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
    float *rho,
    int myid);


void parametrization_layer_el_iso_tti(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    inter_t &interfaces,
    float *rho,
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
    int myid);

bool getInterfaceIntersection(
  int nx, int ny, float dx, float dy, 
  float x0, float y0, const float *elevation,
  const Point3 &v1, const Point3 &v2,
  Point3 &interP);

void getInersectList(int nx, int ny, float dx, float dy, float x0, float y0, 
                const float *elev, Mesh3 M, int edgeState, 
                std::set<Point3> &intersectionList);

void cal_iso_tti_parametrization_val(
    int i, int j, int k,
    inter_t &interfaces, Point3 *SubGrid, 
    const Vector3 &a, const Vector3 &b, const Vector3 &c, 
    float &c11, float &c12, float &c13, float &c14, float &c15, float &c16, 
    float &c22, float &c23, float &c24, float &c25, float &c26,
    float &c33, float &c34, float &c35, float &c36,
    float &c44, float &c45, float &c46,
    float &c55, float &c56, float &c66, float &rho3d);

void cal_ort_represent_val(
    int indx_i, int indx_j, int indx_k,
    inter_t &interfaces, Point3 *SubGrid,
    float &c11, float &c12, float &c13, 
    float &c22, float &c23, float &c33,
    float &c44, float &c55, float &c66, float &rho3d);


#endif /* __MEDID_LAYER2MODEL__ */
