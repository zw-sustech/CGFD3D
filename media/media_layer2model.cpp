/***************************************************************************
 *
 * This function is used for medium parameterization.
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    06/2021: Created by Luqian Jiang 
 *
 ***************************************************************************/
#include <iostream>
#include <vector>
#include <string.h>
#include <cmath>
#include "media_geometry3d.hpp"
#include "media_layer2model.hpp"
#include "media_read_file.hpp"
#include "media_utility.hpp"

//using namespace std;

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

int AssignLayerMediaPara2Point(
    int ix, int iy, int iz,         /* To print error messages */ 
    Point3 A,  
    inter_t interfaces,
    int media_type,                /* the type can be found in media_utility.hpp */ 
    std::vector<float> &var)
{
    size_t  NI = interfaces.NI;
    size_t  NX = interfaces.NX;
    size_t  NY = interfaces.NY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;

    size_t interface_slice = NX*NY;

    /* If out of the INTERFACE MESH area, exit! */
    PrintIsPointOutOfInterfaceRange(A, ix, iy, iz, MINX, MAXX, MINY, MAXY);

    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> elevation(NI, -FLT_MAX);
    for (size_t i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (size_t i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }

    /* 
     * For each interface, interpolate to get the elevation of 
     *  point A at the projection position of the interface mesh.
     */
    for (int ni = 0; ni < NI; ni++) {
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*interface_slice, A.x, A.y);
    }

    /* Find which material_index to use */
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    if (mi == -1) {
        fprintf(stderr,"Warning: z-location of Grid(%d, %d, %d) = (%f, %f, %f) is higher than the given elevation in " \
            "the interfaces file, it assigned by the top medium! \n", ix, iy, iz, A.x, A.y, A.z);
        fflush(stderr);
    }
    
    CalPointValue(media_type, interfaces, interface_slice, XVEC, YVEC, A, elevation, mi, var);

    return 0;
}

/* Calculate the value of the point for different media type (to avoid multiple geometric calculations) */
void CalPointValue(int media_type, 
                   inter_t interfaces,
                   size_t slice, 
                   std::vector<float> xvec,  /* interface mesh */
                   std::vector<float> yvec,
                   Point3 A,
                   std::vector<float> elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var)
{
    float dz = elevation[mi] - A.z;

    switch(media_type)
    {
    case ONE_COMPONENT: /*0. var*/
        /* If grid_z > elevation of top_interface, it given by the medium of top non-zero thickness layer */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.var + mi*slice, A.x, A.y); 
        } else {
            float var0      = BilinearInterpolation(xvec, yvec, interfaces.var      + mi*slice , A.x, A.y);
            float var_grad  = BilinearInterpolation(xvec, yvec, interfaces.var_grad + mi*slice , A.x, A.y);
            float var_pow   = BilinearInterpolation(xvec, yvec, interfaces.var_pow  + mi*slice , A.x, A.y);
            var[0]  = var0  + pow(dz, var_pow)* var_grad;
        }
    break;

    case ELASTIC_ISOTROPIC: /*1. vp, vs, rho*/
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vs  + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice , A.x, A.y);
        } else {
            float vp_grad  = BilinearInterpolation(xvec, yvec, interfaces.vp_grad + mi*slice, A.x, A.y);
            float vp0      = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
            float vp_pow   = BilinearInterpolation(xvec, yvec, interfaces.vp_pow  + mi*slice, A.x, A.y);
            float vs0      = BilinearInterpolation(xvec, yvec, interfaces.vs      + mi*slice, A.x, A.y);
            float vs_grad  = BilinearInterpolation(xvec, yvec, interfaces.vs_grad + mi*slice, A.x, A.y);
            float vs_pow   = BilinearInterpolation(xvec, yvec, interfaces.vs_pow  + mi*slice, A.x, A.y);
            float rho0     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad+ mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice, A.x, A.y);
            var[0] = vp0  + pow(dz, vp_pow)* vp_grad;
            var[1] = vs0  + pow(dz, vs_pow)* vs_grad;
            var[2] = rho0 + pow(dz,rho_pow)*rho_grad;
        }
    break;

    case ELASTIC_VTI_PREM: /*2. vph, vpv, vsh, vsv, rho, eta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
        } else {
            float vph = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            float vpv = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            float vsh = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            float vsv = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float eta = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
            float vph_grad = BilinearInterpolation(xvec, yvec, interfaces.vph_grad + mi*slice, A.x, A.y);
            float vpv_grad = BilinearInterpolation(xvec, yvec, interfaces.vpv_grad + mi*slice, A.x, A.y);
            float vsh_grad = BilinearInterpolation(xvec, yvec, interfaces.vsh_grad + mi*slice, A.x, A.y);
            float vsv_grad = BilinearInterpolation(xvec, yvec, interfaces.vsv_grad + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
            float eta_grad = BilinearInterpolation(xvec, yvec, interfaces.eta_grad + mi*slice, A.x, A.y);
            float vph_pow  = BilinearInterpolation(xvec, yvec, interfaces.vph_pow  + mi*slice, A.x, A.y);
            float vpv_pow  = BilinearInterpolation(xvec, yvec, interfaces.vpv_pow  + mi*slice, A.x, A.y);
            float vsh_pow  = BilinearInterpolation(xvec, yvec, interfaces.vsh_pow  + mi*slice, A.x, A.y);
            float vsv_pow  = BilinearInterpolation(xvec, yvec, interfaces.vsv_pow  + mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
            float eta_pow  = BilinearInterpolation(xvec, yvec, interfaces.eta_pow  + mi*slice, A.x, A.y);
            var[0] = vph + pow(dz, vph_pow)* vph_grad;
            var[1] = vpv + pow(dz, vpv_pow)* vpv_grad;
            var[2] = vsh + pow(dz, vsh_pow)* vsh_grad;
            var[3] = vsv + pow(dz, vsv_pow)* vsv_grad;
            var[4] = rho + pow(dz, rho_pow)* rho_grad;
            var[5] = eta + pow(dz, eta_pow)* eta_grad;
        }   
    break;

    case ELASTIC_VTI_THOMSEN: /*3. vpv, vsv, epsilon, delta, gamma, rho */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.vpv    + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vsv    + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.epsilon+ mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.delta  + mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.gamma  + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.rho    + mi*slice , A.x, A.y);
        } else {
            float vpv     = BilinearInterpolation(xvec, yvec, interfaces.vpv     + mi*slice, A.x, A.y);
            float vsv     = BilinearInterpolation(xvec, yvec, interfaces.vsv     + mi*slice, A.x, A.y);
            float epsilon = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);
            float delta   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            float rho     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float vpv_grad     = BilinearInterpolation(xvec, yvec, interfaces.vpv_grad     + mi*slice, A.x, A.y);
            float vsv_grad     = BilinearInterpolation(xvec, yvec, interfaces.vsv_grad     + mi*slice, A.x, A.y);
            float epsilon_grad = BilinearInterpolation(xvec, yvec, interfaces.epsilon_grad + mi*slice, A.x, A.y);
            float delta_grad   = BilinearInterpolation(xvec, yvec, interfaces.delta_grad   + mi*slice, A.x, A.y);
            float gamma_grad   = BilinearInterpolation(xvec, yvec, interfaces.gamma_grad   + mi*slice, A.x, A.y);
            float rho_grad     = BilinearInterpolation(xvec, yvec, interfaces.rho_grad     + mi*slice, A.x, A.y);
            float vpv_pow      = BilinearInterpolation(xvec, yvec, interfaces.vpv_pow      + mi*slice, A.x, A.y);
            float vsv_pow      = BilinearInterpolation(xvec, yvec, interfaces.vsv_pow      + mi*slice, A.x, A.y);
            float epsilon_pow  = BilinearInterpolation(xvec, yvec, interfaces.epsilon_pow  + mi*slice, A.x, A.y);
            float delta_pow    = BilinearInterpolation(xvec, yvec, interfaces.delta_pow    + mi*slice, A.x, A.y);
            float gamma_pow    = BilinearInterpolation(xvec, yvec, interfaces.gamma_pow    + mi*slice, A.x, A.y);
            float rho_pow      = BilinearInterpolation(xvec, yvec, interfaces.rho_pow      + mi*slice, A.x, A.y);
            var[0] = vpv     + pow(dz, vpv_pow)    * vpv_grad;
            var[1] = vsv     + pow(dz, vsv_pow)    * vsv_grad;
            var[2] = epsilon + pow(dz, epsilon_pow)* epsilon_grad;
            var[3] = delta   + pow(dz, delta_pow)  * delta_grad;
            var[4] = gamma   + pow(dz, gamma_pow)  * gamma_grad;
            var[5] = rho     + pow(dz, rho_pow)    * rho_grad;
        }   
    break;

    case ELASTIC_VTI_CIJ: /*4. c11 c33 c44 c66 c13 rho*/
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
        } else {
            float c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c44 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            float c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice, A.x, A.y);
            float c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice, A.x, A.y);
            float c44_grad = BilinearInterpolation(xvec, yvec, interfaces.c44_grad + mi*slice, A.x, A.y);
            float c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice, A.x, A.y);
            float c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
            float c11_pow  = BilinearInterpolation(xvec, yvec, interfaces.c11_pow  + mi*slice, A.x, A.y);
            float c33_pow  = BilinearInterpolation(xvec, yvec, interfaces.c33_pow  + mi*slice, A.x, A.y);
            float c44_pow  = BilinearInterpolation(xvec, yvec, interfaces.c44_pow  + mi*slice, A.x, A.y);
            float c66_pow  = BilinearInterpolation(xvec, yvec, interfaces.c66_pow  + mi*slice, A.x, A.y);
            float c13_pow  = BilinearInterpolation(xvec, yvec, interfaces.c13_pow  + mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
            var[0] = c11 + pow(dz, c11_pow)* c11_grad;
            var[1] = c33 + pow(dz, c33_pow)* c33_grad;
            var[2] = c44 + pow(dz, c44_pow)* c44_grad;
            var[3] = c66 + pow(dz, c66_pow)* c66_grad;
            var[4] = c13 + pow(dz, c13_pow)* c13_grad;
            var[5] = rho + pow(dz, rho_pow)* rho_grad;
        }   
    break;

    case ELASTIC_TTI_THOMSEN: /*5. vp0, vs0, epsilon, delta, gamma, rho, azimuth, dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vs      + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice , A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice , A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice , A.x, A.y);
        } else {
            float vp0     = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
            float vs0     = BilinearInterpolation(xvec, yvec, interfaces.vs      + mi*slice, A.x, A.y);
            float epsilon = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);
            float delta   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            float rho     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float azimuth = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            float dip     = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
            float vp0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vp_grad      + mi*slice, A.x, A.y);
            float vs0_grad     = BilinearInterpolation(xvec, yvec, interfaces.vs_grad      + mi*slice, A.x, A.y);
            float epsilon_grad = BilinearInterpolation(xvec, yvec, interfaces.epsilon_grad + mi*slice, A.x, A.y);
            float delta_grad   = BilinearInterpolation(xvec, yvec, interfaces.delta_grad   + mi*slice, A.x, A.y);
            float gamma_grad   = BilinearInterpolation(xvec, yvec, interfaces.gamma_grad   + mi*slice, A.x, A.y);
            float rho_grad     = BilinearInterpolation(xvec, yvec, interfaces.rho_grad    + mi*slice, A.x, A.y);
            float azimuth_grad = BilinearInterpolation(xvec, yvec, interfaces.azimuth_grad + mi*slice, A.x, A.y);
            float dip_grad     = BilinearInterpolation(xvec, yvec, interfaces.dip_grad     + mi*slice, A.x, A.y);
            float vp0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vp_pow       + mi*slice, A.x, A.y);
            float vs0_pow      = BilinearInterpolation(xvec, yvec, interfaces.vs_pow       + mi*slice, A.x, A.y);
            float epsilon_pow  = BilinearInterpolation(xvec, yvec, interfaces.epsilon_pow  + mi*slice, A.x, A.y);
            float delta_pow    = BilinearInterpolation(xvec, yvec, interfaces.delta_pow    + mi*slice, A.x, A.y);
            float gamma_pow    = BilinearInterpolation(xvec, yvec, interfaces.gamma_pow    + mi*slice, A.x, A.y);
            float rho_pow      = BilinearInterpolation(xvec, yvec, interfaces.rho_pow     + mi*slice, A.x, A.y);
            float azimuth_pow  = BilinearInterpolation(xvec, yvec, interfaces.azimuth_pow  + mi*slice, A.x, A.y);
            float dip_pow      = BilinearInterpolation(xvec, yvec, interfaces.dip_pow      + mi*slice, A.x, A.y);
            var[0] = vp0     + pow(dz, vp0_pow)    * vp0_grad;
            var[1] = vs0     + pow(dz, vs0_pow)    * vs0_grad;
            var[2] = epsilon + pow(dz, epsilon_pow)* epsilon_grad;
            var[3] = delta   + pow(dz, delta_pow)  * delta_grad;
            var[4] = gamma   + pow(dz, gamma_pow)  * gamma_grad;
            var[5] = rho     + pow(dz, rho_pow)    * rho_grad;
            var[6] = azimuth + pow(dz, azimuth_pow)* azimuth_grad;
            var[7] = dip     + pow(dz, dip_pow)    * dip_grad;
        }   
    break;

    case ELASTIC_TTI_CIJ: /* 7. c11 c12 c13 c14 c15 c16 c22 c23 c24 c25 c26 c33 c34 c35 c36 c44 c45 c46 c55 c56 c66 rho*/
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0]  = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[1]  = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            var[2]  = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[3]  = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            var[4]  = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            var[5]  = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            var[6]  = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            var[7]  = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            var[8]  = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            var[9]  = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            var[10] = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            var[11] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[12] = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            var[13] = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            var[14] = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            var[15] = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            var[16] = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            var[17] = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            var[18] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[19] = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            var[20] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[21] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
        } else {
            float c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c12 = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            float c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float c14 = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            float c15 = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            float c16 = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            float c22 = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            float c23 = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            float c24 = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            float c25 = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            float c26 = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            float c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c34 = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            float c35 = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            float c36 = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            float c44 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            float c45 = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            float c46 = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            float c55 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c56 = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            float c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_pow = BilinearInterpolation(xvec, yvec, interfaces.c11_pow + mi*slice , A.x, A.y);
            float c12_pow = BilinearInterpolation(xvec, yvec, interfaces.c12_pow + mi*slice , A.x, A.y);
            float c13_pow = BilinearInterpolation(xvec, yvec, interfaces.c13_pow + mi*slice , A.x, A.y);
            float c14_pow = BilinearInterpolation(xvec, yvec, interfaces.c14_pow + mi*slice , A.x, A.y);
            float c15_pow = BilinearInterpolation(xvec, yvec, interfaces.c15_pow + mi*slice , A.x, A.y);
            float c16_pow = BilinearInterpolation(xvec, yvec, interfaces.c16_pow + mi*slice , A.x, A.y);
            float c22_pow = BilinearInterpolation(xvec, yvec, interfaces.c22_pow + mi*slice , A.x, A.y);
            float c23_pow = BilinearInterpolation(xvec, yvec, interfaces.c23_pow + mi*slice , A.x, A.y);
            float c24_pow = BilinearInterpolation(xvec, yvec, interfaces.c24_pow + mi*slice , A.x, A.y);
            float c25_pow = BilinearInterpolation(xvec, yvec, interfaces.c25_pow + mi*slice , A.x, A.y);
            float c26_pow = BilinearInterpolation(xvec, yvec, interfaces.c26_pow + mi*slice , A.x, A.y);
            float c33_pow = BilinearInterpolation(xvec, yvec, interfaces.c33_pow + mi*slice , A.x, A.y);
            float c34_pow = BilinearInterpolation(xvec, yvec, interfaces.c34_pow + mi*slice , A.x, A.y);
            float c35_pow = BilinearInterpolation(xvec, yvec, interfaces.c35_pow + mi*slice , A.x, A.y);
            float c36_pow = BilinearInterpolation(xvec, yvec, interfaces.c36_pow + mi*slice , A.x, A.y);
            float c44_pow = BilinearInterpolation(xvec, yvec, interfaces.c44_pow + mi*slice , A.x, A.y);
            float c45_pow = BilinearInterpolation(xvec, yvec, interfaces.c45_pow + mi*slice , A.x, A.y);
            float c46_pow = BilinearInterpolation(xvec, yvec, interfaces.c46_pow + mi*slice , A.x, A.y);
            float c55_pow = BilinearInterpolation(xvec, yvec, interfaces.c55_pow + mi*slice , A.x, A.y);
            float c56_pow = BilinearInterpolation(xvec, yvec, interfaces.c56_pow + mi*slice , A.x, A.y);
            float c66_pow = BilinearInterpolation(xvec, yvec, interfaces.c66_pow + mi*slice , A.x, A.y);
            float rho_pow = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice , A.x, A.y);
            float c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice , A.x, A.y);
            float c12_grad = BilinearInterpolation(xvec, yvec, interfaces.c12_grad + mi*slice , A.x, A.y);
            float c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice , A.x, A.y);
            float c14_grad = BilinearInterpolation(xvec, yvec, interfaces.c14_grad + mi*slice , A.x, A.y);
            float c15_grad = BilinearInterpolation(xvec, yvec, interfaces.c15_grad + mi*slice , A.x, A.y);
            float c16_grad = BilinearInterpolation(xvec, yvec, interfaces.c16_grad + mi*slice , A.x, A.y);
            float c22_grad = BilinearInterpolation(xvec, yvec, interfaces.c22_grad + mi*slice , A.x, A.y);
            float c23_grad = BilinearInterpolation(xvec, yvec, interfaces.c23_grad + mi*slice , A.x, A.y);
            float c24_grad = BilinearInterpolation(xvec, yvec, interfaces.c24_grad + mi*slice , A.x, A.y);
            float c25_grad = BilinearInterpolation(xvec, yvec, interfaces.c25_grad + mi*slice , A.x, A.y);
            float c26_grad = BilinearInterpolation(xvec, yvec, interfaces.c26_grad + mi*slice , A.x, A.y);
            float c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice , A.x, A.y);
            float c34_grad = BilinearInterpolation(xvec, yvec, interfaces.c34_grad + mi*slice , A.x, A.y);
            float c35_grad = BilinearInterpolation(xvec, yvec, interfaces.c35_grad + mi*slice , A.x, A.y);
            float c36_grad = BilinearInterpolation(xvec, yvec, interfaces.c36_grad + mi*slice , A.x, A.y);
            float c44_grad = BilinearInterpolation(xvec, yvec, interfaces.c44_grad + mi*slice , A.x, A.y);
            float c45_grad = BilinearInterpolation(xvec, yvec, interfaces.c45_grad + mi*slice , A.x, A.y);
            float c46_grad = BilinearInterpolation(xvec, yvec, interfaces.c46_grad + mi*slice , A.x, A.y);
            float c55_grad = BilinearInterpolation(xvec, yvec, interfaces.c55_grad + mi*slice , A.x, A.y);
            float c56_grad = BilinearInterpolation(xvec, yvec, interfaces.c56_grad + mi*slice , A.x, A.y);
            float c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice , A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice , A.x, A.y);
            var[0]  = c11 + pow(dz, c11_pow) * c11_grad; 
            var[1]  = c12 + pow(dz, c12_pow) * c12_grad;
            var[2]  = c13 + pow(dz, c13_pow) * c13_grad;
            var[3]  = c14 + pow(dz, c14_pow) * c14_grad;
            var[4]  = c15 + pow(dz, c15_pow) * c15_grad;
            var[5]  = c16 + pow(dz, c16_pow) * c16_grad;
            var[6]  = c22 + pow(dz, c22_pow) * c22_grad;
            var[7]  = c23 + pow(dz, c23_pow) * c23_grad;
            var[8]  = c24 + pow(dz, c24_pow) * c24_grad;
            var[9]  = c25 + pow(dz, c25_pow) * c25_grad;
            var[10] = c26 + pow(dz, c26_pow) * c26_grad;
            var[11] = c33 + pow(dz, c33_pow) * c33_grad;
            var[12] = c34 + pow(dz, c34_pow) * c34_grad;
            var[13] = c35 + pow(dz, c35_pow) * c35_grad;
            var[14] = c36 + pow(dz, c36_pow) * c36_grad;
            var[15] = c44 + pow(dz, c44_pow) * c44_grad;
            var[16] = c45 + pow(dz, c45_pow) * c45_grad;
            var[17] = c46 + pow(dz, c46_pow) * c46_grad;
            var[18] = c55 + pow(dz, c55_pow) * c55_grad;
            var[19] = c56 + pow(dz, c56_pow) * c56_grad;
            var[20] = c66 + pow(dz, c66_pow) * c66_grad;
            var[21] = rho + pow(dz, rho_pow) * rho_grad;
        }     
    break;

    case ELASTIC_TTI_BOND: /* 6. c11 c33 c44 c66 c13 rho azimuth dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
        } else {
            float c11 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c44 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            float c66 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float rho = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_grad = BilinearInterpolation(xvec, yvec, interfaces.c11_grad + mi*slice, A.x, A.y);
            float c33_grad = BilinearInterpolation(xvec, yvec, interfaces.c33_grad + mi*slice, A.x, A.y);
            float c44_grad = BilinearInterpolation(xvec, yvec, interfaces.c44_grad + mi*slice, A.x, A.y);
            float c66_grad = BilinearInterpolation(xvec, yvec, interfaces.c66_grad + mi*slice, A.x, A.y);
            float c13_grad = BilinearInterpolation(xvec, yvec, interfaces.c13_grad + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad + mi*slice, A.x, A.y);
            float c11_pow  = BilinearInterpolation(xvec, yvec, interfaces.c11_pow  + mi*slice, A.x, A.y);
            float c33_pow  = BilinearInterpolation(xvec, yvec, interfaces.c33_pow  + mi*slice, A.x, A.y);
            float c44_pow  = BilinearInterpolation(xvec, yvec, interfaces.c44_pow  + mi*slice, A.x, A.y);
            float c66_pow  = BilinearInterpolation(xvec, yvec, interfaces.c66_pow  + mi*slice, A.x, A.y);
            float c13_pow  = BilinearInterpolation(xvec, yvec, interfaces.c13_pow  + mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow  + mi*slice, A.x, A.y);
            float azimuth      = BilinearInterpolation(xvec, yvec, interfaces.azimuth      + mi*slice, A.x, A.y);
            float dip          = BilinearInterpolation(xvec, yvec, interfaces.dip          + mi*slice, A.x, A.y);
            float azimuth_grad = BilinearInterpolation(xvec, yvec, interfaces.azimuth_grad + mi*slice, A.x, A.y);
            float dip_grad     = BilinearInterpolation(xvec, yvec, interfaces.dip_grad     + mi*slice, A.x, A.y);
            float azimuth_pow  = BilinearInterpolation(xvec, yvec, interfaces.azimuth_pow  + mi*slice, A.x, A.y);
            float dip_pow      = BilinearInterpolation(xvec, yvec, interfaces.dip_pow      + mi*slice, A.x, A.y);
            var[0] = c11 + pow(dz, c11_pow)* c11_grad;
            var[1] = c33 + pow(dz, c33_pow)* c33_grad;
            var[2] = c44 + pow(dz, c44_pow)* c44_grad;
            var[3] = c66 + pow(dz, c66_pow)* c66_grad;
            var[4] = c13 + pow(dz, c13_pow)* c13_grad;
            var[5] = rho + pow(dz, rho_pow)* rho_grad;
            var[6] = azimuth + pow(dz, azimuth_pow)* azimuth_grad;
            var[7] = dip     + pow(dz, dip_pow)    * dip_grad;
        }   
    break;

    case ACOUSTIC_ISOTROPIC: /* 7. vp, rho*/
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
        } else {
            float vp_grad  = BilinearInterpolation(xvec, yvec, interfaces.vp_grad + mi*slice, A.x, A.y);
            float vp0      = BilinearInterpolation(xvec, yvec, interfaces.vp      + mi*slice, A.x, A.y);
            float vp_pow   = BilinearInterpolation(xvec, yvec, interfaces.vp_pow  + mi*slice, A.x, A.y);
            float rho0     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float rho_grad = BilinearInterpolation(xvec, yvec, interfaces.rho_grad+ mi*slice, A.x, A.y);
            float rho_pow  = BilinearInterpolation(xvec, yvec, interfaces.rho_pow + mi*slice, A.x, A.y);
            var[0] = vp0  + pow(dz, vp_pow)* vp_grad;
            var[1] = rho0 + pow(dz,rho_pow)*rho_grad;
        }
    break;

    default: // for self-check
        fprintf(stderr,"Error: Unknow meida, please check the code! (for code check, please contact Luqian Jiang)");
        fflush(stderr);
        exit(1);

    }     
} 


/* 
 * For half-grid point, marked the materials number.  
 */
int MediaNumberAtPoint(
        Point3 A,  
        inter_t interfaces) 
{

    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    int     NX = interfaces.NX;
    int     NY = interfaces.NY;
    int     NI = interfaces.NI;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;
    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> elevation(NI, -FLT_MAX);

    for (int i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (int i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }


    for (int ni = 0; ni < NI; ni++) {
        /* Get the elevation for the corresponding location */
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*NX*NY, A.x, A.y);
    }


    // use which material 
    int mi = findLastGreaterEqualIndex(A.z, elevation);
    return mi;
}

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
    float *var3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                
                std::vector<float> var(1, 0.0);

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }


                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ONE_COMPONENT, var);
    
                var3d[indx] = var[0];     
            }
        }
    }

}

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
    float *rho3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(3, 0.0); // vp, vs, rho

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_ISOTROPIC, var);

                /*calculate output lambda and mu */
                float vp = var[0], vs = var[1], rho = var[2];
                mu3d[indx]  = vs*vs*rho;
                lam3d[indx] = vp*vp*rho - 2.0*mu3d[indx];
                rho3d[indx] = rho; 
                
            }
        }
    }
}

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
    float *rho_3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(6, 0.0); /* vph, vpv, vsh, vsv, rho, eta */

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_VTI_PREM, var);
                
                /* Calculate cij */
                float vph = var[0], vpv = var[1];
                float vsh = var[2], vsv = var[3];
                float rho = var[4], eta = var[5];

                c11_3d[indx] = vph*vph*rho;
                c33_3d[indx] = vpv*vpv*rho;
                c66_3d[indx] = vsh*vsh*rho;
                c44_3d[indx] = vsv*vsv*rho;
                rho_3d[indx] = rho; 
                c13_3d[indx] = eta*(c11_3d[indx] - 2.0*c44_3d[indx]);
            }
        }
    }
}

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
    float *rho_3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(6, 0.0); /* vp0 (alpha), vs0 (beta), epsilon, delta, gamma, rho */

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_VTI_THOMSEN, var);
                
                /* Calculate cij */
                float vp0 = var[0], vs0 = var[1];
                float epsilon = var[2], delta = var[3];
                float gamma = var[4], rho = var[5];
                float c33 = vp0*vp0*rho, c44 = vs0*vs0*rho;
                rho_3d[indx] = rho; 
                c33_3d[indx] = c33;
                c44_3d[indx] = c44;                
                c66_3d[indx] = 2.0*gamma  *c44 + c44;
                c11_3d[indx] = 2.0*epsilon*c33 + c33;
                c13_3d[indx] = sqrt( 2.0*delta*c33*(c33-c44) + (c33-c44)*(c33-c44) ) - c44;
            }
        }
    }
}

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
    float *rho_3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(6, 0.0); /* c11 c33 c44 c66 c13 rho*/

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }
                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_VTI_CIJ, var);
                
                /* Calculate cij */
                c11_3d[indx] = var[0];
                c33_3d[indx] = var[1];
                c44_3d[indx] = var[2];
                c66_3d[indx] = var[3];
                c13_3d[indx] = var[4];
                rho_3d[indx] = var[5];
            }
        }
    }
}


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
    float *rho_3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(8, 0.0); /* vp0, vs0, epsilon, delta, gamma, rho, azimuth, dip */

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_TTI_THOMSEN, var);
                
                /* Calculate cij */
                float vp0 = var[0], vs0 = var[1];
                float epsilon = var[2], delta = var[3];
                float gamma = var[4], rho = var[5];
                float azimuth = var[6];
                float dip = var[7];
                float c33 = vp0*vp0*rho;
                float c44 = vs0*vs0*rho;
                float c66 = 2.0*gamma  *c44 + c44;
                float c11 = 2.0*epsilon*c33 + c33;
                float c13 = sqrt( 2.0*delta*c33*(c33-c44) + (c33-c44)*(c33-c44) ) - c44;
                rho_3d[indx] = rho; 

                BondTransform(c11, c11-2*c66, c13, 0, 0, 0,  // cij
                              c11, c13, 0, 0, 0,             // c2j: 23-26
                              c33, 0, 0, 0,                  // c3j
                              c44, 0, 0,                     // c4j
                              c44, 0, c66, var[7], var[6],
                              c11_3d[indx], c12_3d[indx], c13_3d[indx], c14_3d[indx], c15_3d[indx], c16_3d[indx],
                              c22_3d[indx], c23_3d[indx], c24_3d[indx], c25_3d[indx], c26_3d[indx], c33_3d[indx],
                              c34_3d[indx], c35_3d[indx], c36_3d[indx], c44_3d[indx], c45_3d[indx], c46_3d[indx],
                              c55_3d[indx], c56_3d[indx], c66_3d[indx]);

            }
        }
    }
}

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
    float *rho_3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(8, 0.0);  /* c11 c33 c44 c66 c13 rho azimuth dip */

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }
                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_TTI_BOND, var);
                
                /* Calculate cij */
                float c11 = var[0];
                float c33 = var[1];
                float c44 = var[2];
                float c66 = var[3];
                float c13 = var[4];
                rho_3d[indx] = var[5];

                BondTransform(c11, c11-2*c66, c13, 0, 0, 0,  // cij
                              c11, c13, 0, 0, 0,             // c2j: 23-26
                              c33, 0, 0, 0,                  // c3j
                              c44, 0, 0,                     // c4j
                              c44, 0, c66, var[7], var[6],
                              c11_3d[indx], c12_3d[indx], c13_3d[indx], c14_3d[indx], c15_3d[indx], c16_3d[indx],
                              c22_3d[indx], c23_3d[indx], c24_3d[indx], c25_3d[indx], c26_3d[indx], c33_3d[indx],
                              c34_3d[indx], c35_3d[indx], c36_3d[indx], c44_3d[indx], c45_3d[indx], c46_3d[indx],
                              c55_3d[indx], c56_3d[indx], c66_3d[indx]);
            }
        }
    }
}

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
    float *rho_3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                /* c11 c12 c13 c14 c15 c16 c22 c23 c24 c25 c26 c33 c34 c35 c36 c44 c45 c46 c55 c56 c66 rho*/
                std::vector<float> var(22, 0.0); 

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }
                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_TTI_CIJ, var);
                
                /* Calculate cij */
                c11_3d[indx] = var[0];
                c12_3d[indx] = var[1];
                c13_3d[indx] = var[2];
                c14_3d[indx] = var[3];
                c15_3d[indx] = var[4];
                c16_3d[indx] = var[5];
                c22_3d[indx] = var[6];
                c23_3d[indx] = var[7];
                c24_3d[indx] = var[8];
                c25_3d[indx] = var[9];
                c26_3d[indx] = var[10];
                c33_3d[indx] = var[11];
                c34_3d[indx] = var[12];
                c35_3d[indx] = var[13];
                c36_3d[indx] = var[14];
                c44_3d[indx] = var[15];
                c45_3d[indx] = var[16];
                c46_3d[indx] = var[17];
                c55_3d[indx] = var[18];
                c56_3d[indx] = var[19];
                c66_3d[indx] = var[20];
                rho_3d[indx] = var[21];
            }
        }
    }
}

/* 8. assign the parameter directly (use the local values): isotropic, acoustic */
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
    float *rho3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                std::vector<float> var(2, 0.0); // vp, rho

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == 1) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == 3) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignLayerMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ACOUSTIC_ISOTROPIC, var);

                /*calculate output lambda and mu */
                float vp = var[0], rho = var[1];
                kappa[indx] = vp*vp*rho;
                rho3d[indx] = rho; 
                
            }
        }
    }
}


// TODO jlq: temporarily only for curvilinear grid!
// assign the parameter by volume arithmetic and harmonic averaging method
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
    float *rho3d) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;

    // assign the local value first.
    parametrization_el_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, lam3d, mu3d, rho3d);
    
    // jlq: modify to a function, in: 2nd pointer, grid_type, return a grid.
    // For equivalent medium parameterization method
    float *Hx        = new float [siz_volume]; 
    float *Hy        = new float [siz_volume]; 
    float *Hz        = new float [siz_volume]; 
    int   *MaterNum  = new int[siz_volume];

    for (size_t i = 0; i < siz_volume; i++) {
        Hx[i] = Gridx[i];
        Hy[i] = Gridy[i];
        Hz[i] = Gridz[i];
    }


    GenerateHalfGrid(nx, ny, nz, Gridx, Gridy, Gridz, Hx, Hy, Hz); 

    // mark the interface number at the half grid.
    for (size_t i = 1; i < siz_volume; i++) {
        MaterNum[i] = MediaNumberAtPoint(Point3(Hx[i], Hy[i], Hz[i]), interfaces); 
    }

     
    for (size_t k = 1; k < nz; k++) {
        for (size_t j = 1; j < ny; j++) {
            for (size_t i = 1; i < nx; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 
                // check if the mesh have different values;
                std::vector<int> v(8);
                // clockwise: conducive to debugging
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                // There is more than one medium value in the half-grid mesh, 
                if ( NumOfValues(v, NI) > 1) {
                    Mesh3 M( Point3( Hx[indx-1-siz_line-siz_slice], Hy[indx-1-siz_line-siz_slice], Hz[indx-1-siz_line-siz_slice] ),
                             Point3( Hx[indx  -siz_line-siz_slice], Hy[indx  -siz_line-siz_slice], Hz[indx  -siz_line-siz_slice] ),
                             Point3( Hx[indx           -siz_slice], Hy[indx           -siz_slice], Hz[indx           -siz_slice] ),
                             Point3( Hx[indx-1         -siz_slice], Hy[indx-1         -siz_slice], Hz[indx-1         -siz_slice] ),
                             Point3( Hx[indx-1-siz_line], Hy[indx-1-siz_line], Hz[indx-1-siz_line] ),
                             Point3( Hx[indx  -siz_line], Hy[indx  -siz_line], Hz[indx  -siz_line] ),
                             Point3( Hx[indx           ], Hy[indx           ], Hz[indx           ] ),
                             Point3( Hx[indx-1         ], Hy[indx-1         ], Hz[indx-1         ] ));

                    // recalculate the material value of the point
                    float vol_rho   = 0.0;
                    float har_kappa = 0.0;
                    float har_mu    = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3);
                        AssignLayerMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var);

                        float vp = var[0], vs = var[1], rho = var[2];
                        float mu = vs*vs*rho;
                        float lambda = vp*vp*rho - 2.0*mu;

                        vol_rho   += rho;
                        har_kappa += 1.0/(lambda + 2.0/3.0*mu);
                        har_mu    += 1.0/mu;
                    }

                    har_mu    = nsg*1.0/har_mu;
                    har_kappa = nsg*1.0/har_kappa;
                    vol_rho   = vol_rho/nsg; 

                    lam3d[indx] = har_kappa - 2.0/3.0*har_mu;
                    mu3d[indx]  = har_mu;
                    rho3d[indx] = vol_rho;

                    delete []SubGrid;
                }

            }
        }
    }

    delete [] Hx;
    delete [] Hy;
    delete [] Hz;
    delete [] MaterNum;
}

/*=================== for C call =======================*/
//---- 0. one component
// grid type 1
int media_layer2model_onecmp(float *var3d,
                             const float *x3d, // 1d
                             const float *y3d, // 1d
                             const float *z3d, // 1d
                             int nx,
                             int ny, 
                             int nz,
                             int grid_type, 
                             const char *in_var_file,
                             const char *average_method) 
{
    inter_t interfaces;
    bool first_read = true;

    /*Read interface file*/
    read_interface_file(in_var_file, first_read, &interfaces, 
        &interfaces.var, &interfaces.var_grad, &interfaces.var_pow);
   
    if (strcmp(average_method, "loc") == 0) {
        parametrization_oncmp_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, interfaces, var3d);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        isotropic_har(nx, ny, nz, x3d, y3d, z3d, interfaces, var3d);
//    } else if (strcmp(equivalent_medium_method, "ari") == 0) {      //arithemtic
//        isotropic_ari(nx, ny, nz, x3d, y3d, z3d, interfaces, var3d);
    } else {                                                        // default = loc
        parametrization_oncmp_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, interfaces, var3d);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.var;
    delete [] interfaces.var_pow;
    delete [] interfaces.var_grad;

    return 0;
}

//---- 1. elastic isotropic
// grid type 1
int media_layer2model_el_iso(
        float *lam3d,
        float *mu3d,
        float *rho3d,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
        int nx,
        int ny,
        int nz,
        int grid_type, 
        const char *in_rho_file,
        const char *in_vp_file,
        const char *in_vs_file,
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vp, rho, vs */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_vp_file, !first_read, &interfaces,
        &interfaces.vp, &interfaces.vp_grad, &interfaces.vp_pow);

    read_interface_file(in_vs_file, !first_read, &interfaces,
        &interfaces.vs, &interfaces.vs_grad, &interfaces.vs_pow);      

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d);
    } else if (strcmp(equivalent_medium_method, "har") == 0 && grid_type == 3) {
        parametrization_el_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d);
    } else { //default
        parametrization_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.vp;
    delete [] interfaces.vs;
    delete [] interfaces.rho_grad;
    delete [] interfaces.vp_grad;
    delete [] interfaces.vs_grad;
    delete [] interfaces.rho_pow;
    delete [] interfaces.vp_pow;
    delete [] interfaces.vs_pow;   
 
    return 0; 
}

//--- 2. elastic vti_prem
// grid type 1
int media_layer2model_el_vti_prem(
        float *rho,
        float *c11,
        float *c33,
        float *c44,
        float *c66,
        float *c13,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: rho, vpv/h, vsv/h, vs */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_vph_file, !first_read, &interfaces,
        &interfaces.vph, &interfaces.vph_grad, &interfaces.vph_pow);

    read_interface_file(in_vsh_file, !first_read, &interfaces,
        &interfaces.vsh, &interfaces.vsh_grad, &interfaces.vsh_pow);      

    read_interface_file(in_vpv_file, !first_read, &interfaces,
        &interfaces.vpv, &interfaces.vpv_grad, &interfaces.vpv_pow);

    read_interface_file(in_vsv_file, !first_read, &interfaces,
        &interfaces.vsv, &interfaces.vsv_grad, &interfaces.vsv_pow); 

    read_interface_file(in_eta_file, !first_read, &interfaces,
        &interfaces.eta, &interfaces.eta_grad, &interfaces.eta_pow); 

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_vti_prem_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c44, c66, c13, rho);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        parametrization_el_vti_prem_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, c11, c33, c44, c66, c13, rho);
    } else { //default
        parametrization_el_vti_prem_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c44, c66, c13, rho);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.vpv;
    delete [] interfaces.vsv;
    delete [] interfaces.vpv;
    delete [] interfaces.vsv;
    delete [] interfaces.eta;
    delete [] interfaces.rho_grad;
    delete [] interfaces.vpv_grad;
    delete [] interfaces.vsv_grad;
    delete [] interfaces.vpv_grad;
    delete [] interfaces.vsv_grad;
    delete [] interfaces.eta_grad;  
    delete [] interfaces.rho_pow;
    delete [] interfaces.vpv_pow;
    delete [] interfaces.vsv_pow;
    delete [] interfaces.vpv_pow;
    delete [] interfaces.vsv_pow;
    delete [] interfaces.eta_pow;

    return 0; 
}

//--- 3. elastic vti_thomsen
// grid type 1
int media_layer2model_vti_thomsen(
        float *rho,
        float *c11,
        float *c33,
        float *c44,
        float *c66,
        float *c13,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vpv, vsv, epsilon, delta, gamma, rho */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_vpv_file, !first_read, &interfaces,
        &interfaces.vpv, &interfaces.vpv_grad, &interfaces.vpv_pow);

    read_interface_file(in_vsv_file, !first_read, &interfaces,
        &interfaces.vsv, &interfaces.vsv_grad, &interfaces.vsv_pow); 

    read_interface_file(in_epsilon_file, !first_read, &interfaces,
        &interfaces.epsilon, &interfaces.epsilon_grad, &interfaces.epsilon_pow); 

    read_interface_file(in_delta_file, !first_read, &interfaces,
        &interfaces.delta, &interfaces.delta_grad, &interfaces.delta_pow); 

    read_interface_file(in_gamma_file, !first_read, &interfaces,
        &interfaces.gamma, &interfaces.gamma_grad, &interfaces.gamma_pow); 

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_vti_thomsen_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c44, c66, c13, rho);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        parametrization_el_vti_prem_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, c11, c33, c44, c66, c13, rho);
    } else { //default
        parametrization_el_vti_thomsen_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c44, c66, c13, rho);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.vpv;
    delete [] interfaces.vsv;
    delete [] interfaces.rho_grad;
    delete [] interfaces.vpv_grad;
    delete [] interfaces.vsv_grad;
    delete [] interfaces.rho_pow;
    delete [] interfaces.vpv_pow;
    delete [] interfaces.vsv_pow;  
    delete [] interfaces.gamma;
    delete [] interfaces.gamma_pow;
    delete [] interfaces.gamma_grad;
    delete [] interfaces.epsilon;
    delete [] interfaces.epsilon_pow;
    delete [] interfaces.epsilon_grad;
    delete [] interfaces.delta;
    delete [] interfaces.delta_grad;
    delete [] interfaces.delta_pow;

    return 0; 
}

//--- 4. elastic vti_cij
// grid type 1
int media_layer2model_vti_cij(
        float *rho,
        float *c11,
        float *c33,
        float *c44,
        float *c66,
        float *c13,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vpv, vsv, epsilon, delta, gamma, rho */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_c11_file, !first_read, &interfaces,
        &interfaces.c11, &interfaces.c11_grad, &interfaces.c11_pow);

    read_interface_file(in_c33_file, !first_read, &interfaces,
        &interfaces.c33, &interfaces.c33_grad, &interfaces.c33_pow); 

    read_interface_file(in_c44_file, !first_read, &interfaces,
        &interfaces.c44, &interfaces.c44_grad, &interfaces.c44_pow); 

    read_interface_file(in_c66_file, !first_read, &interfaces,
        &interfaces.c66, &interfaces.c66_grad, &interfaces.c66_pow); 

    read_interface_file(in_c13_file, !first_read, &interfaces,
        &interfaces.c13, &interfaces.c13_grad, &interfaces.c13_pow); 

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_vti_cij_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c44, c66, c13, rho);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        parametrization_el_vti_cij_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, c11, c33, c44, c66, c13, rho);
    } else { //default
        parametrization_el_vti_cij_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c33, c44, c66, c13, rho);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.c11;
    delete [] interfaces.c33;
    delete [] interfaces.c44;
    delete [] interfaces.c66;
    delete [] interfaces.c13;
    delete [] interfaces.rho_grad;
    delete [] interfaces.c11_grad;
    delete [] interfaces.c33_grad;
    delete [] interfaces.c44_grad;
    delete [] interfaces.c66_grad;
    delete [] interfaces.c13_grad;
    delete [] interfaces.rho_pow;
    delete [] interfaces.c11_pow;
    delete [] interfaces.c33_pow;
    delete [] interfaces.c44_pow;
    delete [] interfaces.c66_pow;
    delete [] interfaces.c13_pow;

    return 0; 
}

//--- 5. elastic tti_thomsen
// grid type 1
int media_layer2model_tti_thomsen(
        float *rho,
        float *c11, float *c12, float *c13,
        float *c14, float *c15, float *c16,
        float *c22, float *c23, float *c24,
        float *c25, float *c26, float *c33,
        float *c34, float *c35, float *c36,
        float *c44, float *c45, float *c46,
        float *c55, float *c56, float *c66,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vpv, vsv, epsilon, delta, gamma, rho */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_vpv_file, !first_read, &interfaces,
        &interfaces.vpv, &interfaces.vpv_grad, &interfaces.vpv_pow);

    read_interface_file(in_vsv_file, !first_read, &interfaces,
        &interfaces.vsv, &interfaces.vsv_grad, &interfaces.vsv_pow); 

    read_interface_file(in_epsilon_file, !first_read, &interfaces,
        &interfaces.epsilon, &interfaces.epsilon_grad, &interfaces.epsilon_pow); 

    read_interface_file(in_delta_file, !first_read, &interfaces,
        &interfaces.delta, &interfaces.delta_grad, &interfaces.delta_pow); 

    read_interface_file(in_gamma_file, !first_read, &interfaces,
        &interfaces.gamma, &interfaces.gamma_grad, &interfaces.gamma_pow); 

    read_interface_file(in_azimuth_file, !first_read, &interfaces,
        &interfaces.azimuth, &interfaces.azimuth_grad, &interfaces.azimuth_pow); 

    read_interface_file(in_dip_file, !first_read, &interfaces,
        &interfaces.dip, &interfaces.dip_grad, &interfaces.dip_pow); 

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_tti_thomsen_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//       parametrization_el_tti_thomsen_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
//            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
    } else { //default
                parametrization_el_tti_thomsen_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.vpv;
    delete [] interfaces.vsv;
    delete [] interfaces.rho_grad;
    delete [] interfaces.vpv_grad;
    delete [] interfaces.vsv_grad;
    delete [] interfaces.rho_pow;
    delete [] interfaces.vpv_pow;
    delete [] interfaces.vsv_pow;  
    delete [] interfaces.gamma;
    delete [] interfaces.gamma_pow;
    delete [] interfaces.gamma_grad;
    delete [] interfaces.epsilon;
    delete [] interfaces.epsilon_pow;
    delete [] interfaces.epsilon_grad;
    delete [] interfaces.delta;
    delete [] interfaces.delta_grad;
    delete [] interfaces.delta_pow;
    delete [] interfaces.azimuth;
    delete [] interfaces.azimuth_pow;
    delete [] interfaces.azimuth_grad;
    delete [] interfaces.dip;
    delete [] interfaces.dip_pow;
    delete [] interfaces.dip_grad;

    return 0; 
}

//--- 6. elastic tti_bond
// grid type 1
int media_layer2model_tti_bond(
        float *rho,
        float *c11, float *c12, float *c13,
        float *c14, float *c15, float *c16,
        float *c22, float *c23, float *c24,
        float *c25, float *c26, float *c33,
        float *c34, float *c35, float *c36,
        float *c44, float *c45, float *c46,
        float *c55, float *c56, float *c66,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vpv, vsv, epsilon, delta, gamma, rho */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_c11_file, !first_read, &interfaces,
        &interfaces.c11, &interfaces.c11_grad, &interfaces.c11_pow);

    read_interface_file(in_c33_file, !first_read, &interfaces,
        &interfaces.c33, &interfaces.c33_grad, &interfaces.c33_pow); 

    read_interface_file(in_c44_file, !first_read, &interfaces,
        &interfaces.c44, &interfaces.c44_grad, &interfaces.c44_pow); 

    read_interface_file(in_c66_file, !first_read, &interfaces,
        &interfaces.c66, &interfaces.c66_grad, &interfaces.c66_pow); 

    read_interface_file(in_c13_file, !first_read, &interfaces,
        &interfaces.c13, &interfaces.c13_grad, &interfaces.c13_pow); 

    read_interface_file(in_azimuth_file, !first_read, &interfaces,
        &interfaces.azimuth, &interfaces.azimuth_grad, &interfaces.azimuth_pow); 

    read_interface_file(in_dip_file, !first_read, &interfaces,
        &interfaces.dip, &interfaces.dip_grad, &interfaces.dip_pow); 

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_tti_bond_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        parametrization_el_vti_cij_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, c11, c33, c44, c66, c13, rho);
    } else { //default
        parametrization_el_tti_bond_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.c11;
    delete [] interfaces.c33;
    delete [] interfaces.c44;
    delete [] interfaces.c66;
    delete [] interfaces.c13;
    delete [] interfaces.rho_grad;
    delete [] interfaces.c11_grad;
    delete [] interfaces.c33_grad;
    delete [] interfaces.c44_grad;
    delete [] interfaces.c66_grad;
    delete [] interfaces.c13_grad;
    delete [] interfaces.rho_pow;
    delete [] interfaces.c11_pow;
    delete [] interfaces.c33_pow;
    delete [] interfaces.c44_pow;
    delete [] interfaces.c66_pow;
    delete [] interfaces.c13_pow;
    delete [] interfaces.azimuth;
    delete [] interfaces.azimuth_pow;
    delete [] interfaces.azimuth_grad;
    delete [] interfaces.dip;
    delete [] interfaces.dip_pow;
    delete [] interfaces.dip_grad;

    return 0; 
}

//--- 7. elastic tti_cij
// grid type 1
int media_layer2model_tti_cij(
        float *rho,
        float *c11, float *c12, float *c13,
        float *c14, float *c15, float *c16,
        float *c22, float *c23, float *c24,
        float *c25, float *c26, float *c33,
        float *c34, float *c35, float *c36,
        float *c44, float *c45, float *c46,
        float *c55, float *c56, float *c66,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vpv, vsv, epsilon, delta, gamma, rho */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    // c1j
    read_interface_file(in_c11_file, !first_read, &interfaces,
        &interfaces.c11, &interfaces.c11_grad, &interfaces.c11_pow);

    read_interface_file(in_c12_file, !first_read, &interfaces,
        &interfaces.c12, &interfaces.c12_grad, &interfaces.c12_pow);

    read_interface_file(in_c13_file, !first_read, &interfaces,
        &interfaces.c13, &interfaces.c13_grad, &interfaces.c13_pow);

    read_interface_file(in_c14_file, !first_read, &interfaces,
        &interfaces.c14, &interfaces.c14_grad, &interfaces.c14_pow);

    read_interface_file(in_c15_file, !first_read, &interfaces,
        &interfaces.c15, &interfaces.c15_grad, &interfaces.c15_pow);

    read_interface_file(in_c16_file, !first_read, &interfaces,
        &interfaces.c16, &interfaces.c16_grad, &interfaces.c16_pow);

    // c2j
    read_interface_file(in_c22_file, !first_read, &interfaces,
        &interfaces.c22, &interfaces.c22_grad, &interfaces.c22_pow);

    read_interface_file(in_c23_file, !first_read, &interfaces,
        &interfaces.c23, &interfaces.c23_grad, &interfaces.c23_pow);

    read_interface_file(in_c24_file, !first_read, &interfaces,
        &interfaces.c24, &interfaces.c24_grad, &interfaces.c24_pow);

    read_interface_file(in_c25_file, !first_read, &interfaces,
        &interfaces.c25, &interfaces.c25_grad, &interfaces.c25_pow);

    read_interface_file(in_c26_file, !first_read, &interfaces,
        &interfaces.c26, &interfaces.c26_grad, &interfaces.c26_pow);

    // c3j
    read_interface_file(in_c33_file, !first_read, &interfaces,
        &interfaces.c33, &interfaces.c33_grad, &interfaces.c33_pow);

    read_interface_file(in_c34_file, !first_read, &interfaces,
        &interfaces.c34, &interfaces.c34_grad, &interfaces.c34_pow);

    read_interface_file(in_c35_file, !first_read, &interfaces,
        &interfaces.c35, &interfaces.c35_grad, &interfaces.c35_pow);

    read_interface_file(in_c36_file, !first_read, &interfaces,
        &interfaces.c36, &interfaces.c36_grad, &interfaces.c36_pow);

    // c4j
    read_interface_file(in_c44_file, !first_read, &interfaces,
        &interfaces.c44, &interfaces.c44_grad, &interfaces.c44_pow);

    read_interface_file(in_c45_file, !first_read, &interfaces,
        &interfaces.c45, &interfaces.c45_grad, &interfaces.c45_pow);

    read_interface_file(in_c46_file, !first_read, &interfaces,
        &interfaces.c46, &interfaces.c46_grad, &interfaces.c46_pow);

    // c5j
    read_interface_file(in_c55_file, !first_read, &interfaces,
        &interfaces.c55, &interfaces.c55_grad, &interfaces.c55_pow);

    read_interface_file(in_c56_file, !first_read, &interfaces,
        &interfaces.c56, &interfaces.c56_grad, &interfaces.c56_pow);

    read_interface_file(in_c66_file, !first_read, &interfaces,
        &interfaces.c66, &interfaces.c66_grad, &interfaces.c66_pow);

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_el_tti_bond_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        parametrization_el_vti_cij_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, c11, c33, c44, c66, c13, rho);
    } else { //default
        parametrization_el_tti_cij_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
            c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.c11;
    delete [] interfaces.c12;
    delete [] interfaces.c13;
    delete [] interfaces.c14;
    delete [] interfaces.c15;
    delete [] interfaces.c16;
    delete [] interfaces.c22;
    delete [] interfaces.c23;
    delete [] interfaces.c24;
    delete [] interfaces.c25;
    delete [] interfaces.c26;
    delete [] interfaces.c33;
    delete [] interfaces.c34;
    delete [] interfaces.c35;
    delete [] interfaces.c36;
    delete [] interfaces.c44;
    delete [] interfaces.c45;
    delete [] interfaces.c46;
    delete [] interfaces.c55;
    delete [] interfaces.c56;
    delete [] interfaces.c66;
    delete [] interfaces.rho_grad;
    delete [] interfaces.c11_grad;
    delete [] interfaces.c12_grad;
    delete [] interfaces.c13_grad;
    delete [] interfaces.c14_grad;
    delete [] interfaces.c15_grad;
    delete [] interfaces.c16_grad;
    delete [] interfaces.c22_grad;
    delete [] interfaces.c23_grad;
    delete [] interfaces.c24_grad;
    delete [] interfaces.c25_grad;
    delete [] interfaces.c26_grad;
    delete [] interfaces.c33_grad;
    delete [] interfaces.c34_grad;
    delete [] interfaces.c35_grad;
    delete [] interfaces.c36_grad;
    delete [] interfaces.c44_grad;
    delete [] interfaces.c45_grad;
    delete [] interfaces.c46_grad;
    delete [] interfaces.c55_grad;
    delete [] interfaces.c56_grad;
    delete [] interfaces.c66_grad;
    delete [] interfaces.rho_grad;
    delete [] interfaces.c11_pow;
    delete [] interfaces.c12_pow;
    delete [] interfaces.c13_pow;
    delete [] interfaces.c14_pow;
    delete [] interfaces.c15_pow;
    delete [] interfaces.c16_pow;
    delete [] interfaces.c22_pow;
    delete [] interfaces.c23_pow;
    delete [] interfaces.c24_pow;
    delete [] interfaces.c25_pow;
    delete [] interfaces.c26_pow;
    delete [] interfaces.c33_pow;
    delete [] interfaces.c34_pow;
    delete [] interfaces.c35_pow;
    delete [] interfaces.c36_pow;
    delete [] interfaces.c44_pow;
    delete [] interfaces.c45_pow;
    delete [] interfaces.c46_pow;
    delete [] interfaces.c55_pow;
    delete [] interfaces.c56_pow;
    delete [] interfaces.c66_pow;
    delete [] interfaces.rho_pow;

    return 0; 
}

//---- 7. acoustic isotropic
// grid type 1
int media_layer2model_ac_iso(
        float *rho3d,
        float *kappa3d,
        const float *x3d, // 1d
        const float *y3d, // 1d
        const float *z3d, // 1d
        int nx,
        int ny,
        int nz,
        int grid_type, 
        const char *in_rho_file,
        const char *in_vp_file,
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;

    /* Read interface file: vp, rho, vs */
    read_interface_file(in_rho_file, first_read, &interfaces,
        &interfaces.rho, &interfaces.rho_grad, &interfaces.rho_pow);

    read_interface_file(in_vp_file, !first_read, &interfaces,
        &interfaces.vp, &interfaces.vp_grad, &interfaces.vp_pow);

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_ac_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, kappa3d, rho3d);
//    } else if (strcmp(equivalent_medium_method, "har") == 0) {
//        parametrization_el_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
//            interfaces, lam3d, mu3d, rho3d);
    } else { //default
        parametrization_ac_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, kappa3d, rho3d);
    }

    delete [] interfaces.elevation;
    delete [] interfaces.rho;
    delete [] interfaces.vp;
    delete [] interfaces.rho_grad;
    delete [] interfaces.vp_grad;
    delete [] interfaces.rho_pow;
    delete [] interfaces.vp_pow;
 
    return 0; 
}
