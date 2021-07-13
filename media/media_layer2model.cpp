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
#include "media_geometry3d.hpp"
#include "media_layer2model.hpp"
#include "media_read_interface_file.hpp"
#include "media_utility.hpp"

//using namespace std;


int AssignMediaPara2Point(
    int ix, int iy, int iz, 
    Point3 A, 
    int NI, 
    Interfaces *interfaces,
    float &vp, 
    float &vs,
    float &rho)
{
    float MINX = interfaces[0].MINX;
    float MINY = interfaces[0].MINY;
    float   DX = interfaces[0].DX;
    float   DY = interfaces[0].DY;
    size_t  NX = interfaces[0].NX;
    size_t  NY = interfaces[0].NY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;
    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> altitude(NI, -FLT_MAX);
    for (size_t i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (size_t i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }

    // for each interface, interpolate the altitude of the position
    //  where the Point A is located.
    for (int ni = 0; ni < NI; ni++) {
        // Out of the INTERFACE MESH area. 
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            fprintf(stderr,"\033[47;31mGrid(%d, %d, %d) is out of the interfaces mesh range, "\
                           "please check the interfaces file!\033[0m\n", ix, iy, iz);
            fflush(stderr);
            exit(1);
        }


        // Get the altitude for the corresponding location
        altitude[ni] = BilinearInterpolation(XVEC, YVEC, interfaces[ni].altitude, A.x, A.y);
    }

    // use which material 
    int mi = findNearestGreaterIndex(A.z, altitude);

    if (mi == -1) {
        // z-axis is positive downward
        fprintf(stderr,"\033[47;31mz-location of Grid(%d, %d, %d) is smaller than the given altitude in " \
            "the interfaces file, please check the interfaces file or the grid!\033[0m\n", ix, iy, iz);
        fflush(stderr);
        exit(1);    
    }


    // Bilinear
    float vp0      = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vp      , A.x, A.y);
    float vp_grad  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vp_grad , A.x, A.y);
    float vs0      = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vs      , A.x, A.y);
    float vs_grad  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vs_grad , A.x, A.y);
    float rho0     = BilinearInterpolation(XVEC, YVEC, interfaces[mi].rho     , A.x, A.y);
    float rho_grad = BilinearInterpolation(XVEC, YVEC, interfaces[mi].rho_grad, A.x, A.y);

    vp  = vp0  + (altitude[mi] - A.z)* vp_grad;
    vs  = vs0  + (altitude[mi] - A.z)* vs_grad;
    rho = rho0 + (altitude[mi] - A.z)*rho_grad;

    return 0;
}

/* 
 * For half-grid point, marked the materials number.  
 */
int MediaNumberAtPoint(
    Point3 A, 
    int NI, 
    Interfaces *interfaces) 
{

    float MINX = interfaces[0].MINX;
    float MINY = interfaces[0].MINY;
    float   DX = interfaces[0].DX;
    float   DY = interfaces[0].DY;
    int     NX = interfaces[0].NX;
    int     NY = interfaces[0].NY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;
    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> altitude(NI, -FLT_MAX);

    for (int i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (int i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }


    for (int ni = 0; ni < NI; ni++) {
        // Out of the INTERFACE MESH area.
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            continue;
        }
        // Get the altitude for the corresponding location
        altitude[ni] = BilinearInterpolation(XVEC, YVEC, interfaces[ni].altitude, A.x, A.y);
    }


    // use which material 
    int mi = findNearestGreaterIndex(A.z, altitude);
    return mi;
}


// assign the parameter directly (use the local values)
void isotropic_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NI,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d )
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                float vp = 0.0, vs = 0.0, rho = 1.0;
                size_t indx =  i + j * siz_line + k * siz_slice;

                AssignMediaPara2Point(i, j, k,
                    Point3(Gridx[indx], Gridy[indx], Gridz[indx]),
                    NI, interfaces, vp, vs, rho);

                mu3d[indx]  = vs*vs*rho;
                lam3d[indx] = vp*vp*rho - 2.0*mu3d[indx];
                rho3d[indx] = rho; 
                
            }
        }
    }

}

// assign the parameter by volume arithmetic and harmonic averaging method
void isotropic_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NI,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;

    // assign the local value first.
    isotropic_loc(nx,ny, nz, Gridx, Gridy, Gridz,
        NI, interfaces, lam3d, mu3d, rho3d);
    
    // For equivalent medium parameterization method
    float *Hx = new float [siz_volume]; 
    float *Hy = new float [siz_volume]; 
    float *Hz = new float [siz_volume]; 
    int   *MaterNum  = new int[siz_volume];

    for (size_t i = 0; i < siz_volume; i++) {
        Hx[i] = Gridx[i];
        Hy[i] = Gridy[i];
        Hz[i] = Gridz[i];
    }

    GenerateHalfGrid(nx, ny, nz, Gridx, Gridy, Gridz, Hx, Hy, Hz); 

    // mark the interface number at the half grid.
    for (size_t i = 1; i < siz_volume; i++) {
        MaterNum[i] = MediaNumberAtPoint(Point3(Hx[i], Hy[i], Hz[i]), NI, interfaces); 
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

                        float vp = 0.0, vs = 0.0, rho = 0.0;
                        AssignMediaPara2Point(i, j, k, SubGrid[isg],
                            NI, interfaces, vp, vs, rho);
                        float mu =  vs*vs*rho;
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


    delete []Hx;
    delete []Hy;
    delete []Hz;
    delete []MaterNum;
}

// read interface file and medium parameterization
void media_el_iso_layer2model(
    float *lam3d,
    float *mu3d,
    float *rho3d,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    const char *interfaces_file,
    const char *equivalent_medium_method) 
{

    int NI = 0;
    Interfaces *interfaces = nullptr;

    /*Read interface file*/
    read_interface_file(interfaces_file, NI, &interfaces); 

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        isotropic_loc(nx, ny, nz, x3d, y3d, z3d, NI,
                      interfaces, lam3d, mu3d, rho3d);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        isotropic_har(nx, ny, nz, x3d, y3d, z3d, NI,
                      interfaces, lam3d, mu3d, rho3d);
    } else { //default
        isotropic_loc(nx, ny, nz, x3d, y3d, z3d, NI,
                      interfaces, lam3d, mu3d, rho3d);
    }


    for (int ni = 0; ni < NI; ni++) {
        delete []interfaces[ni].altitude;
        delete []interfaces[ni].rho;
        delete []interfaces[ni].vp;
        delete []interfaces[ni].vs;
        delete []interfaces[ni].rho_grad;
        delete []interfaces[ni].vp_grad;
        delete []interfaces[ni].vs_grad;
    }

    delete []interfaces;
    
}

// duplicated fuction of ../forward/md_el_iso.c file
//void md_el_iso_rho_to_slow(float *rho, size_t siz_volume) {
//    for (size_t i = 0; i < siz_volume; i++) {
//        if (rho[i] > EPS) {
//            rho[i] = 1.0 / rho[i];
//        } else {
//            rho[i] = 0.0;
//        }    
//    }
//}
