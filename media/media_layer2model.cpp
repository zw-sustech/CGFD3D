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
#include "media_interpolation.hpp"
#include "media_read_interface_file.hpp"

//using namespace std;

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
    // half grid point is inside the grid, 
    // so the number of grid point in one direction is np-1
    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;

    for (size_t k = 0; k < nz-1; k++) {
        for (size_t j = 0; j < ny-1; j++) {
            for (size_t i = 0; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice;
                // clockwise: conducive to debugging
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

// How many different values in the vector, 
//  NI is the upper limitation of the v[i].
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


int AssignMediaPara2Point(
    int ix, int iy, int iz, 
    Point3 A, 
    int NI, 
    Interfaces *interfaces,
    float &vp, 
    float &vs,
    float &rho)
{
    int flag = 0;
    std::vector<float> depth(NI, -FLT_MAX);

    // for each interface, interpolate the depth of the position
    //  where the Point A is located.
    for (int ni = 0; ni < NI; ni++) {
        float MINX = interfaces[ni].MINX;
        float MINY = interfaces[ni].MINY;
        float   DX = interfaces[ni].DX;
        float   DY = interfaces[ni].DY;
        size_t  NX = interfaces[ni].NX;
        size_t  NY = interfaces[ni].NY;
        float MAXX = MINX + (NX-1)*DX;  
        float MAXY = MINY + (NY-1)*DY;

        // Out of the INTERFACE MESH area.
        // TODO: For we change the interface file by one mesh,
        //       an error will be reported if out of range.  
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            fprintf(stderr,"\033[47;31mGrid(%d, %d, %d) is out of the interfaces mesh range, "\
                           "please check the interfaces file!\033[0m\n", ix, iy, iz);
            fflush(stderr);
            exit(1);
            //flag++;
            //continue;
        }

        std::vector<float> XVEC(NX), YVEC(NY);
        for (size_t i = 0; i < NX; i++)
            XVEC[i] = MINX + DX*i;
        for (size_t i = 0; i < NY; i++)
            YVEC[i] = MINY + DY*i;

        // Geth the depth for the corresponding location
        depth[ni] = BilinearInterpolation(XVEC, YVEC, interfaces[ni].depth, A.x, A.y);
    }

 /*
    if (flag == NI) {
        fprintf(stderr,"\033[47;31mGrid(%d, %d, %d) is out of the interfaces mesh range, "\
            "please check the interfaces file!\033[0m\n", ix, iy, iz);
        fflush(stderr);
        exit(1);    
    }
*/
    // use which material 
    int mi = findNearestGreaterIndex(A.z, depth);

    if (mi == -1) {
        // z-axis is positive downward
        fprintf(stderr,"\033[47;31mz-location of Grid(%d, %d, %d) is smaller than the given depth in " \
            "the interfaces file, please check the interfaces file or the grid!\033[0m\n", ix, iy, iz);
        fflush(stderr);
        exit(1);    
    }

    // Get the Material info of this point from interfaces[mi]
    std::vector<float> XVEC(interfaces[mi].NX), YVEC(interfaces[mi].NY);
    for (size_t i = 0; i < interfaces[mi].NX; i++)
        XVEC[i] = interfaces[mi].MINX + interfaces[mi].DX*i;
    for (size_t i = 0; i < interfaces[mi].NY; i++)
        YVEC[i] = interfaces[mi].MINY + interfaces[mi].DY*i;

    // Bilinear or trilinear ???
    float vp0      = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vp      , A.x, A.y);
    float vp_grad  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vp_grad , A.x, A.y);
    float vs0      = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vs      , A.x, A.y);
    float vs_grad  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vs_grad , A.x, A.y);
    float rho0     = BilinearInterpolation(XVEC, YVEC, interfaces[mi].rho     , A.x, A.y);
    float rho_grad = BilinearInterpolation(XVEC, YVEC, interfaces[mi].rho_grad, A.x, A.y);

    vp  = vp0  + (A.z - depth[mi])* vp_grad;
    vs  = vs0  + (A.z - depth[mi])* vs_grad;
    rho = rho0 + (A.z - depth[mi])*rho_grad;

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

    std::vector<float> depth(NI, -FLT_MAX);
    for (int ni = 0; ni < NI; ni++) {
        float MINX = interfaces[ni].MINX;
        float MINY = interfaces[ni].MINY;
        float   DX = interfaces[ni].DX;
        float   DY = interfaces[ni].DY;
        int     NX = interfaces[ni].NX;
        int     NY = interfaces[ni].NY;
        float MAXX = MINX + (NX-1)*DX;  
        float MAXY = MINY + (NY-1)*DY;

        // Out of the INTERFACE MESH area.
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            continue;
        }

        std::vector<float> XVEC(NX), YVEC(NY);
        for (int i = 0; i < NX; i++)
            XVEC[i] = MINX + DX*i;
        for (int i = 0; i < NY; i++)
            YVEC[i] = MINY + DY*i;

        // Geth the depth for the corresponding location
        depth[ni] = BilinearInterpolation(XVEC, YVEC, interfaces[ni].depth, A.x, A.y);
    }


    // use which material 
    int mi = findNearestGreaterIndex(A.z, depth);
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
                size_t indx = indx =  i + j * siz_line + k * siz_slice;

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

     
    for (size_t k = 1; k < nz-1; k++) {
        for (size_t j = 1; j < ny-1; j++) {
            for (size_t i = 1; i < nx-1; i++) {
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
    size_t siz_line,
    size_t siz_slice,
    size_t siz_volume, 
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
