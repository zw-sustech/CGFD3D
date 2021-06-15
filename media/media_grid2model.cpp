/******************************************************************************
 *
 * This function is used to discretize a given grid model to a calculation grid.
 *
 * Authors: Luqian Jiang <jianglq@mail.ustc.edu.cn>
 *          Wei Zhang <zhangwei@sustech.edu.cn>
 *
 * Copyright (c) 2021 zwlab
 *
 * ChangeLog: 
 *    06/2021: Created by Luqian Jiang 
 *
 *******************************************************************************/
#include <iostream>
#include <vector>
#include <string.h>
#include <math.h>  
#include "media_geometry3d.hpp"
#include "media_grid2model.hpp"
#include "media_read_interface_file.hpp"
#include "media_utility.hpp"

//using namespace std;


int AssignGridMediaPara2Point(
    int ix, int iy, int iz, 
    Point3 A, 
    int NI, 
    Interfaces *interfaces,
    float &vp, 
    float &vs,
    float &rho)
{

    std::vector<float> altitude(NI, -FLT_MAX);
    /* All the given grid mesh is the same */
    float MINX = interfaces[0].MINX;
    float MINY = interfaces[0].MINY;
    float   DX = interfaces[0].DX;
    float   DY = interfaces[0].DY;
    size_t  NX = interfaces[0].NX;
    size_t  NY = interfaces[0].NY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;

    std::vector<float> XVEC(NX), YVEC(NY);
    for (size_t i = 0; i < NX; i++)
        XVEC[i] = MINX + DX*i;
    for (size_t i = 0; i < NY; i++)
        YVEC[i] = MINY + DY*i;

    /* 
     * For each interface, interpolate the altitude of the position
     *   to get where the Point A is located in xoy plane.
     */
    for (int ni = 0; ni < NI; ni++) {

        /* Out of the INTERFACE MESH area.*/
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            fprintf(stderr,"\033[47;31mGrid(%d, %d, %d) is out of the interfaces mesh range, "\
                           "please check the interfaces file!\033[0m\n", ix, iy, iz);
            fflush(stderr);
            exit(1);
        }

        /* Get the altitude for the corresponding xoy location */
        altitude[ni] = BilinearInterpolation(XVEC, YVEC, interfaces[ni].altitude, A.x, A.y);
    }

    /* which material is used, the z-grid is given from top to bottom */
    int mi = findNearestGreaterIndex(A.z, altitude);

    if (mi == -1) {   
        /* The altitude is above the surface, assigned by the 1st layer para. */
        vp  = BilinearInterpolation(XVEC, YVEC, interfaces[0].vp  , A.x, A.y);
        vs  = BilinearInterpolation(XVEC, YVEC, interfaces[0].vs  , A.x, A.y);
        rho = BilinearInterpolation(XVEC, YVEC, interfaces[0].rho , A.x, A.y);        
    } else if (mi == NI-1) { 
        /* The altitude is below the last layer, assigned by the last para */    
        vp  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vp  , A.x, A.y);
        vs  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vs  , A.x, A.y);
        rho = BilinearInterpolation(XVEC, YVEC, interfaces[mi].rho , A.x, A.y);         
    } else { 
        /* Linear interpolation in z-direction */    

        float vp0  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vp  , A.x, A.y);
        float vs0  = BilinearInterpolation(XVEC, YVEC, interfaces[mi].vs  , A.x, A.y);
        float rho0 = BilinearInterpolation(XVEC, YVEC, interfaces[mi].rho , A.x, A.y); 

        float vp1  = BilinearInterpolation(XVEC, YVEC, interfaces[mi+1].vp  , A.x, A.y);
        float vs1  = BilinearInterpolation(XVEC, YVEC, interfaces[mi+1].vs  , A.x, A.y);
        float rho1 = BilinearInterpolation(XVEC, YVEC, interfaces[mi+1].rho , A.x, A.y); 

        /* Init: if A.z is equal to the value in altitude, 
         *  assign the average value.        
         */
        float dis_r = 1.0/2.0;

        if ( fabs(altitude[mi+1]-altitude[mi]) > 1e-6 ) {
            dis_r = (A.z - altitude[mi])/(altitude[mi+1] - altitude[mi]);
        }

        vp  = vp0  + (vp1-vp0)   * dis_r;
        vs  = vs0  + (vs1-vs0)   * dis_r;
        rho = rho0 + (rho1-rho0) * dis_r;
    }

    return 0;
}

/* 
 * For half-grid point, mark the interface number.  
 */
int LayerNumberAtPoint(
    Point3 A, 
    int NL,
    std::vector<int> NGz, 
    Interfaces *interfaces) 
{
    std::vector<float> altitude(NL, -FLT_MAX);
    float MINX = interfaces[0].MINX;
    float MINY = interfaces[0].MINY;
    float   DX = interfaces[0].DX;
    float   DY = interfaces[0].DY;
    int     NX = interfaces[0].NX;
    int     NY = interfaces[0].NY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;

    std::vector<float> XVEC(NX), YVEC(NY);
    for (int i = 0; i < NX; i++)
        XVEC[i] = MINX + DX*i;
    for (int i = 0; i < NY; i++)
        YVEC[i] = MINY + DY*i;

    /* 
     * Mark the interfaces[0] and interfaces[accumulation of NGz[nl]], 
     * nl is form 0 to NL-2.
     */
    int ni = 0;
    for (int nl = 0; nl < NL; nl++) {
        /* 
         * Out of the INTERFACE MESH area. 
         * The half-grid is inside the integer grid,
         * If out of the integer-grid, the program will stop,
         *  so this judgement is useless.
         */
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            continue;
        }

        /* Get the altitude for the corresponding location */
        altitude[nl] = BilinearInterpolation(XVEC, YVEC, interfaces[ni].altitude, A.x, A.y);
        ni += NGz[nl];
    }

    /* Mark the layer number */
    int mi = findNearestGreaterIndex(A.z, altitude);

    /* If the altitude is above the surface, it is assigned by the first layer para. */
    if (mi == -1) mi = 0;

    return mi;
}


/* assign the parameter directly (use the local values) */
void iso_grid_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NL, 
    std::vector <int> NGz,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    int NI = 0;
    for (int i = 0; i < NL; i++) {
        NI += NGz[i];
    }

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                float vp = 0.0, vs = 0.0, rho = 1.0;
                size_t indx = i + j * siz_line + k * siz_slice;

                AssignGridMediaPara2Point(i, j, k,
                    Point3(Gridx[indx], Gridy[indx], Gridz[indx]),
                    NI, interfaces, vp, vs, rho);

                mu3d[indx]  = vs*vs*rho;
                lam3d[indx] = vp*vp*rho - 2.0*vs*vs*rho;
                rho3d[indx] = rho; 
                
            }
        }
    }

}

/* assign the parameter by volume arithmetic and harmonic averaging method */
void iso_grid_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int NL, 
    std::vector <int> NGz,
    Interfaces *interfaces,
    float *lam3d,
    float *mu3d,
    float *rho3d) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;

    int NI = 0;
    for (int i = 0; i < NL; i++) {
        NI += NGz[i];
    }

    /* assign the local value first.*/
    iso_grid_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        NL, NGz, interfaces, lam3d, mu3d, rho3d);
    
    /* For equivalent medium parameterization method */
    float *Hx = new float [siz_volume]; 
    float *Hy = new float [siz_volume]; 
    float *Hz = new float [siz_volume]; 
    int   *MaterNum  = new int[siz_volume];

    /* init the Hx, Hy, and Hz */
    for (size_t i = 0; i < siz_volume; i++) {
        Hx[i] = Gridx[i];
        Hy[i] = Gridy[i];
        Hz[i] = Gridz[i];
    }
    GenerateHalfGrid(nx, ny, nz, Gridx, Gridy, Gridz, Hx, Hy, Hz); 

    /* 
     * Mark the layer number at the half grid,
     *  for the grid media format, just 0 - NL-1 are marked.
     */
    for (size_t i = 0; i < siz_volume; i++) {
        MaterNum[i] = LayerNumberAtPoint(Point3(Hx[i], Hy[i], Hz[i]), NL, NGz, interfaces); 
    }

    /* Loop integer-grid points */ 
    for (size_t k = 1; k < nz; k++) {
        for (size_t j = 1; j < ny; j++) {
            for (size_t i = 1; i < nx; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice; 

                /* Check if the corresponding the half-grid mesh have different values */
                std::vector<int> v(8);
                /* clockwise: conducive to debugging */
                v[0] = MaterNum[indx-1-siz_line-siz_slice];
                v[1] = MaterNum[indx  -siz_line-siz_slice];
                v[2] = MaterNum[indx           -siz_slice];
                v[3] = MaterNum[indx-1         -siz_slice];
                v[4] = MaterNum[indx-1-siz_line];
                v[5] = MaterNum[indx  -siz_line];
                v[6] = MaterNum[indx           ];
                v[7] = MaterNum[indx-1         ];

                /* 
                 * There is more than one medium value in the half-grid mesh, 
                 *  subdivide the mesh.
                 */
                if ( NumOfValues(v, NL) > 1) {
                    Mesh3 M( Point3( Hx[indx-1-siz_line-siz_slice], Hy[indx-1-siz_line-siz_slice], Hz[indx-1-siz_line-siz_slice] ),
                             Point3( Hx[indx  -siz_line-siz_slice], Hy[indx  -siz_line-siz_slice], Hz[indx  -siz_line-siz_slice] ),
                             Point3( Hx[indx           -siz_slice], Hy[indx           -siz_slice], Hz[indx           -siz_slice] ),
                             Point3( Hx[indx-1         -siz_slice], Hy[indx-1         -siz_slice], Hz[indx-1         -siz_slice] ),
                             Point3( Hx[indx-1-siz_line], Hy[indx-1-siz_line], Hz[indx-1-siz_line] ),
                             Point3( Hx[indx  -siz_line], Hy[indx  -siz_line], Hz[indx  -siz_line] ),
                             Point3( Hx[indx           ], Hy[indx           ], Hz[indx           ] ),
                             Point3( Hx[indx-1         ], Hy[indx-1         ], Hz[indx-1         ] ));

                    /* recalculate the material value of the point */
                    float vol_rho   = 0.0;
                    float har_kappa = 0.0;
                    float har_mu    = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    for (int isg = 0; isg < nsg; isg++) {

                        float vp = 0.0, vs = 0.0, rho = 0.0;

                        AssignGridMediaPara2Point(i, j, k, SubGrid[isg],
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
 
                    delete[] SubGrid;
                }

            }
        }
    }


    delete[] Hx;
    delete[] Hy;
    delete[] Hz;
    delete[] MaterNum;
}

/* read interface file and medium parameterization */
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
    const char *equivalent_medium_method) // equivalent medium method: just for 
{

    // NL: n-layer
    int NL = 0; // N-layer
    std::vector<int> NGz;
    Interfaces *interfaces = nullptr;

    /*Read interface file*/
    read_grid_file(grid_file, Xmin, Xmax, Ymin, Ymax, NL, NGz, &interfaces);


    if (strcmp(equivalent_medium_method, "loc") == 0) {

        iso_grid_loc(nx, ny, nz, x3d, y3d, z3d, NL, NGz,
                          interfaces, lam3d, mu3d, rho3d);


    } else if (strcmp(equivalent_medium_method, "har") == 0) {

        iso_grid_har(nx, ny, nz, x3d, y3d, z3d, NL, NGz,
                     interfaces, lam3d, mu3d, rho3d);

    } else { //default

        iso_grid_loc(nx, ny, nz, x3d, y3d, z3d, NL, NGz,
                          interfaces, lam3d, mu3d, rho3d);

    }


    int Ng = 0;
    for (int nl = 0; nl < NL; nl++) {
        Ng += NGz[nl];
    }

    for (int ni = 0; ni < Ng; ni++) {
        delete[] interfaces[ni].altitude;
        delete[] interfaces[ni].rho;
        delete[] interfaces[ni].vp;
        delete[] interfaces[ni].vs;
    }

    delete[] interfaces;
    
}


