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

int AssignLayerMediaPara2Point_onecmp(
    int ix, int iy, int iz, 
    Point3 A,  
    inter_t interfaces,
    float &var) // use interface.var
{
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    size_t  NX = interfaces.NX;
    size_t  NY = interfaces.NY;
    size_t  NI = interfaces.NI;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;

    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> elevation(NI, -FLT_MAX);
    for (size_t i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (size_t i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }

    // for each interface, interpolate the elevation of the position
    //  where the Point A is located.
    for (int ni = 0; ni < NI; ni++) {
        // Out of the INTERFACE MESH area. 
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            fprintf(stderr,"Error: Grid(%d, %d, %d) is out of the interfaces mesh range, "\
                           "please check the interfaces file!\n", ix, iy, iz);
            fflush(stderr);
            exit(1);
        }

        // Get the elevation for the corresponding location
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*NX*NY, A.x, A.y);
    }

    // use which material 
    int mi = findNearestGreaterIndex(A.z, elevation);

    if (mi == -1) {
        // z-axis is positive downward
        fprintf(stderr,"Warning: z-location of Grid(%d, %d, %d) is smaller than the given elevation in " \
            "the interfaces file, please check the interfaces file or the grid! \n", ix, iy, iz);
        fflush(stderr);

        var = BilinearInterpolation(XVEC, YVEC, interfaces.var, A.x, A.y);

    } else {
        // Bilinear
        float var0      = BilinearInterpolation(XVEC, YVEC, interfaces.var      + mi*NX*NY , A.x, A.y);
        float var_grad  = BilinearInterpolation(XVEC, YVEC, interfaces.var_grad + mi*NX*NY , A.x, A.y);
        float var_pow   = BilinearInterpolation(XVEC, YVEC, interfaces.var_pow  + mi*NX*NY , A.x, A.y);
    
        var  = var0  + pow((elevation[mi] - A.z), var_pow)* var_grad;
    }
    return 0;
}

int AssignLayerMediaPara2Point_el_iso(
    int ix, int iy, int iz, 
    Point3 A, 
    inter_t interfaces,
    float &vp, float &vs, float &rho)
{
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    size_t  NX = interfaces.NX;
    size_t  NY = interfaces.NY;
    size_t  NI = interfaces.NI;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;

    std::vector<float> XVEC(NX), YVEC(NY);
    std::vector<float> elevation(NI, -FLT_MAX);
    for (size_t i = 0; i < NX; i++) {
        XVEC[i] = MINX + DX*i;
    }
    for (size_t i = 0; i < NY; i++) {
        YVEC[i] = MINY + DY*i;
    }

    // for each interface, interpolate the elevation of the position
    //  where the Point A is located.
    for (int ni = 0; ni < NI; ni++) {
        // Out of the INTERFACE MESH area. 
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            fprintf(stderr,"Error: Grid(%d, %d, %d) is out of the interfaces mesh range, "\
                           "please check the interfaces file!\n", ix, iy, iz);
            fflush(stderr);
            exit(1);
        }


        // Get the elevation for the corresponding location
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*NX*NY, A.x, A.y);

    }

    // use which material 
    int mi = findNearestGreaterIndex(A.z, elevation);

    if (mi == -1) {
        // z-axis is positive downward
        fprintf(stderr,"Warning: z-location of Grid(%d, %d, %d) is higher than the given elevation in " \
            "the interfaces file, assigned by the top-interface value! \n", ix, iy, iz);
        fflush(stderr);

        vp  = BilinearInterpolation(XVEC, YVEC, interfaces.vp  , A.x, A.y);
        vs  = BilinearInterpolation(XVEC, YVEC, interfaces.vs  , A.x, A.y);
        rho = BilinearInterpolation(XVEC, YVEC, interfaces.rho , A.x, A.y);


    } else {
        // Bilinear
        float vp_grad  = BilinearInterpolation(XVEC, YVEC, interfaces.vp_grad + mi*NX*NY, A.x, A.y);
        float vp0      = BilinearInterpolation(XVEC, YVEC, interfaces.vp +    + mi*NX*NY, A.x, A.y);
        float vp_pow   = BilinearInterpolation(XVEC, YVEC, interfaces.vp_pow  + mi*NX*NY, A.x, A.y);
        float vs0      = BilinearInterpolation(XVEC, YVEC, interfaces.vs      + mi*NX*NY, A.x, A.y);
        float vs_grad  = BilinearInterpolation(XVEC, YVEC, interfaces.vs_grad + mi*NX*NY, A.x, A.y);
        float vs_pow   = BilinearInterpolation(XVEC, YVEC, interfaces.vs_pow  + mi*NX*NY, A.x, A.y);
        float rho0     = BilinearInterpolation(XVEC, YVEC, interfaces.rho     + mi*NX*NY, A.x, A.y);
        float rho_grad = BilinearInterpolation(XVEC, YVEC, interfaces.rho_grad+ mi*NX*NY, A.x, A.y);
        float rho_pow  = BilinearInterpolation(XVEC, YVEC, interfaces.rho_pow + mi*NX*NY, A.x, A.y);
    
        vp  = vp0  + pow((elevation[mi]-A.z), vp_pow)* vp_grad;
        vs  = vs0  + pow((elevation[mi]-A.z), vs_pow)* vs_grad;
        rho = rho0 + pow((elevation[mi]-A.z),rho_pow)*rho_grad;
    }
    return 0;
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
        // Out of the INTERFACE MESH area.
        if (A.x < MINX || A.x > MAXX || A.y < MINY || A.y > MAXY) {
            continue;
        }
        // Get the elevation for the corresponding location
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*NX*NY, A.x, A.y);
    }


    // use which material 
    int mi = findNearestGreaterIndex(A.z, elevation);
    return mi;
}

// assign the parameter directly (use the local values)
void parametrization_oncmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, // 1: xyz1d; 2: xy1d, z3d; 3: xyz3d.
    inter_t interfaces,
    float *var3d)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = grid_type==1? k : indx;
                size_t indx_y = grid_type==3? indx : j;
                size_t indx_x = grid_type==3? indx : i;

                AssignLayerMediaPara2Point_onecmp(i, j, k,
                    Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                    interfaces, var3d[indx]);

            }
        }
    }
}

// assign the parameter directly (use the local values)
void parametrization_el_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, // 1: xyz1d; 2: xy1d, z3d; 3: xyz3d.
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
    
                float vp = 0.0, vs = 0.0, rho = 1.0;
                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = grid_type==1? k : indx;
                size_t indx_y = grid_type==3? indx : j;
                size_t indx_x = grid_type==3? indx : i;

                AssignLayerMediaPara2Point_el_iso(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, vp, vs, rho);
    
                mu3d[indx]  = vs*vs*rho;
                lam3d[indx] = vp*vp*rho - 2.0*mu3d[indx];
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
    parametrization_el_iso_loc(nx,ny, nz, Gridx, Gridy, Gridz,
        grid_type, interfaces, lam3d, mu3d, rho3d);
    
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

                        float vp = 0.0, vs = 0.0, rho = 0.0;
                        AssignLayerMediaPara2Point_el_iso(i, j, k, SubGrid[isg],
                            interfaces, vp, vs, rho);
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
        const char *equivalent_medium_method) 
{

    inter_t interfaces;
    bool first_read = true;
    int grid_type = 3;

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
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_el_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d);
    } else { //default
        parametrization_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
            interfaces, lam3d, mu3d, rho3d);
    }


    
    delete []interfaces.elevation;
    delete []interfaces.rho;
    delete []interfaces.vp;
    delete []interfaces.vs;
    delete []interfaces.rho_grad;
    delete []interfaces.vp_grad;
    delete []interfaces.vs_grad;
    delete []interfaces.rho_pow;
    delete []interfaces.vp_pow;
    delete []interfaces.vs_pow;   

    return 0; 
}


/* Not complete yet! */
int media_layer2model_curv_onecmp(float *var3d,
                                  const float *x3d, // nx*ny*nz
                                  const float *y3d,
                                  const float *z3d,
                                  int nx,
                                  int ny, 
                                  int nz,
                                  const char *in_var_file,
                                  const char *average_method) 
{
    inter_t interfaces;
    bool first_read = true;

    /*Read interface file*/
    read_interface_file(in_var_file, first_read, &interfaces, 
        &interfaces.var, &interfaces.var_grad, &interfaces.var_pow);

/*   Not complete yet
    if (strcmp(average_method, "loc") == 0) {
        isotropic_loc(nx, ny, nz, x3d, y3d, z3d, interfaces, var3d);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        isotropic_har(nx, ny, nz, x3d, y3d, z3d, interfaces, var3d);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {      //arithemtic
        isotropic_ari(nx, ny, nz, x3d, y3d, z3d, interfaces, var3d);
    } else {                                                        // default = loc
        isotropic_loc(nx, ny, nz, x3d, y3d, z3d, interfaces, var3d);
    }
*/

    delete []interfaces.elevation;
    delete []interfaces.var;
    delete []interfaces.var_pow;
    delete []interfaces.var_grad;

    return 0;
}