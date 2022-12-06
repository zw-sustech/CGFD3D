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
#include <set>
#include "media_geometry3d.hpp"
#include "media_utility.hpp"
#include "media_grid2model.hpp"
#include "media_read_file.hpp"


std::map<int, std::string> md2str = create_md2str_map();

/*============================ for C call ====================================*/
//--- 0. one component
int media_grid2model_onecmp(
    float *var3d,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method,
    int myid)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz;
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, Ymin, Ymax, NL, NGz, &interfaces);

    if (interfaces.media_type != ONE_COMPONENT) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_grid2model_onecmp() only supports one_component,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);        
        fflush(stderr);
        exit(1);
     }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_onecmp_loc(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, var3d, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_onecmp_har(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, var3d, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_onecmp_ari(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, var3d, myid);
    } else { //default
        fprintf(stderr,"Error: Wrong average method %s for one_component.\n", equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }

    return 0;
}

//--- 1. acoustic isotropic
int media_grid2model_ac_iso(
    float *rho3d,
    float *kappa3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method,
    int myid)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, Ymin, Ymax, NL, NGz, &interfaces);

    if (interfaces.media_type != ACOUSTIC_ISOTROPIC) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_grid2model_ac_iso() only supports acoustic_isotropic,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);        
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_ac_iso_loc(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, rho3d, kappa3d, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_ac_iso_har(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, rho3d, kappa3d, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_ac_iso_ari(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, rho3d, kappa3d, myid);
    } else { //default
        fprintf(stderr,"Error: Wrong average method %s for acoustic_isotropic media.\n", 
            equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }

    return 0;
}

//--- 2. elastic isotropic
int media_grid2model_el_iso(
    float *rho3d,
    float *lam3d,
    float *mu3d, 
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method,
    int myid)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, Ymin, Ymax, NL, NGz, &interfaces);

    if (interfaces.media_type != ELASTIC_ISOTROPIC) {
        fprintf(stderr, "Error: media_type=%s is not supported,\n"\
                        "       media_grid2model_el_iso() only supports elastic_isotropic,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);         
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, rho3d, lam3d, mu3d, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_el_iso_har(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, rho3d, lam3d, mu3d, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_el_iso_ari(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, rho3d, lam3d, mu3d, myid);
    } else { //default
        fprintf(stderr,"Error: Wrong parametrization method %s in media_grid2model_el_iso(), " \
                       "       if you want to use tti equivalent medium method, " \
                       "       please call the media_grid2model_el_aniso().\n", 
                       equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }

    return 0;
}


//--- 3. elastic vti
int media_grid2model_el_vti(
    float *rho,
    float *c11, 
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    const float *x3d,
    const float *y3d,
    const float *z3d,
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method,
    int myid)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    /* Read grid media file*/
    read_grid_file(in_media_file, Xmin, Xmax, Ymin, Ymax, NL, NGz, &interfaces);
    int md_type = interfaces.media_type;

    if (md_type != ELASTIC_VTI_PREM && md_type != ELASTIC_VTI_THOMSEN && md_type != ELASTIC_VTI_CIJ) {
        fprintf(stderr, "Error: media_type=%s is not supported in media_grid2model_el_vti(),\n"\
                        "       it only supports elastic_vti_prem, elastic_vti_thomsen, elastic_vti_cij,\n"\
                        "       please check the media_type in %s! \n", 
                md2str[interfaces.media_type].c_str(), in_media_file);    
        fflush(stderr);
        exit(1);
    }

    if (strcmp(equivalent_medium_method, "loc") == 0) {
        parametrization_grid_el_vti_loc(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
    } else if (strcmp(equivalent_medium_method, "har") == 0) {
        parametrization_grid_el_vti_har(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
    } else if (strcmp(equivalent_medium_method, "ari") == 0) {
        parametrization_grid_el_vti_ari(nx, ny, nz, x3d, y3d, z3d, 
           grid_type, NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
    } else { //default
        fprintf(stderr,"Error: Wrong parametrization method %s in media_grid2model_el_vti(), " \
                       "       if you want to use tti equivalent medium method, " \
                       "       please call the media_grid2model_el_aniso().\n", 
                       equivalent_medium_method);
        fflush(stderr);
        exit(1);
    }

    return 0;
}

//--- 4. elastic anisotropic/TTI
int media_grid2model_el_aniso(
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
    int nx,
    int ny,
    int nz,
    float Xmin, float Xmax,
    float Ymin, float Ymax, 
    int grid_type,
    const char *in_media_file,
    const char *equivalent_medium_method,
    int myid)  
{

    int NL = 0; // N-layer
    std::vector<int> NGz; // N z-grid in each layer
    inter_t interfaces;

    size_t siz_volume = nx*ny*nz;

    // Read grid media file
    read_grid_file(in_media_file, Xmin, Xmax, Ymin, Ymax, NL, NGz, &interfaces);

    int md_type = interfaces.media_type;

    // the function does not support one component and acoustic wave
    if (md_type == ONE_COMPONENT){
        fprintf(stderr, "Error: media_type=one_component is not supported in media_grid2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_media_file);        
        fflush(stderr);
        exit(1);
    } else if (md_type == ACOUSTIC_ISOTROPIC){
        fprintf(stderr, "Error: media_type=acoustic_isotropic is not supported in media_grid2model_el_aniso(),\n"\
                        "       please check the media_type of %s! \n", in_media_file);         
        fflush(stderr);
        exit(1);
    }

    //- isotropic: loc, har, tti
    if (md_type == ELASTIC_ISOTROPIC) {
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_grid_el_iso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, rho, c13, c44, myid);
            for (size_t i = 0; i < siz_volume; ++i) {
                c11[i] = c13[i] + 2.0*c44[i]; 
                c22[i] = c11[i]; c33[i] = c11[i]; 
                c12[i] = c13[i]; c23[i] = c13[i];
                c55[i] = c44[i]; c66[i] = c44[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_grid_el_iso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, rho, c13, c44, myid);
            for (size_t i = 0; i < siz_volume; ++i) {
                c11[i] = c13[i] + 2.0*c44[i]; 
                c22[i] = c11[i]; c33[i] = c11[i]; 
                c12[i] = c13[i]; c23[i] = c13[i];
                c55[i] = c44[i]; c66[i] = c44[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_grid_el_iso_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, rho, c13, c44, myid);
            for (size_t i = 0; i < siz_volume; ++i) {
                c11[i] = c13[i] + 2.0*c44[i]; 
                c22[i] = c11[i]; c33[i] = c11[i]; 
                c12[i] = c13[i]; c23[i] = c13[i];
                c55[i] = c44[i]; c66[i] = c44[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
// TODO: tti equivalent medium method for iso
        } else {
            fprintf(stderr, "Error: no such equivalent_medium_method: %s!\n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    // vti: loc, har, ari    
    } else if (md_type == ELASTIC_VTI_PREM || md_type == ELASTIC_VTI_THOMSEN || md_type == ELASTIC_VTI_CIJ) {

        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_grid_el_vti_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
            for (size_t i = 0; i < siz_volume; i++) {
                c12[i] = c11[i]-2.0*c66[i];
                c22[i] = c11[i];
                c23[i] = c13[i];
                c44[i] = c55[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }

        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_grid_el_vti_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
            for (size_t i = 0; i < siz_volume; i++) {
                c12[i] = c11[i]-2.0*c66[i];
                c22[i] = c11[i];
                c23[i] = c13[i];
                c44[i] = c55[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }

        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_grid_el_vti_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);

            for (size_t i = 0; i < siz_volume; i++) {
                c12[i] = c11[i]-2.0*c66[i];
                c22[i] = c11[i];
                c23[i] = c13[i];
                c44[i] = c55[i];
                c14[i] = 0.0; c24[i] = 0.0;  c34[i] = 0.0;  
                c15[i] = 0.0; c25[i] = 0.0;  c35[i] = 0.0; 
                c16[i] = 0.0; c26[i] = 0.0;  c36[i] = 0.0;
                c45[i] = 0.0; c46[i] = 0.0;  c56[i] = 0.0;
            }

        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
// TODO: tti equivalent medium for vti
        } else {
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    } else {
        if (strcmp(equivalent_medium_method, "loc") == 0) {
            parametrization_grid_el_aniso_loc(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, myid);
        } else if (strcmp(equivalent_medium_method, "har") == 0) {
            parametrization_grid_el_aniso_har(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, myid);
        } else if (strcmp(equivalent_medium_method, "ari") == 0) {
            parametrization_grid_el_aniso_ari(nx, ny, nz, x3d, y3d, z3d, grid_type, 
                NL, NGz, interfaces, c11, c12, c13, c14, c15, c16, c22, c23, c24, c25, c26, 
                c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, rho, myid);
        } else if(strcmp(equivalent_medium_method,"tti") == 0) {
// TODO: tti equivalent medium for tti
        } else { //default
            fprintf(stderr, "Error: no such equivalent_medium_method: %s! \n", equivalent_medium_method);        
            fflush(stderr);
            exit(1);        
        }
    }

    return 0;
}
/*=======================================================================================================*/

int AssignGridMediaPara2Point(
    size_t ix, size_t iy, size_t iz, 
    Point3 A, 
    inter_t &interfaces,
    int media_type,
    std::vector<float> &var, 
    std::vector<int> &NGz, int NL)
{
    /* All the given grid mesh is the same */
    size_t  NI = interfaces.NI;
    size_t  NX = interfaces.NX;
    size_t  NY = interfaces.NY;
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    float MAXX = MINX + (NX-1)*DX;  
    float MAXY = MINY + (NY-1)*DY;
    int isPointOnInter = 0;

    // if out of the MEDIA GRID MESH area, exit !
    PrintIsPointOutOfInterfaceRange(A, ix, iy, iz, MINX, MAXX, MINY, MAXY);

    std::set<int> layer_indx;
    int indx = 0;
    for (int i = 0; i < NL-1; i++) {
        indx += NGz[i];
        layer_indx.insert(indx);
    }

    std::vector<float> XVEC(NX), YVEC(NY);
    for (size_t i = 0; i < NX; i++)
        XVEC[i] = MINX + DX*i;
    for (size_t i = 0; i < NY; i++)
        YVEC[i] = MINY + DY*i;

    /* 
     * For each interface, interpolate the elevation of the position
     *   to get where the Point A is located in xoy plane.
     */
    size_t inter_slice = NX*NY;
    std::vector<float> elevation(NI, -FLT_MAX);
    for (int ni = 0; ni < NI; ni++) {
        /* Get the elevation for the corresponding xoy location */
        elevation[ni] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*inter_slice, A.x, A.y);
    }

    /* find which material is used, the z-grid is given from top to bottom */
    //int mi = findNearestGreaterIndex(A.z, elevation);
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    if (layer_indx.count(mi) && isEqual(A.z, elevation[mi])) 
        isPointOnInter = 1;

    CalPointValue_grid(media_type, interfaces, inter_slice, XVEC, YVEC, NI, layer_indx, A, elevation, mi, var);

    return isPointOnInter;
}


//- Calculate the value of the point for different media type (to avoid multiple geometric calculations) 
//   for grid2model
void CalPointValue_grid(int media_type, 
                   inter_t &interfaces,
                   size_t slice, 
                   std::vector<float> &xvec,  /* interface mesh */
                   std::vector<float> &yvec,
                   int NI,
                   std::set<int> &layer_indx, 
                   Point3 &A,
                   std::vector<float> &elevation, /*the elevation of point A at the projection position of the interface mesh. */
                   int mi,
                   std::vector<float> &var)
{
    float dis_r = 1.0/2.0;
    if ( fabs(elevation[mi+1]-elevation[mi]) > 1e-6 ) {
        dis_r = (A.z - elevation[mi])/(elevation[mi+1] - elevation[mi]);
    }

    switch(media_type)
    {
    case ONE_COMPONENT: /* 0. var */
        /* If grid_z > elevation of top_interface, it given by the medium of top non-zero thickness layer */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.var + mi*slice, A.x, A.y); 
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.var + mi*slice, A.x, A.y); 
        } else { 
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float var0 = BilinearInterpolation(xvec, yvec, interfaces.var + mi*slice, A.x, A.y);
            float var1 = BilinearInterpolation(xvec, yvec, interfaces.var + (mi+1)*slice, A.x, A.y);
            var[0]  = var0  + (var1-var0) * dis_r;
        }
    break;

    case ELASTIC_ISOTROPIC: /*1. rho, vp, vs*/
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs  + mi*slice, A.x, A.y);
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs  + mi*slice, A.x, A.y);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            // special treatment: point is on the media grid, average of upper and lower media
            float vp0  = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            float vs0  = BilinearInterpolation(xvec, yvec, interfaces.vs  + mi*slice, A.x, A.y);
            float rho0 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            mi += 1;
            float vp1  = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            float vs1  = BilinearInterpolation(xvec, yvec, interfaces.vs  + mi*slice, A.x, A.y);
            float rho1 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = vp0  + (vp1-vp0)*dis_r;
            var[2] = vs0  + (vs1-vs0)*dis_r;
            var[0] = rho0 + (rho1-rho0)*dis_r;
            if (var[1] <= sqrt(2) * var[2]) {
              fprintf(stderr, "Error: vp must larger than sqrt(2)*vs!\n");
              fprintf(stderr, "       please check the md3grd file!\n");
              fflush(stderr); exit(1);
            }
        }
    break;

    case ELASTIC_VTI_PREM: /*2. rho, vph, vpv, vsh, vsv, eta */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);     
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
        }
        else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float vph0 = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            float vpv0 = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            float vsh0 = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            float vsv0 = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            float rho0 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float eta0 = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
            mi += 1;
            float vph1 = BilinearInterpolation(xvec, yvec, interfaces.vph + mi*slice, A.x, A.y);
            float vpv1 = BilinearInterpolation(xvec, yvec, interfaces.vpv + mi*slice, A.x, A.y);
            float vsh1 = BilinearInterpolation(xvec, yvec, interfaces.vsh + mi*slice, A.x, A.y);
            float vsv1 = BilinearInterpolation(xvec, yvec, interfaces.vsv + mi*slice, A.x, A.y);
            float rho1 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float eta1 = BilinearInterpolation(xvec, yvec, interfaces.eta + mi*slice, A.x, A.y);
            var[0] = rho0 + (rho1-rho0)*dis_r;
            var[1] = vph0 + (vph1-vph0)*dis_r;
            var[2] = vpv0 + (vpv1-vpv0)*dis_r;
            var[3] = vsh0 + (vsh1-vsh0)*dis_r;
            var[4] = vsv0 + (vsv1-vsv0)*dis_r;
            var[5] = eta0 + (eta1-eta0)*dis_r;
        }   
    break;

    case ELASTIC_VTI_THOMSEN: /*3. rho, vp0, vs0, epsilon, delta, gamma */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho    + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp0    + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs0    + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.epsilon+ mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.delta  + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.gamma  + mi*slice , A.x, A.y);    
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho    + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp0    + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs0    + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.epsilon+ mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.delta  + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.gamma  + mi*slice , A.x, A.y);    
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float rho0     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float vp00     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
            float vs00     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
            float epsil0   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);  // epsilon
            float delta0  = BilinearInterpolation(xvec, yvec, interfaces.delta    + mi*slice, A.x, A.y);
            float gamma0   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            mi += 1;
            float rho1     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float vp01     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
            float vs01     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
            float epsil1   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);  // epsilon
            float delta1   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma1   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);            
            var[0] = rho0   + ( rho1   - rho0  )*dis_r;
            var[1] = vp00   + ( vp01   - vp00  )*dis_r;
            var[2] = vs00   + ( vs01   - vs00  )*dis_r;
            var[3] = epsil0 + ( epsil1 - epsil0)*dis_r;
            var[4] = delta0 + ( delta1 - delta0)*dis_r;
            var[5] = gamma0 + ( gamma1 - gamma0)*dis_r;
        }   
    break;

    case ELASTIC_VTI_CIJ: /*4. rho c11 c33 c55 c66 c13 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);   
        } else  if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);   
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float rho_0 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_0 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33_0 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c55_0 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c66_0 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13_0 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            mi += 1;
            float rho_1 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float c11_1 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33_1 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c55_1 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c66_1 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13_1 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = c11_0 + (c11_1 - c11_0)*dis_r;
            var[2] = c33_0 + (c33_1 - c33_0)*dis_r;
            var[3] = c55_0 + (c55_1 - c55_0)*dis_r;
            var[4] = c66_0 + (c66_1 - c66_0)*dis_r;
            var[5] = c13_0 + (c13_1 - c13_0)*dis_r;
        }   
    break;

    case ELASTIC_TTI_THOMSEN: /*5. rho, vp0, vs0, epsilon, delta, gamma, azimuth, dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice , A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice , A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice , A.x, A.y);
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice , A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice , A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice , A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice , A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice , A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice , A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice , A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice , A.x, A.y);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float rho_0     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float vp0_0     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
            float vs0_0     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
            float dip_0     = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
            float epsil_0   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);
            float delta_0   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma_0   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            float azimu_0   = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            mi += 1;
            float rho_1     = BilinearInterpolation(xvec, yvec, interfaces.rho     + mi*slice, A.x, A.y);
            float vp0_1     = BilinearInterpolation(xvec, yvec, interfaces.vp0     + mi*slice, A.x, A.y);
            float vs0_1     = BilinearInterpolation(xvec, yvec, interfaces.vs0     + mi*slice, A.x, A.y);
            float dip_1     = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
            float epsil_1   = BilinearInterpolation(xvec, yvec, interfaces.epsilon + mi*slice, A.x, A.y);
            float delta_1   = BilinearInterpolation(xvec, yvec, interfaces.delta   + mi*slice, A.x, A.y);
            float gamma_1   = BilinearInterpolation(xvec, yvec, interfaces.gamma   + mi*slice, A.x, A.y);
            float azimu_1   = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);            

            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = vp0_0 + (vp0_1 - vp0_0)*dis_r;
            var[2] = vs0_0 + (vs0_1 - vs0_0)*dis_r;
            var[3] = epsil_0 + (epsil_1 - epsil_0)*dis_r;
            var[4] = delta_0 + (delta_1 - delta_0)*dis_r;
            var[5] = gamma_0 + (gamma_1 - gamma_0)*dis_r;
            var[6] = azimu_0 + (azimu_1 - azimu_0)*dis_r;
            var[7] = dip_0 + (dip_1 - dip_0)*dis_r;    
        }   
    break;

    case ELASTIC_ANISO_CIJ: /* 7. rho c11 c12 c13 c14 c15 c16 c22 c23 c24 c25 c26 c33 c34 c35 c36 c44 c45 c46 c55 c56 c66 */
         if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0]  = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1]  = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2]  = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            var[3]  = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[4]  = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            var[5]  = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            var[6]  = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            var[7]  = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            var[8]  = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            var[9]  = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            var[10] = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            var[11] = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            var[12] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[13] = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            var[14] = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            var[15] = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            var[16] = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            var[17] = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            var[18] = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            var[19] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[20] = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            var[21] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
        } else if (mi == NI-1) {
            var[0]  = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1]  = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2]  = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            var[3]  = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[4]  = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            var[5]  = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            var[6]  = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            var[7]  = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            var[8]  = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            var[9]  = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            var[10] = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            var[11] = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            var[12] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[13] = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            var[14] = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            var[15] = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            var[16] = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            var[17] = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            var[18] = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            var[19] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[20] = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            var[21] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float c11_0 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c12_0 = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            float c13_0 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float c14_0 = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            float c15_0 = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            float c16_0 = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            float c22_0 = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            float c23_0 = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            float c24_0 = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            float c25_0 = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            float c26_0 = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            float c33_0 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c34_0 = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            float c35_0 = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            float c36_0 = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            float c44_0 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            float c45_0 = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            float c46_0 = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            float c55_0 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c56_0 = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            float c66_0 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float rho_0 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            mi += 1;
            float c11_1 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c12_1 = BilinearInterpolation(xvec, yvec, interfaces.c12 + mi*slice, A.x, A.y);
            float c13_1 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float c14_1 = BilinearInterpolation(xvec, yvec, interfaces.c14 + mi*slice, A.x, A.y);
            float c15_1 = BilinearInterpolation(xvec, yvec, interfaces.c15 + mi*slice, A.x, A.y);
            float c16_1 = BilinearInterpolation(xvec, yvec, interfaces.c16 + mi*slice, A.x, A.y);
            float c22_1 = BilinearInterpolation(xvec, yvec, interfaces.c22 + mi*slice, A.x, A.y);
            float c23_1 = BilinearInterpolation(xvec, yvec, interfaces.c23 + mi*slice, A.x, A.y);
            float c24_1 = BilinearInterpolation(xvec, yvec, interfaces.c24 + mi*slice, A.x, A.y);
            float c25_1 = BilinearInterpolation(xvec, yvec, interfaces.c25 + mi*slice, A.x, A.y);
            float c26_1 = BilinearInterpolation(xvec, yvec, interfaces.c26 + mi*slice, A.x, A.y);
            float c33_1 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c34_1 = BilinearInterpolation(xvec, yvec, interfaces.c34 + mi*slice, A.x, A.y);
            float c35_1 = BilinearInterpolation(xvec, yvec, interfaces.c35 + mi*slice, A.x, A.y);
            float c36_1 = BilinearInterpolation(xvec, yvec, interfaces.c36 + mi*slice, A.x, A.y);
            float c44_1 = BilinearInterpolation(xvec, yvec, interfaces.c44 + mi*slice, A.x, A.y);
            float c45_1 = BilinearInterpolation(xvec, yvec, interfaces.c45 + mi*slice, A.x, A.y);
            float c46_1 = BilinearInterpolation(xvec, yvec, interfaces.c46 + mi*slice, A.x, A.y);
            float c55_1 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c56_1 = BilinearInterpolation(xvec, yvec, interfaces.c56 + mi*slice, A.x, A.y);
            float c66_1 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float rho_1 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[0]  = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1]  = c11_0 + (c11_1 - c11_0)*dis_r; 
            var[2]  = c12_0 + (c12_1 - c12_0)*dis_r;
            var[3]  = c13_0 + (c13_1 - c13_0)*dis_r;
            var[4]  = c14_0 + (c14_1 - c14_0)*dis_r;
            var[5]  = c15_0 + (c15_1 - c15_0)*dis_r;
            var[6]  = c16_0 + (c16_1 - c16_0)*dis_r;
            var[7]  = c22_0 + (c22_1 - c22_0)*dis_r;
            var[8]  = c23_0 + (c23_1 - c23_0)*dis_r;
            var[9]  = c24_0 + (c24_1 - c24_0)*dis_r;
            var[10] = c25_0 + (c25_1 - c25_0)*dis_r;
            var[11] = c26_0 + (c26_1 - c26_0)*dis_r;
            var[12] = c33_0 + (c33_1 - c33_0)*dis_r;
            var[13] = c34_0 + (c34_1 - c34_0)*dis_r;
            var[14] = c35_0 + (c35_1 - c35_0)*dis_r;
            var[15] = c36_0 + (c36_1 - c36_0)*dis_r;
            var[16] = c44_0 + (c44_1 - c44_0)*dis_r;
            var[17] = c45_0 + (c45_1 - c45_0)*dis_r;
            var[18] = c46_0 + (c46_1 - c46_0)*dis_r;
            var[19] = c55_0 + (c55_1 - c55_0)*dis_r;
            var[20] = c56_0 + (c56_1 - c56_0)*dis_r;
            var[21] = c66_0 + (c66_1 - c66_0)*dis_r;
        }     
    break;

    case ELASTIC_TTI_BOND: /* 6. rho c11 c33 c55 c66 c13 azimuth dip */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            var[2] = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            var[3] = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            var[4] = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            var[5] = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            var[6] = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y);
            var[7] = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
        }
        else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float c11_0 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33_0 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c55_0 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c66_0 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13_0 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float rho_0 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float azi_0 = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y); // azimuth
            float dip_0 = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
            mi += 1;
            float c11_1 = BilinearInterpolation(xvec, yvec, interfaces.c11 + mi*slice, A.x, A.y);
            float c33_1 = BilinearInterpolation(xvec, yvec, interfaces.c33 + mi*slice, A.x, A.y);
            float c55_1 = BilinearInterpolation(xvec, yvec, interfaces.c55 + mi*slice, A.x, A.y);
            float c66_1 = BilinearInterpolation(xvec, yvec, interfaces.c66 + mi*slice, A.x, A.y);
            float c13_1 = BilinearInterpolation(xvec, yvec, interfaces.c13 + mi*slice, A.x, A.y);
            float rho_1 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            float azi_1 = BilinearInterpolation(xvec, yvec, interfaces.azimuth + mi*slice, A.x, A.y); // azimuth
            float dip_1 = BilinearInterpolation(xvec, yvec, interfaces.dip     + mi*slice, A.x, A.y);
            var[0] = rho_0 + (rho_1 - rho_0)*dis_r;
            var[1] = c11_0 + (c11_1 - c11_0)*dis_r;
            var[2] = c33_0 + (c33_1 - c33_0)*dis_r;
            var[3] = c55_0 + (c55_1 - c55_0)*dis_r;
            var[4] = c66_0 + (c66_1 - c66_0)*dis_r;
            var[5] = c13_0 + (c13_1 - c13_0)*dis_r;
            var[6] = azi_0 + (azi_1 - azi_0)*dis_r;
            var[7] = dip_0 + (dip_1 - dip_0)*dis_r;
        }   
    break;

    case ACOUSTIC_ISOTROPIC: /* 7. rho vp */
        if (mi == -1) {
            mi = findLastGreaterEqualIndex(elevation[0], elevation);
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
        } else if (mi == NI-1) {
            var[0] = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[1] = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
        } else {
            // special treatment: point is on the media grid, average of upper and lower media
            if (mi > 0 && isEqual(elevation[mi], A.z) && layer_indx.count(mi)) {
                mi--;
                dis_r = 1.0/2.0;
            }
            float vp0  = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            float rho0 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            mi+=1;
            float vp1  = BilinearInterpolation(xvec, yvec, interfaces.vp  + mi*slice, A.x, A.y);
            float rho1 = BilinearInterpolation(xvec, yvec, interfaces.rho + mi*slice, A.x, A.y);
            var[0] = rho0 + (rho1 - rho0) * dis_r;
            var[1] = vp0  + (vp1 - vp0 ) * dis_r;
        }
    break;

    default: // for self-check
        fprintf(stderr,"Error: Unknow media! (for code check, please contact Luqian Jiang)\n");
        fflush(stderr);
        exit(1);

    }     
} 

//- 0. assign the parameter directly (use the local values): one component 
void parametrization_grid_onecmp_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type,
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *var3d, int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    float slow_k = 1.0/(nz-1); // for print progress
    if (myid == 0)
        std::cout << "- discrete model by local values:\n\n";

    for (size_t k = 0; k < nz; ++k) {

        if (myid ==0) printProgress(slow_k*k);

        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                std::vector<float> var(1, 0.0);

                size_t indx =  i + j * siz_line + k * siz_slice;
                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }
                
                AssignGridMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ONE_COMPONENT, var, NGz, NL);
    
                var3d[indx] = var[0];     
            }
        }
    }

}

//- 1. assign the parameter directly (use the local values): isotropic, acoustic 
void parametrization_grid_ac_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *kappa,
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    float slow_k = 1.0/(nz-1); // for print progress
    if (myid == 0)
        std::cout << "- discrete model by local values:\n\n";
    for (size_t k = 0; k < nz; ++k) {

        if (myid == 0) printProgress(slow_k*k);

        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                std::vector<float> var(2, 0.0); // rho, vp

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignGridMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ACOUSTIC_ISOTROPIC, var, NGz, NL);

                /*calculate output kappa */
                float vp = var[1], rho = var[0];
                kappa[indx] = vp*vp*rho;
                rho3d[indx] = rho; 
            }
        }
    }
}


//- 2. assign the parameter directly (use the local values): elastic isotropic 
void parametrization_grid_el_iso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *rho3d, 
    float *lam3d,
    float *mu3d, 
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    float slow_k = 1.0/(nz-1); // for print progress
    if (myid == 0)
        std::cout << "- discrete model by local values:\n\n";

    for (size_t k = 0; k < nz; ++k) {

        if (myid==0) printProgress(slow_k*k);

        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                std::vector<float> var(3, 0.0); // rho, vp, vs

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                AssignGridMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);

                /*calculate output lambda and mu */
                float vp = var[1], vs = var[2], rho = var[0];
                mu3d[indx]  = vs*vs*rho;
                lam3d[indx] = vp*vp*rho - 2.0*mu3d[indx];
                rho3d[indx] = rho; 
                
            }
        }
    }
}


//- 3. Assign the parameter directly (use the local values): elastic vti
void parametrization_grid_el_vti_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;
    int media_type = interfaces.media_type;

    float slow_k = 1.0/(nz-1); // for print progress
    if (myid == 0)
        std::cout << "- discrete model by local values:\n\n";
    for (size_t k = 0; k < nz; ++k) {

        if (myid == 0) printProgress(slow_k*k);

        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y
                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                std::vector<float> var(6, 0.0); 
                AssignGridMediaPara2Point(i, j, k,
                        Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                        interfaces, media_type, var, NGz, NL);
                
                para2vti(var, media_type, 
                    c11[indx], c33[indx], c55[indx], c66[indx], c13[indx], rho[indx]);
                
            }
        }
    }
}

//- 4. Assign the parameter directly (use the local values): elastic tti
void parametrization_grid_el_aniso_loc(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, std::vector<int> &NGz,
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
    int myid)
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;
    int media_type = interfaces.media_type;

    float slow_k = 1.0/(nz-1); // for print progress
    if (myid == 0)
        std::cout << "- discrete model by local values:\n\n";

    for (size_t k = 0; k < nz; ++k) {
        if (myid == 0) printProgress(slow_k*k);

        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {

                size_t indx =  i + j * siz_line + k * siz_slice;

                size_t indx_z = indx;          // for vmap and curv: z
                size_t indx_x = i, indx_y = j; // for vmap and cart: x, y

                if (grid_type == GRID_CART) {
                    indx_z = k;                // for cart: z
                } else if (grid_type == GRID_CURV) {
                    indx_x = indx;
                    indx_y = indx;             // for curv: x, y
                }

                std::vector<float> var(22, 0.0); 

                AssignGridMediaPara2Point(i, j, k,
                            Point3(Gridx[indx_x], Gridy[indx_y], Gridz[indx_z]), 
                            interfaces, media_type, var, NGz, NL);

                para2tti(var, media_type, // return cij
                    c11[indx], c12[indx], c13[indx], c14[indx], c15[indx], c16[indx],
                    c22[indx], c23[indx], c24[indx], c25[indx], c26[indx],
                    c33[indx], c34[indx], c35[indx], c36[indx],
                    c44[indx], c45[indx], c46[indx], 
                    c55[indx], c56[indx], c66[indx], rho[indx]);
            }
        }
    }
}


//======================= averaging/equivalent medium method =========================

/* 
 * For half-grid point, mark the interface number.  
 */
int LayerNumberAtPoint(
    Point3 A, 
    int NL,
    std::vector<int> &NGz, 
    inter_t &interfaces) 
{
    std::vector<float> elevation(NL, -FLT_MAX);
    float MINX = interfaces.MINX;
    float MINY = interfaces.MINY;
    float   DX = interfaces.DX;
    float   DY = interfaces.DY;
    int     NX = interfaces.NX;
    int     NY = interfaces.NY;
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
        /* Get the elevation for the corresponding location */
        elevation[nl] = BilinearInterpolation(XVEC, YVEC, interfaces.elevation + ni*NX*NY, A.x, A.y);
        ni += NGz[nl];
    }

    /* Mark the layer number */
    int mi = findLastGreaterEqualIndex(A.z, elevation);

    /* If the elevation is above the surface, it is assigned by the first layer para. */
    if (mi == -1) 
        mi = findLastGreaterEqualIndex(elevation[0], elevation);

    return mi;
}

void MarkLayerNumber(
    int grid_type, 
    float *Hx, float *Hy, float *Hz,
    size_t nx, size_t ny, size_t nz,
    int NL, std::vector<int> &NGz,
    int *MaterNum,
    inter_t &interfaces) 
{
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;
    /* 
     * Mark the layer number at the half grid,
     *  for the grid media format, just 0 - NL-1 are marked.
     */
    if (grid_type == GRID_CART) {
        for (size_t k = 0; k < nz; k++) {
            for (size_t j = 0; j < ny; j++) {
                for (size_t i = 0; i < nx; i++) {
                    size_t indx =  i + j * siz_line + k * siz_slice; 
                    MaterNum[indx] = LayerNumberAtPoint(Point3(Hx[i], Hy[j], Hz[k]), 
                                                        NL, NGz, interfaces);
                }
            }
        }
    } else if (grid_type == GRID_VMAP) {
        for (size_t k = 0; k < nz; k++) {
            for (size_t j = 0; j < ny; j++) {
                for (size_t i = 0; i < nx; i++) {
                    size_t indx =  i + j * siz_line + k * siz_slice; 
                    MaterNum[indx] = LayerNumberAtPoint(Point3(Hx[i], Hy[j], Hz[indx]), 
                                                        NL, NGz, interfaces);
                }
            }
        }
    } else if (grid_type == GRID_CURV) {
        for (size_t i = 0; i < siz_volume; i++) {
            MaterNum[i] = LayerNumberAtPoint(Point3(Hx[i], Hy[i], Hz[i]), NL, NGz, interfaces);
        }
    }

}

//- 0.1 assign the parameter by volume harmonic averaging
//- one component
int parametrization_grid_onecmp_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *var3d, int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_onecmp_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, var3d, myid);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied. \n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic average method):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /* recalculate the material value of the point */
                    
                    float har_var = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg; 
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(1, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, ONE_COMPONENT, var, NGz, NL);
                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            har_var += (1.0/var[0]);
                        }
                    }

                    var3d[indx] = num_dis*1.0/har_var;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

//- 0.2 assign the parameter by volume arithmetic averaging
//- one component
int parametrization_grid_onecmp_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *var3d, int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_onecmp_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, var3d, myid);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied. \n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic integral average):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /* recalculate the material value of the point */
                    
                    float ari_var = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(1, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, ONE_COMPONENT, var, NGz, NL);
                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            ari_var += (var[0]);
                        }
                    }

                    var3d[indx] = ari_var/num_dis;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

//- 1.1 assign the parameter by volume arithmetic and harmonic averaging method
//- acoustic isotropic 
int parametrization_grid_ac_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *kappa,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_ac_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, rho3d, kappa, myid);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied. \n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /* recalculate the material value of the point */
                    float ari_rho   = 0.0;
                    float har_kappa = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, ACOUSTIC_ISOTROPIC, var, NGz, NL);
                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            float rho = var[0], vp = var[1];

                            ari_rho   += rho;
                            har_kappa += (1.0/(vp*vp*rho));
                        }
                    }

                    rho3d[indx] = ari_rho/num_dis;
                    kappa[indx] = 1.0*num_dis/har_kappa;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

//- 1.2 assign the parameter by volume arithmetic averaging method
//- acoustic isotropic 
int parametrization_grid_ac_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *kappa,
    int myid) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_ac_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, rho3d, kappa, myid);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /* recalculate the material value of the point */
                    float ari_rho   = 0.0;
                    float ari_kappa = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, ACOUSTIC_ISOTROPIC, var, NGz, NL);
                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            float rho = var[0], vp = var[1];
                            ari_rho   += rho;
                            ari_kappa += (vp*vp*rho);
                        }
                    }

                    rho3d[indx] = ari_rho/num_dis;
                    kappa[indx] = ari_kappa/num_dis;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}



//- 2.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic isotropic
// Moczo et al., 2002 
int parametrization_grid_el_iso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *lam3d,
    float *mu3d,
    int myid ) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_el_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, rho3d, lam3d, mu3d, myid);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium (volume arithmetic and harmonic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /* recalculate the material value of the point */
                    float ari_rho   = 0.0;
                    float har_kappa = 0.0;
                    float har_mu    = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);

                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            float rho = var[0], vp = var[1], vs = var[2];

                            float mu =  vs*vs*rho;
                            float lambda = vp*vp*rho - 2.0*mu;
    
                            ari_rho   += rho;
                            har_kappa += (1.0/(lambda + 2.0/3.0*mu));
                            har_mu    += (1.0/mu);
                        }
                    }

                    har_mu    = num_dis*1.0/har_mu;
                    har_kappa = num_dis*1.0/har_kappa;
                    ari_rho   = ari_rho/num_dis; 

                    lam3d[indx] = har_kappa - 2.0/3.0*har_mu;
                    mu3d[indx]  = har_mu;
                    rho3d[indx] = ari_rho;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

//- 2.1 assign the parameter by volume arithmetic averaging method
//- elastic isotropic
int parametrization_grid_el_iso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *rho3d,
    float *lam3d,
    float *mu3d,
    int myid ) 
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    
    /* assign the local value first.*/
    parametrization_grid_el_iso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, rho3d, lam3d, mu3d, myid);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);
    
    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                            siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float ari_lam = 0.0;
                    float ari_mu  = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);

                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(3, 0.0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, SubGrid[isg],
                            interfaces, ELASTIC_ISOTROPIC, var, NGz, NL);

                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            float vp = var[1], vs = var[2], rho = var[0];
                            float mu = vs*vs*rho;
                            float lambda = vp*vp*rho - 2.0*mu;
    
                            ari_rho += rho;
                            ari_lam += lambda;
                            ari_mu  += mu;
                        }
                    }

                    lam3d[indx] = ari_lam/(num_dis*1.0);
                    mu3d[indx]  = ari_mu /(num_dis*1.0);
                    rho3d[indx] = ari_rho/(num_dis*1.0);

                    if(SubGrid != nullptr) delete []SubGrid;

                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


//- 3.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic vti
int parametrization_grid_el_vti_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid)
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_vti_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
   
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";

    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float har_c11 = 0.0;
                    float har_c33 = 0.0;
                    float har_c55 = 0.0;
                    float har_c66 = 0.0;
                    float har_c13 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(6, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, SubGrid[isg],
                            interfaces, media_type, var, NGz, NL);

                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            // for every sub-point, transfer para to cij
                            float c11_p = 0.0, c33_p = 0.0;
                            float c55_p = 0.0, c66_p = 0.0;
                            float c13_p = 0.0, rho_p = 0.0;
                            para2vti(var, media_type, 
                                c11_p, c33_p, c55_p, c66_p, c13_p, rho_p);
    
                            har_c11 += (1.0/c11_p);
                            har_c13 += (1.0/c13_p);
                            har_c33 += (1.0/c33_p);
                            har_c55 += (1.0/c55_p);
                            har_c66 += (1.0/c66_p);
                            ari_rho += rho_p;
                        }
                    }

                    c11[indx] = 1.0*num_dis/har_c11;
                    c13[indx] = 1.0*num_dis/har_c13;
                    c33[indx] = 1.0*num_dis/har_c33;
                    c55[indx] = 1.0*num_dis/har_c55;
                    c66[indx] = 1.0*num_dis/har_c66;
                    rho[indx] = ari_rho/num_dis;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


//- 3.1 assign the parameter by volume arithmetic and harmonic averaging method
//- elastic vti
int parametrization_grid_el_vti_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces,
    float *c11,
    float *c33,
    float *c55,
    float *c66,
    float *c13,
    float *rho,
    int myid)
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_vti_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, c11, c33, c55, c66, c13, rho, myid);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    // recalculate the material value of the point
                    float ari_rho = 0.0;
                    float ari_c11 = 0.0;
                    float ari_c33 = 0.0;
                    float ari_c55 = 0.0;
                    float ari_c66 = 0.0;
                    float ari_c13 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(6, 0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, media_type, var, NGz, NL);
                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            // for every sub-point, transfer para to cij
                            float c11_p = 0.0, c33_p = 0.0;
                            float c55_p = 0.0, c66_p = 0.0;
                            float c13_p = 0.0, rho_p = 0.0;
                            para2vti(var, media_type, 
                                c11_p, c33_p, c55_p, c66_p, c13_p, rho_p);
    
                            ari_c11 += c11_p;
                            ari_c13 += c13_p;
                            ari_c33 += c33_p;
                            ari_c55 += c55_p;
                            ari_c66 += c66_p;
                            ari_rho += rho_p;
                        }
                    }

                    c11[indx] = ari_c11/num_dis;
                    c13[indx] = ari_c13/num_dis;
                    c33[indx] = ari_c33/num_dis;
                    c55[indx] = ari_c55/num_dis;
                    c66[indx] = ari_c66/num_dis;
                    rho[indx] = ari_rho/num_dis;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

//- 4.1 assign the parameter by volume arithmetic averaging method
//- elastic tti/anisotropic
int parametrization_grid_el_aniso_har(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
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
    int myid)
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_aniso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, c11, c12, c13, c14, c15, c16,
        c22, c23, c24, c25, c26, c33, c34, c35, c36, 
        c44, c45, c46, c55, c56, c66, rho, myid);
    
    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }

    /* For equivalent medium parameterization method */
    float *Hx = nullptr, *Hy = nullptr, *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    if (myid == 0)
        std::cout << "- equivalent medium parametrization (harmonic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {
        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /// recalculate the material value of the point
                    float ari_rho = 0.0;
                    float har_c11 = 0.0;
                    float har_c12 = 0.0;
                    float har_c13 = 0.0;
                    float har_c14 = 0.0;
                    float har_c15 = 0.0;
                    float har_c16 = 0.0;
                    float har_c22 = 0.0;
                    float har_c23 = 0.0;
                    float har_c24 = 0.0;
                    float har_c25 = 0.0;
                    float har_c26 = 0.0;
                    float har_c33 = 0.0;
                    float har_c34 = 0.0;
                    float har_c35 = 0.0;
                    float har_c36 = 0.0;
                    float har_c44 = 0.0;
                    float har_c45 = 0.0;
                    float har_c46 = 0.0;
                    float har_c55 = 0.0;
                    float har_c56 = 0.0;
                    float har_c66 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(22, 0.0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k, 
                            SubGrid[isg], interfaces, media_type, var, NGz, NL);

                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            // for every sub-point, transfer para to cij
                            float c11_p = 0.0, c12_p = 0.0, c13_p = 0.0;
                            float c14_p = 0.0, c15_p = 0.0, c16_p = 0.0;
                            float c22_p = 0.0, c23_p = 0.0, c24_p = 0.0;
                            float c25_p = 0.0, c26_p = 0.0, c33_p = 0.0;
                            float c34_p = 0.0, c35_p = 0.0, c36_p = 0.0;
                            float c44_p = 0.0, c45_p = 0.0, c46_p = 0.0;
                            float c55_p = 0.0, c56_p = 0.0, c66_p = 0.0;
                            float rho_p = 0.0;
                            para2tti(var, media_type, 
                                c11_p, c12_p, c13_p, c14_p, c15_p, c16_p,
                                c22_p, c23_p, c24_p, c25_p, c26_p,
                                c33_p, c34_p, c35_p, c36_p,
                                c44_p, c45_p, c46_p, 
                                c55_p, c56_p, c66_p, rho_p);
    
                            har_c11 += (1.0/c11_p);
                            har_c12 += (1.0/c12_p);
                            har_c13 += (1.0/c13_p);
                            har_c14 += (1.0/c14_p);
                            har_c15 += (1.0/c15_p);
                            har_c16 += (1.0/c16_p);
                            har_c22 += (1.0/c22_p);
                            har_c23 += (1.0/c23_p);
                            har_c24 += (1.0/c24_p);
                            har_c25 += (1.0/c25_p);
                            har_c26 += (1.0/c26_p);
                            har_c33 += (1.0/c33_p);
                            har_c34 += (1.0/c34_p);
                            har_c35 += (1.0/c35_p);
                            har_c36 += (1.0/c36_p);
                            har_c44 += (1.0/c44_p);
                            har_c45 += (1.0/c45_p);
                            har_c46 += (1.0/c46_p);
                            har_c55 += (1.0/c55_p);
                            har_c56 += (1.0/c56_p);
                            har_c66 += (1.0/c66_p);
                            ari_rho += rho_p;
                        }
                    }

                    c11[indx] = (1.0*num_dis)/har_c11;
                    c12[indx] = (1.0*num_dis)/har_c12;
                    c13[indx] = (1.0*num_dis)/har_c13;
                    c14[indx] = (1.0*num_dis)/har_c14;
                    c15[indx] = (1.0*num_dis)/har_c15;
                    c16[indx] = (1.0*num_dis)/har_c16;
                    c22[indx] = (1.0*num_dis)/har_c22;
                    c23[indx] = (1.0*num_dis)/har_c23;
                    c24[indx] = (1.0*num_dis)/har_c24;
                    c25[indx] = (1.0*num_dis)/har_c25;
                    c26[indx] = (1.0*num_dis)/har_c26;
                    c33[indx] = (1.0*num_dis)/har_c33;
                    c34[indx] = (1.0*num_dis)/har_c34;
                    c35[indx] = (1.0*num_dis)/har_c35;
                    c36[indx] = (1.0*num_dis)/har_c36;
                    c44[indx] = (1.0*num_dis)/har_c44;
                    c45[indx] = (1.0*num_dis)/har_c45;
                    c46[indx] = (1.0*num_dis)/har_c46;
                    c55[indx] = (1.0*num_dis)/har_c55;
                    c56[indx] = (1.0*num_dis)/har_c56;
                    c66[indx] = (1.0*num_dis)/har_c66;
                    
                    rho[indx] = ari_rho/num_dis;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}

//- 4.2 assign the parameter by volume arithmetic averaging method
//- elastic tti
int parametrization_grid_el_aniso_ari(
    size_t nx, 
    size_t ny, 
    size_t nz,
    const float *Gridx, 
    const float *Gridy, 
    const float *Gridz,
    int grid_type, 
    int NL, 
    std::vector <int> &NGz,
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
    int myid)
{

    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;
    size_t siz_volume = nz * siz_slice;
    int NI = interfaces.NI;
    int media_type = interfaces.media_type;

    /* assign the local value first.*/
    parametrization_grid_el_aniso_loc(nx, ny, nz, Gridx, Gridy, Gridz,
        grid_type, NL, NGz, interfaces, c11, c12, c13, c14, c15, c16,
        c22, c23, c24, c25, c26, c33, c34, c35, c36, 
        c44, c45, c46, c55, c56, c66, rho, myid);

    if (NL < 2) {
        fprintf(stdout, "There is %d layer, can not apply equivalent medium method, "\
            "the loc method is applied.\n", NL);        
        fflush(stdout);
        return 0;
    }
    
    /* For equivalent medium parameterization method */
    float *Hx = nullptr; 
    float *Hy = nullptr; 
    float *Hz = nullptr; 
    GenerateHalfGrid(nx, ny, nz, grid_type, Gridx, Gridy, Gridz, &Hx, &Hy, &Hz); 

    // mark the interface number
    int *MaterNum = new int[siz_volume];
    MarkLayerNumber(grid_type, Hx, Hy, Hz, nx, ny, nz, 
            NL, NGz, MaterNum, interfaces);

    /* Loop integer-grid points */ 
    float slow_k = 1.0/(ny-2);
    std::cout << "- equivalent medium parametrization (arithmetic averaging):\n\n";
    for (size_t j = 1; j < ny-1; j++) {

        if (myid == 0) printProgress(slow_k*j);

        for (size_t k = 1; k < nz-1; k++) {
            for (size_t i = 1; i < nx-1; i++) {
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

                    Mesh3 M = GenerateHalfMesh(grid_type, i, j, k, indx, 
                        siz_line, siz_slice, siz_volume, Hx, Hy, Hz);

                    /// recalculate the material value of the point
                    float ari_rho = 0.0;
                    float ari_c11 = 0.0;
                    float ari_c12 = 0.0;
                    float ari_c13 = 0.0;
                    float ari_c14 = 0.0;
                    float ari_c15 = 0.0;
                    float ari_c16 = 0.0;
                    float ari_c22 = 0.0;
                    float ari_c23 = 0.0;
                    float ari_c24 = 0.0;
                    float ari_c25 = 0.0;
                    float ari_c26 = 0.0;
                    float ari_c33 = 0.0;
                    float ari_c34 = 0.0;
                    float ari_c35 = 0.0;
                    float ari_c36 = 0.0;
                    float ari_c44 = 0.0;
                    float ari_c45 = 0.0;
                    float ari_c46 = 0.0;
                    float ari_c55 = 0.0;
                    float ari_c56 = 0.0;
                    float ari_c66 = 0.0;

                    Point3 *SubGrid = MeshSubdivide(M);
                    int nsg = (NG+1)*(NG+1)*(NG+1);
                    int num_dis = nsg;
                    for (int isg = 0; isg < nsg; isg++) {

                        std::vector<float> var(22, 0.0);

                        int isPointOnInter = AssignGridMediaPara2Point(i, j, k,
                            SubGrid[isg], interfaces, media_type, var, NGz, NL);

                        if (isPointOnInter == 1) {
                            num_dis--;
                        } else{
                            // for every sub-point, transfer para to cij
                            float c11_p = 0.0, c12_p = 0.0, c13_p = 0.0;
                            float c14_p = 0.0, c15_p = 0.0, c16_p = 0.0;
                            float c22_p = 0.0, c23_p = 0.0, c24_p = 0.0;
                            float c25_p = 0.0, c26_p = 0.0, c33_p = 0.0;
                            float c34_p = 0.0, c35_p = 0.0, c36_p = 0.0;
                            float c44_p = 0.0, c45_p = 0.0, c46_p = 0.0;
                            float c55_p = 0.0, c56_p = 0.0, c66_p = 0.0;
                            float rho_p = 0.0;
                            para2tti(var, media_type, 
                                c11_p, c12_p, c13_p, c14_p, c15_p, c16_p,
                                c22_p, c23_p, c24_p, c25_p, c26_p,
                                c33_p, c34_p, c35_p, c36_p,
                                c44_p, c45_p, c46_p, 
                                c55_p, c56_p, c66_p, rho_p);
    
                            ari_c11 += c11_p;
                            ari_c12 += c12_p;
                            ari_c13 += c13_p;
                            ari_c14 += c14_p;
                            ari_c15 += c15_p;
                            ari_c16 += c16_p;
                            ari_c22 += c22_p;
                            ari_c23 += c23_p;
                            ari_c24 += c24_p;
                            ari_c25 += c25_p;
                            ari_c26 += c26_p;
                            ari_c33 += c33_p;
                            ari_c34 += c34_p;
                            ari_c35 += c35_p;
                            ari_c36 += c36_p;
                            ari_c44 += c44_p;
                            ari_c45 += c45_p;
                            ari_c46 += c46_p;
                            ari_c55 += c55_p;
                            ari_c56 += c56_p;
                            ari_c66 += c66_p;
                            ari_rho += rho_p;
                        }
                    }

                    c11[indx] = ari_c11/num_dis;
                    c12[indx] = ari_c12/num_dis;
                    c13[indx] = ari_c13/num_dis;
                    c14[indx] = ari_c14/num_dis;
                    c15[indx] = ari_c15/num_dis;
                    c16[indx] = ari_c16/num_dis;
                    c22[indx] = ari_c22/num_dis;
                    c23[indx] = ari_c23/num_dis;
                    c24[indx] = ari_c24/num_dis;
                    c25[indx] = ari_c25/num_dis;
                    c26[indx] = ari_c26/num_dis;
                    c33[indx] = ari_c33/num_dis;
                    c34[indx] = ari_c34/num_dis;
                    c35[indx] = ari_c35/num_dis;
                    c36[indx] = ari_c36/num_dis;
                    c44[indx] = ari_c44/num_dis;
                    c45[indx] = ari_c45/num_dis;
                    c46[indx] = ari_c46/num_dis;
                    c55[indx] = ari_c55/num_dis;
                    c56[indx] = ari_c56/num_dis;
                    c66[indx] = ari_c66/num_dis;
                    
                    rho[indx] = ari_rho/num_dis;
 
                    if (SubGrid != nullptr) delete[] SubGrid;
                }

            }
        }
    }


    if (Hx != nullptr) delete[] Hx;
    if (Hy != nullptr) delete[] Hy;
    if (Hz != nullptr) delete[] Hz;
    if (MaterNum != nullptr) delete[] MaterNum;
    return 0;
}


