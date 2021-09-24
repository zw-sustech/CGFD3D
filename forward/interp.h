#ifndef _HEADER_interp
#define _HEADER_interp

float 
LagInterp_2d(float *x, float *y, float *z, int ni, int nj, float xi, float yi);


float 
LagInterp_3d(float *x, float *y, float *z, float *f,
             int ni, int nj, int nk, 
             float xi, float yi, float zi);

float 
LagrangeInt3d_WF( float *Don_var, float *Don_xcoor, float *Don_ycoor, float *Don_zcoor, 
                  int Don_nx, int Don_ny, int Don_nz,
                  int Donxcoor, int Donycoor, int Donzcoor,
                  float Tar_xcoor, float Tar_ycoor, float Tar_zcoor,
                  int InterOrder);

float 
LagInterp_Piecewise_1d(float *x, float *z, int ni, int order, float t_start, float dt, float ti);

#endif
