#ifndef FDLIB_MATH_H
#define FDLIB_MATH_H

void
fdlib_math_invert3x3(float m[][3]);

void
fdlib_math_matmul3x3(float A[][3], float B[][3], float C[][3]);

void
fdlib_math_cross_product(float *A, float *B, float *C);

float
fdlib_math_dot_product(float *A, float *B);

float
fdlib_math_dist_point2plane(float x0[3], float x1[3], float x2[3], float x3[3]);

void
fdlib_math_bubble_sort(float a[], int index[], int n);

void
fdlib_math_bubble_sort_int(int a[], int index[], int n);

int
fdlib_math_isPoint2InQuad(float px, float py, const float *vertx, const float *verty);

float
fdlib_math_rdinterp_2d(float x, float z, 
                  int num_points,
                  float *points_x, // x coord 
                  float *points_z, // z coord
                  float *points_v);
                  
int 
fdlib_math_isPoint2InTri(float px, float py, const float *vx, const float *vy);                  
 
int
fdlib_math_Point2InWhichTri(float px, float py, const float *vertx, const float *verty);      


float
fdlib_math_computeZfromTri(float x, float y, 
													 float x1, float y1, float z1,
													 float x2, float y2, float z2,
													 float x3, float y3, float z3);


#endif
