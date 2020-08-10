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
#endif
