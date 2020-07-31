#include <math.h>
#include "common.h"

#ifdef ARCH_x86
void invert3x3(float m[][3]){
  float inv[3][3];
  float det;
  int i, j;

  inv[0][0] = m[1][1]*m[2][2] - m[2][1]*m[1][2];
  inv[0][1] = m[2][1]*m[0][2] - m[0][1]*m[2][2];
  inv[0][2] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
  inv[1][0] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
  inv[1][1] = m[0][0]*m[2][2] - m[2][0]*m[0][2];
  inv[1][2] = m[1][0]*m[0][2] - m[0][0]*m[1][2];
  inv[2][0] = m[1][0]*m[2][1] - m[1][1]*m[2][0];
  inv[2][1] = m[2][0]*m[0][1] - m[0][0]*m[2][1];
  inv[2][2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

  det = inv[0][0] * m[0][0] 
      + inv[0][1] * m[1][0] 
      + inv[0][2] * m[2][0];

  det = 1.0f / det;

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      m[i][j] = inv[i][j] * det;

  return ;
}

void matmul3x3(float A[][3], float B[][3], float C[][3]){
  
  int i, j, k;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
      C[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
    }

  return ;
}
#endif

void cross_product(float *A, float *B, float *C){
  C[0] = A[1] * B[2] - A[2] * B[1];
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
  return;
}

float dot_product(float *A, float *B){
  int i;
  float result = 0.0;
  for (i = 0; i < 3; i++)
    result += A[i] * B[i];

  return result;
}

float dist_point2plane(float x0[3], float x1[3], float x2[3], float x3[3]){

  float x12[3], x13[3], p[3];
  x12[0] = x2[0] - x1[0];
  x12[1] = x2[1] - x1[1];
  x12[2] = x2[2] - x1[2];
  x13[0] = x3[0] - x1[0];
  x13[1] = x3[1] - x1[1];
  x13[2] = x3[2] - x1[2];

  cross_product(x12, x13, p);
  float d = dot_product(p, x1);
  float L = (float)fabs( dot_product(p, x0) - d);
  L = L/sqrtf(dot_product(p, p));

  return L;
}

void bubble_sort(float a[], int index[], int n){
  int i, j;float tmp;

  for (i = 0; i < n; i++) index[i] = i;
  int tmpi;

  for (j = 0; j < n-1; j++)
    for (i = 0; i < n-1-j; i++){
      if(a[i] > a[i+1]){
        tmp = a[i];
        a[i] = a[i+1];
        a[i+1] = tmp;
        tmpi = index[i];
        index[i] = index[i+1];
        index[i+1] = tmpi;
      }
    }
}

void bubble_sort_int(int a[], int index[], int n){
  int i, j, tmp;

  for (i = 0; i < n; i++) index[i] = i;
  int tmpi;

  for (j = 0; j < n-1; j++)
    for (i = 0; i < n-1-j; i++){
      if(a[i] > a[i+1]){
        tmp = a[i];
        a[i] = a[i+1];
        a[i+1] = tmp;
        tmpi = index[i];
        index[i] = index[i+1];
        index[i+1] = tmpi;
      }
    }
}

