#include <math.h>

#include "fdlib_math.h"

void
fdlib_math_invert3x3(float m[][3])
{
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

void
fdlib_math_matmul3x3(float A[][3], float B[][3], float C[][3])
{
  int i, j, k;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++){
      C[i][j] = 0.0;
      for (k = 0; k < 3; k++)
        C[i][j] += A[i][k] * B[k][j];
    }

  return ;
}

void
fdlib_math_cross_product(float *A, float *B, float *C)
{
  C[0] = A[1] * B[2] - A[2] * B[1];
  C[1] = A[2] * B[0] - A[0] * B[2];
  C[2] = A[0] * B[1] - A[1] * B[0];
  return;
}

float
fdlib_math_dot_product(float *A, float *B)
{
  int i;
  float result = 0.0;
  for (i = 0; i < 3; i++)
    result += A[i] * B[i];

  return result;
}

float
fdlib_math_dist_point2plane(float x0[3], float x1[3], float x2[3], float x3[3])
{

  float x12[3], x13[3], p[3];
  x12[0] = x2[0] - x1[0];
  x12[1] = x2[1] - x1[1];
  x12[2] = x2[2] - x1[2];
  x13[0] = x3[0] - x1[0];
  x13[1] = x3[1] - x1[1];
  x13[2] = x3[2] - x1[2];

  fdlib_math_cross_product(x12, x13, p);
  float d = fdlib_math_dot_product(p, x1);
  float L = (float)fabs( fdlib_math_dot_product(p, x0) - d);
  L = L/sqrtf(fdlib_math_dot_product(p, p));

  return L;
}

void
fdlib_math_bubble_sort(float a[], int index[], int n)
{
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

void
fdlib_math_bubble_sort_int(int a[], int index[], int n)
{
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

/*
 * Input: vertx, verty are the FOUR vertexes of the quadrilateral
 *
 *    ↑ +z    
 *    |       
 *            2----3 
 *            |    |
 *            |    | 
 *            0----1
 *
 * The quadrilateral must be convex!
 *
 * (anti-clockwise + Judge whether the v->p is on the left side of the edge)
 */
int
fdlib_math_isPoint2InQuad(float px, float py, const float *vertx, const float *verty)
{
    float vx[4], vy[4]; // for colockwise
    vx[0] = vertx[0]; vy[0] = verty[0];
    vx[1] = vertx[1]; vy[1] = verty[1];
    vx[2] = vertx[3]; vy[2] = verty[3];
    vx[3] = vertx[2]; vy[3] = verty[2];

    for (int i = 0; i < 4; i++) {
        // vector: edge
        float vecte_x = i==3?(vx[0]-vx[3]):(vx[i+1]-vx[i]);
        float vecte_y = i==3?(vy[0]-vy[3]):(vy[i+1]-vy[i]);
        // vector: vertex to point
        float vectp_x = px - vx[i];
        float vectp_y = py - vy[i];

        float sign = vecte_x*vectp_y - vecte_y*vectp_x;
        sign /= (sqrt(vecte_x*vecte_x + vecte_y*vecte_y));
        // normalization
        
        if (sign < -1e-6)
            return 0;
    }

    return 1;
}

/* 
 * interp value in 2D grid using inverse distance interp
 */

float
fdlib_math_rdinterp_2d(float x, float z, 
                  int num_points,
                  float *points_x, // x coord 
                  float *points_z, // z coord
                  float *points_v)
{
  float val = 0.0; // return value

  float weight[num_points];
  float total_weight = 0.0 ;

  // cal weight
  int at_point_indx = -1;
  for (int i=0; i<num_points; i++)
  {
    float dist = sqrtf ((x - points_x[i]) * (x - points_x[i])
        + (z - points_z[i]) * (z - points_z[i])
        );
    if (dist < 1e-9) {
      at_point_indx = i;
    } else {
      weight[i]   = 1.0 / dist;
      total_weight += weight[i];
    }
  }
  // if at a point
  if (at_point_indx > 0) {
    total_weight = 1.0;
    // other weight 0
    for (int i=0; i<num_points; i++) {
      weight[i] = 0.0;
    }
    // point weight 1
    weight[at_point_indx] = 1.0;
  }

  // interp

  for (int i=0; i<num_points; i++)
  {
    weight[i] *= 1.0 / total_weight ;

    val += weight[i] * points_v[i];  

    //fprintf(stdout,"---- i=%d,weight=%f,points_v=%f\n",
    //    i,weight[i],points_v[i]);
  }

  return val;
}
