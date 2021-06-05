//#include <bits/stdc++.h>
#include <iostream>
#include <vector>
#include <cfloat>
#include "media_interpolation.hpp"
//using namespace std;

/* find the index of the nearest value (x[i] <= value), 
   the x vector can be arbitrary
*/
int findNearestNeighborIndex(
    float value, 
    std::vector<float> &x)
{
    float dist, newDist;
    int idx = -1, nx = x.size();
    // expolation: INF(idx = -1) or

    dist = FLT_MAX;
    for(int i = 0; i < nx; i++) {
        if (x[i] == -FLT_MAX) continue;
        newDist = value - x[i];
        if (newDist >= 0 && newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

/* find the index of the nearest value and no less than the value (x[i] >= value), 
   the x vector can be arbitrary 
*/
int findNearestGreaterIndex(
    float value, 
    std::vector<float> &x)
{
    float dist, newDist;
    int idx;

    idx = -1;    // expolation: INF(idx = -1) or
    dist = FLT_MAX;
    for(size_t i = 0; i < x.size(); i++) {
        if (x[i] == -FLT_MAX) continue;
        newDist = x[i] - value;
        if (newDist >= 0 && newDist < dist) {
            dist = newDist;
            idx = i;
        }
    }

    return idx;
}

/* 
 * BilinearInterpolate to the (xq, vq) point 
 * It's just for the problum: (x,y) is a grid and the x, y vector are increment 
 * If x, y are arbitrary, we need to re-detrimine points for interpolation
 */
float BilinearInterpolation(
    std::vector<float> &x, 
    std::vector<float> &y, 
    float *v,
    float xq,
    float yq )
{    

    float vq;
    int NX = x.size();
    int NY = y.size();

    int ix = -1; //findNearestNeighborIndex(xq, x);
    int iy = -1; //findNearestNeighborIndex(yq, y);
    ix = (xq-x[0])/(x[1]-x[0]);  // just for the given uniform interface mash.
    iy = (yq-y[0])/(y[1]-y[0]);

    if (ix >= 0 && iy >= 0) {
        int indx = iy*NX + ix;
        if (ix == NX-1 || iy == NY-1) {
            vq = v[indx];
        }
        else {
            float area = (x[ix+1] - x[ix]) * (y[iy+1] - y[iy]); 

            vq =  *(v+indx)        * (x[ix+1] - xq) * (y[iy+1] - yq) 
                + *(v+(indx+1   )) * (xq -   x[ix]) * (y[iy+1] - yq)
                + *(v+(indx+NX  )) * (x[ix+1] - xq) * (yq - y[iy]  )
                + *(v+(indx+1+NX)) * (xq -   x[ix]) * (yq - y[iy]  );
            vq /= area;
        }
    }
    else
        vq = FLT_MAX;

    return vq;
}


/* 
 * BilinearInterpolate to the vector xq, yq.
 */
//void
//interp2(
//    vector<float> &x, 
//    vector<float> &y, 
//    float **v,
//    vector<float> &xq,
//    vector<float> &yq,
//    float **vq )
//{
//    int len_xq = xq.size();
//    int len_yq = yq.size();
//    for (int i = 0; i < len_xq; i++) {
//        for (int j = 0; j < len_yq; j++) {
//            vq[i][j] = BilinearInterpolation(x, y, v, xq[i], yq[j]);
//        }
//    }
//}
//

/*
// Lagrange interpolation
void interp1(float *x, int x_len, float *v, float *xq, int xq_len, float *vq)
{
    float dx, dy, *slope, *intercept;
    int i, indiceEnVector;

    slope     = (float*) malloc( x_len*sizeof(float) );
    intercept = (float*) malloc( x_len*sizeof(float) );

    for (i = 0; i < x_len; i++) {
        if (i < x_len-1) {
            dx = x[i+1] - x[i];
            dy = v[i+1] - v[i];
            slope[i] = dy / dx;
            intercept[i] = v[i] - x[i] * slope[i];
        } else {
            slope[i] = slope[i-1];
            intercept[i] = intercept[i-1];
        }
    }

    for (i = 0; i < xq_len; i++) {
        indiceEnVector = findNearestNeighborIndex(xq[i], x, x_len);
        //indiceEnVector = (int) xq[i] ;  // For our problem
        if (indiceEnVector != -1)
            vq[i] = slope[indiceEnVector] * xq[i] + intercept[indiceEnVector];
        else
            vq[i] = DBL_MAX;
    }

    free(slope);
    free(intercept);
}
*/