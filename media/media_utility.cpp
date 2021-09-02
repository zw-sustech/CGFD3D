#include <iostream>
#include <vector>
#include <cfloat>
#include <cmath>
#include <string.h>
#include "media_utility.hpp"

bool isEqual(float a, float b) {
    return abs(a-b) < FLT_EPSILON;
}

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
    /* 
     * half grid point is inside the grid, 
     * so the number of grid point in one direction is np-1.
     */
    size_t siz_line = nx;
    size_t siz_slice = ny * siz_line;

    for (size_t k = 0; k < nz-1; k++) {
        for (size_t j = 0; j < ny-1; j++) {
            for (size_t i = 0; i < nx-1; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice;

                /* clockwise: conducive to debug */
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


/* 
 * How many different values in the vector, 
 *  NI is the upper limitation of the v[i].
 */
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

/* 
 * Find the index of the nearest value (x[i] <= value), 
 *  the x vector can be arbitrary
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

/* 
 * Find the index of the nearest value and no less than the value (value <= x[i]), 
 *    Just for the ordered array x is from largest to smallest.
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

        /* 
         * If the value is equal to the value in x, return the first index of the value,
         * If the value is not in the array x, return the last index which x[index] > value.
         */ 
        newDist = x[i] - value;   
        if (newDist >= 0 && newDist <= dist) {
            if (newDist <= 1e-6) 
                return i;
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

            vq =  v[indx]      * (x[ix+1] - xq) * (y[iy+1] - yq) 
                + v[indx+1]    * (xq -   x[ix]) * (y[iy+1] - yq)
                + v[indx+NX]   * (x[ix+1] - xq) * (yq - y[iy]  )
                + v[indx+1+NX] * (xq -   x[ix]) * (yq - y[iy]  );

            vq /= area;
        }
    }
    else
        vq = FLT_MAX;

    return vq;
}
