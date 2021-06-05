#ifndef __INTERPOLATION_H__
#define __INTERPOLATION_H__

#include <iostream>
#include <vector>
//using namespace std;

int findNearestNeighborIndex(
    float value, std::vector<float> &x);

int findNearestGreaterIndex(
    float value, std::vector<float> &x);

float BilinearInterpolation(
    std::vector<float> &x, 
    std::vector<float> &y, 
    float *v,
    float xq,
    float yq );

//void interp2(
//    vector<float> &x, 
//    vector<float> &y, 
//    float **v,
//    vector<float> &xq,
//    vector<float> &yq,
//    float **vq );



#endif