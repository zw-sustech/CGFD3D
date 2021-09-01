#ifndef __MEDIA_READ_INTERFACE_FILE__
#define __MEDIA_READ_INTERFACE_FILE__
#define MAX_BUF_LEN 1024

#include <iostream>
#include <vector>
#include "media_utility.hpp"

FILE *gfopen(const char *filename, const char *mode);

void read_interface_file(
    const char *interface_file,
    bool first_read, // Is it the first para to be read
    inter_t *interfaces,
    float **var, float **var_grad, float **var_pow);

void read_grid_file(
    const char *grid_file,
    // the calculation grid range
    float Xmin, float Xmax,
    float Ymin, float Ymax,
    int &NL,
    std::vector<int> &NGz, // how many z-grid in each layer
    inter_t **interfaces);

#endif