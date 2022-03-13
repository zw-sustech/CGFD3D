#ifndef _MEDIA_READ_INTERFACE_FILE_
#define _MEDIA_READ_INTERFACE_FILE_
#define MAX_BUF_LEN 1024

#include <iostream>
#include <vector>
#include "media_utility.hpp"

FILE *gfopen(const char *filename, const char *mode);

void read_interface_file(
    const char *interface_file,
    inter_t *interfaces);

/* 
 * Just read the grid data within the given
 *  [Xmin, Xmax]\times[Ymin, Ymax] domain.
 */
void read_grid_file(
    const char *grid_file,
    // the given grid
    float Xmin, float Xmax,
    float Ymin, float Ymax,
    int &NL,
    std::vector<int> &NGz, // how many z-grid in each layer
    inter_t *interfaces);

// check whether the elevation[ng[i]-1] == elevation[ng[i]] 
int checkGridData(int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces, 
    const char *grid_file); 

void read_bin_file(
    const char *bin_file,
    float *var,
    int dimx, 
    int dimy, 
    int dimz,
    int *bin_start, 
    int *bin_end, 
    int *bin_size, 
    size_t bin_line, 
    size_t bin_slice);

#endif