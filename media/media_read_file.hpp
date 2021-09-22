#ifndef __MEDIA_READ_INTERFACE_FILE__
#define __MEDIA_READ_INTERFACE_FILE__
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

#endif