#ifndef __READ_INTERFACE_FILE__
#define __READ_INTERFACE_FILE__

#include "media_layer2model.hpp"
#define MAX_BUF_LEN 1024

FILE *gfopen(const char *filename, const char *mode);

void read_interface_file(
    const char *interface_file,
    int &NI,
    Interfaces **interfaces);

#endif