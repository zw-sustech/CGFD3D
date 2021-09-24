#ifndef ISPOINTINHEXAHEDRON_H
#define ISPOINTINHEXAHEDRON_H

#include <stdbool.h>
// for C code call
bool isPointInHexahedron(float px, float py, float pz,
                         float *vx, float *vy, float *vz);

#endif
