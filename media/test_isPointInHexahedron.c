#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "isPointInHexahedron.h"

int main()
{
    /*
     * Input: vx, vy, vz are the EIGHT vertexes of the hexahedron 
     *
     *    ↑ +z       4----6
     *    |         /|   /|
     *             / 0--/-2
     *            5----7 /
     *            |/   |/
     *            1----3
     */
    float x0 = 0.0, dx = 2.0;
    float y0 = 0.0, dy = 1.0;
    float z0 = 0.0, dz = 10.0;

    float vx[8], vy[8], vz[8];

    for (size_t k = 0; k < 2; k++) {
        for (size_t j = 0; j < 2; j++) {
            for (size_t i = 0; i < 2; i++) {
                size_t indx =  i + j * 2 + k * 4;
                vx[indx] = x0 + i * dx;
                vz[indx] = z0 + k * dz;
                vy[indx] = y0 + j * dy;
            }
        }
    }

//    for (int i = 0; i < 8; i++) {
//        printf("%f  %f %f\n", vx[i], vy[i], vz[i]);
//    }

// USAGE !!
    if ( isPointInHexahedron(2.0,0.5,5, vx, vy, vz) ) 
        printf("\nInside\n");
    else 
        printf("\nOutside\n");

    return 0;
}