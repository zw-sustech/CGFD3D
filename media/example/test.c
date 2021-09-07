#include <stdio.h>
#include <stdlib.h>
#include "media_discrete_model.h"
int main()
{

    size_t nx = 501, ny = 501, nz = 450;
    float  x0 = -2000.0, y0 = 0.0, z0 = -4500.0;
    float  dx = 10.0, dy = 10.0, dz = 10.0;
    int grid_type = MEDIA_USE_CURV; // just for test

    // x dimention varies first
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    float *x1d = (float*) malloc(nx*sizeof(float));
    float *y1d = (float*) malloc(ny*sizeof(float));
    float *z1d = (float*) malloc(nz*sizeof(float));
    float *x3d = (float*) malloc(siz_volume*sizeof(float));
    float *y3d = (float*) malloc(siz_volume*sizeof(float));
    float *z3d = (float*) malloc(siz_volume*sizeof(float));
    float *lam3d = (float*) malloc(siz_volume*sizeof(float));
    float *mu3d  = (float*) malloc(siz_volume*sizeof(float));
    float *rho3d = (float*) malloc(siz_volume*sizeof(float));

    for (size_t k = 0; k < nz; k++) {
        z1d[k] = z0 + k*dz;
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice;
                x3d[indx] = x0 + i * dx;
                z3d[indx] = z0 + k * dz;
                y3d[indx] = y0 + j * dy;
            }
        }
    }
    for (size_t j = 0; j < ny; j++) {
        y1d[j] = y0 + j * dy;
    }
    for (size_t i = 0; i < nx; i++) {
        x1d[i] = x0 + i * dx;
    }


 //   cout << Gridx[0] << " " << Gridx[siz_volume-1] << endl;
 //   cout << Gridy[0] << " " << Gridy[siz_volume-1] << endl;

    if (grid_type == 1)
        media_layer2model_el_iso( lam3d, mu3d, rho3d,
            x1d, y1d, z1d, nx,ny, nz, MEDIA_USE_CART, 
            "can4_rho.md3lay",
            "can4_vp.md3lay",
            "can4_vs.md3lay", 
            "har"); 
    else if (grid_type == 2)
        media_layer2model_el_iso( lam3d, mu3d, rho3d,
            x1d, y1d, z3d, nx, ny, nz, MEDIA_USE_VMAP,
            "can4_rho.md3lay",
            "can4_vp.md3lay",
            "can4_vs.md3lay", 
            "har"); 
    else
        media_layer2model_el_iso( lam3d, mu3d, rho3d,
            x3d, y3d, z3d, nx,ny, nz, MEDIA_USE_CURV,
            "can4_rho.md3lay",
            "can4_vp.md3lay",
            "can4_vs.md3lay", 
            "har"); 


    // import
    FILE *fp;
    fp = fopen("rho.dat","w");
    for (size_t i = 0; i < siz_volume; i++)
        fwrite(&rho3d[i],sizeof(float),1,fp);  
    fclose(fp);

    fp = fopen("lam.dat","w");
    for (size_t i = 0; i < siz_volume; i++)
        fwrite(&lam3d[i],sizeof(float),1,fp);  
    fclose(fp);

    fp = fopen("mu.dat","w");
    for (size_t i = 0; i < siz_volume; i++)
        fwrite(&mu3d[i],sizeof(float),1,fp);  
    fclose(fp);

    free(lam3d); free(mu3d); free(rho3d);
    free(x1d); free(y1d); free(z1d);
    free(x3d); free(y3d); free(z3d);

    return 0;
}
