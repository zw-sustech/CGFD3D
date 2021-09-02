#include <iostream>
#include "media_layer2model.hpp"
#include "media_grid2model.hpp"
using namespace std;

int main()
{

    size_t nx = 501, ny = 501, nz = 300;
    float  x0 = -2500.0, y0 = 0.0, z0 = -3000.0;
    float  dx = 10.0, dy = 10.0, dz = 10.0;

    // x dimention varies first
    size_t siz_line   = nx; 
    size_t siz_slice  = nx * ny; 
    size_t siz_volume = nx * ny * nz;

    float *Gridx = (float*) malloc(siz_volume*sizeof(float));
    float *Gridy = (float*) malloc(siz_volume*sizeof(float));
    float *Gridz = (float*) malloc(siz_volume*sizeof(float));
    float *lam3d = (float*) malloc(siz_volume*sizeof(float));
    float *mu3d  = (float*) malloc(siz_volume*sizeof(float));
    float *rho3d = (float*) malloc(siz_volume*sizeof(float));

    for (size_t k = 0; k < nz; k++) {
        for (size_t j = 0; j < ny; j++) {
            for (size_t i = 0; i < nx; i++) {
                size_t indx =  i + j * siz_line + k * siz_slice;
                Gridx[indx] = x0 + i * dx;
                Gridz[indx] = z0 + k * dz;
                Gridy[indx] = y0 + j * dy;
            }
        }
    }


 //   cout << Gridx[0] << " " << Gridx[siz_volume-1] << endl;
 //   cout << Gridy[0] << " " << Gridy[siz_volume-1] << endl;

    media_layer2model_curv_el_iso( lam3d, mu3d, rho3d,
        Gridx, Gridy, Gridz, nx,ny, nz, 
        "can4_rho.md3lay",
        "can4_vp.md3lay",
        "can4_vs.md3lay", 
        "har"); 


    // import
    FILE *fp;
    fp = gfopen("rho.dat","w");
    for (size_t i = 0; i < siz_volume; i++)
        fwrite(&rho3d[i],sizeof(float),1,fp);  
    fclose(fp);

    fp = gfopen("lam.dat","w");
    for (size_t i = 0; i < siz_volume; i++)
        fwrite(&lam3d[i],sizeof(float),1,fp);  
    fclose(fp);

    fp = gfopen("mu.dat","w");
    for (size_t i = 0; i < siz_volume; i++)
        fwrite(&mu3d[i],sizeof(float),1,fp);  
    fclose(fp);

    free(lam3d); free(mu3d); free(rho3d);
    free(Gridx); free(Gridy); free(Gridz);

    return 0;
}
