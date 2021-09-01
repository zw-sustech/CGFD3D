#include <iostream>
#include "media_read_file.hpp"

FILE *gfopen(const char *filename, const char *mode)
{
    FILE *fp;
    if ((fp = fopen(filename,mode)) == NULL) {
        fprintf(stderr, "Cannot open %s, " \
            "please check your file path and run-directory.\n",filename);
        exit(1);
    }
    return fp;
}


void read_interface_file(
    const char *interface_file,
    bool first_read, // Is it the first para to be read
    inter_t *interfaces,
    float **var, float **var_grad, float **var_pow)   // value need to read
{
    char line[MAX_BUF_LEN];
    FILE *file = gfopen(interface_file, "r");
    FILE *tmp_file = tmpfile();

    int   NI, NX, NY;
    float DX, DY, MINX, MINY;

    /* Read every line in the file. */
    while(fgets(line, MAX_BUF_LEN, file) != NULL)
    {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        fputs(line,tmp_file);
    } 

    /* Set the file pointer at the beginning of the stream */
    rewind(tmp_file);

    /* Read the temporary data file, Can not use the feof() */
    while(feof(tmp_file) != EOF)
    {
        if (fscanf(tmp_file, "%d", &NI) < 1) {
            fprintf(stderr,"Error: The INTERFACES FILE is wrong, please give a number of layers!\n");
            fflush(stderr);
            exit(1);
        }

        if (!first_read && NI != interfaces->NI) {
            fprintf(stderr, "Error: NI of %s is not equal to the existing one!\n", interface_file);
            fflush(stderr);
            exit(1);
        }

        if (NI < 1) {
            fprintf(stderr, "Error: No enough interfaes (minimum is 1)!\n");
            fflush(stderr);
            exit(1);
        }

        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"Error: The INTERFACES FILE is wrong, " \
                "please check the given interfaces mesh!\n");
            fflush(stderr);
            exit(1);
        }

        if (!first_read && ( NX != interfaces->NX || NY != interfaces->NY ||
            !isEqual(MINX, interfaces->MINX) || !isEqual(MINY, interfaces->MINY) ||
            !isEqual(DX, interfaces->DX) || !isEqual(DY, interfaces->DY) ) )
        {
            fprintf(stderr, "Error: the interface mesh of %s is not equal to the existing one!\n", interface_file);
            fflush(stderr);
            exit(1);
        }
        
        if (NX < 2 || NY < 2) {
            fprintf(stderr, "Error: No enough point (NX >= 2, NY >= 2)!");
            fflush(stderr);   
            exit(1);         
        }

        if (first_read) {
            interfaces->NI   = NI;
            interfaces->NX   = NX;
            interfaces->NY   = NY;
            interfaces->DX   = DX;
            interfaces->DY   = DY;
            interfaces->MINX = MINX;
            interfaces->MINY = MINY;
        }
        

        /* volume and slice of interfaces info */
        size_t inter_line   =  NX;
        size_t inter_slice  =  NX*NY;
        size_t inter_volume =  NI*inter_slice;
        if (first_read) {
            interfaces->elevation = new float[inter_volume];
        }
        *var       = new float[inter_volume];
        *var_grad  = new float[inter_volume];
        *var_pow   = new float[inter_volume];

        for (size_t ni = 0; ni < NI; ni++) {
            for (size_t j = 0; j < NY; j++) {
                for (size_t i = 0; i < NX; i++) {
                    size_t indx = i + j * inter_line + ni*inter_slice;
                    int num_read = 0;
                    if (first_read) {
                        num_read = fscanf(tmp_file, "%f %f %f %f",
                                    &(interfaces->elevation[indx]), 
                                    (*var)      +indx,
                                    (*var_grad) +indx,
                                    (*var_pow)  +indx );
                    } else {
                        float elevation = 0.0;
                        num_read = fscanf(tmp_file, "%f %f %f %f",
                                    &elevation, 
                                    &((*var)[indx]),
                                    &((*var_grad)[indx]),
                                    &((*var_pow)[indx]) );
                        if (!isEqual(elevation, interfaces->elevation[indx])) {
                            fprintf(stderr, "Error: the elevation of %s is not equal to the existing one!\n",interface_file);
                            fflush(stderr);
                            exit(1);
                        }

                    }
                    if (num_read < 4) {
                        fprintf(stderr,"Error: The INTERFACES FILE is wrong, " \
                            "please check the interfaces data of #%li layer!\n", ni);
                        fflush(stderr);
                        exit(1);
                    }
                    
                }
            }
        }

        // end of the info
        break;
    } 
    fclose(file);   
    fclose(tmp_file);   
}


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
    inter_t **interfaces)
{
    FILE *file = gfopen(grid_file, "r");
    FILE *tmp_file = tmpfile();

    char  line[MAX_BUF_LEN];
    int   NX = 0, NY = 0, NI = 0;
    float DX = 0, DY = 0;
    float MINX = 0, MINY = 0;

    /* Read every line in the file. */
    while(fgets(line, MAX_BUF_LEN, file) != NULL)
    {
        if (line[0] == '#' || line[0] == '\n')
            continue;
        fputs(line,tmp_file);
    } 

    /* Set the file pointer at the beginning of the stream*/
    rewind(tmp_file);

    /* Read the temporary data file, Can not use the feof() */
    while(feof(tmp_file) != EOF)
    {
        if (fscanf(tmp_file, "%d", &NL) < 1) {
            fprintf(stderr,"Error: The MEDIA GRID FILE is wrong, " \
                "please give a number of layers! \n");
            fflush(stderr);
            exit(1);
        }

        if (NL < 1) {
            fprintf(stderr, "Error: No enough layers (minimum is 1)!\n");
            fflush(stderr);
            exit(1);
        }

        for (int i = 0; i < NL; i++) {
            int ng_i = 0;
            if (fscanf(tmp_file, "%d",&ng_i) < 1) {
                fprintf(stderr,"Error: The MEDIA GRID FILE is wrong, " \
                    "please give the number of grids in the %d-th layer! \n", i);
                fflush(stderr);
                exit(1);
            }
            NI += ng_i;
            NGz.push_back(ng_i);
        }

        /* malloc interface */
        (*interfaces) = new inter_t [NI];


        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"Error: The INTERFACES FILE is wrong, " \
                "please check the given interfaces mesh! \n");
            fflush(stderr);
            exit(1);
        }

        if (NX < 2 || NY < 2) {
            fprintf(stderr,"Error: No enough point (NX >= 2, NY >= 2)! \n");
            fflush(stderr); 
            exit(1);           
        }

        /* the range need to read */
        size_t ix0 = (Xmin-MINX)/DX;
        size_t iy0 = (Ymin-MINY)/DY;
        size_t ix1 = (Xmax-MINX)/DX+1;
        size_t iy1 = (Ymax-MINY)/DY+1;
        size_t nx  = ix1-ix0+1;
        size_t ny  = iy1-iy0+1;
        float minx = ix0 * DX + MINX;
        float miny = iy0 * DY + MINY;


        /* The given media domain must bigger than the calculation domain */
        if (ix0 < 0 || ix1 >= NX ) {
            fprintf(stderr,"Error: The given media range is smaller than "\
                "the calculation grid range in x-direction! \n");
            fflush(stderr); 
            exit(1); 
        }
        if (iy0 < 0 || iy0 >= NY) {
            fprintf(stderr, "Error: The given media range is smaller than "\
                "the calculation grid range in y-direction! \n");
            fflush(stderr); 
            exit(1);             
        }

        for (int ni = 0; ni < NI; ni++) {

            (*interfaces)[ni].NX    = nx;
            (*interfaces)[ni].NY    = ny;
            (*interfaces)[ni].DX    = DX;
            (*interfaces)[ni].DY    = DY;
            (*interfaces)[ni].MINX  = minx;
            (*interfaces)[ni].MINY  = miny;
                  
            (*interfaces)[ni].elevation = new float [nx*ny];
            (*interfaces)[ni].vp        = new float [nx*ny];
            (*interfaces)[ni].vs        = new float [nx*ny];
            (*interfaces)[ni].rho       = new float [nx*ny];

            for (size_t iy = 0; iy < NY; iy++) {
                
                for (size_t ix = 0; ix < NX; ix++) {

                    if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {

                        size_t indx = (ix-ix0) + nx * (iy-iy0);

                        int u = fscanf(tmp_file, "%f %f %f %f",
                                       &(*interfaces)[ni].elevation[indx],  // elevation
                                       &(*interfaces)[ni].vp[indx],  
                                       &(*interfaces)[ni].vs[indx],  
                                       &(*interfaces)[ni].rho[indx] );
    
                        if (u < 4) {
                            fprintf(stderr,"Error: The INTERFACES FILE is wrong, " \
                                "please check the grid data of #%d z-grid! \n", ni);
                            fflush(stderr);
                            exit(1);
                        }
                    } else { // read the data no used.
                        float tmp[4];
                        int u = fscanf(tmp_file, "%f %f %f %f", tmp, tmp+1, tmp+2, tmp+3);
                        if (u < 4) {
                            fprintf(stderr,"Error: The INTERFACES FILE is wrong, " \
                                "please check the grid data of #%d z-grid! \n", ni);
                            fflush(stderr);
                            exit(1);
                        }
                    }
                        
                }
            
            }
        }

        // end of the info
        break;
    } 
    fclose(file);   
    fclose(tmp_file);    
}
