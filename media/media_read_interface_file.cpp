#include <iostream>
#include "media_read_interface_file.hpp"

FILE *gfopen(const char *filename, const char *mode)
{
    FILE *fp;
    if ((fp = fopen(filename,mode)) == NULL) {
        fprintf(stderr, "\033[47;31mCannot open %s, " \
            "please check your file path and run-directory.\033[0m\n",filename);
        exit(1);
    }
    return fp;
}


void read_interface_file(
    const char *interface_file,
    int &NI,
    Interfaces **interfaces)
{
    char line[MAX_BUF_LEN];

    FILE *file = gfopen(interface_file, "r");
    FILE *tmp_file = tmpfile();

    int NX, NY;
    float DX, DY, MINX, MINY;

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
        if (fscanf(tmp_file, "%d", &NI) < 1) {
            fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                "please give a number of layers!\033[0m\n");
            fflush(stderr);
            exit(1);
        }

        if (NI < 1) {
            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: No enough interfaes (minimum is 1)!");
            fflush(stderr);
            exit(1);
        }

        /* malloc interface */
        (*interfaces) = new Interfaces [NI];

        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                "please check the given interfaces mesh!\033[0m\n");
            fflush(stderr);
            exit(1);
        }
        
        if (NX < 2 || NY < 2) {
            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: No enough point (NX >= 2, NY >= 2)!");
            fflush(stderr);   
            exit(1);         
        }

        for (int ni = 0; ni < NI; ni++) {

            (*interfaces)[ni].NX   = NX;
            (*interfaces)[ni].NY   = NY;
            (*interfaces)[ni].DX   = DX;
            (*interfaces)[ni].DY   = DY;
            (*interfaces)[ni].MINX = MINX;
            (*interfaces)[ni].MINY = MINY;
                  
            (*interfaces)[ni].altitude = new float [NX*NY];
            (*interfaces)[ni].vp       = new float [NX*NY];
            (*interfaces)[ni].vs       = new float [NX*NY];
            (*interfaces)[ni].rho      = new float [NX*NY];
            (*interfaces)[ni].vp_grad  = new float [NX*NY];
            (*interfaces)[ni].vs_grad  = new float [NX*NY];
            (*interfaces)[ni].rho_grad = new float [NX*NY];

            for (size_t iy = 0; iy < NY; iy++) {
                for (size_t ix = 0; ix < NX; ix++) {
                    size_t indx = iy*NX + ix;

                    int u = fscanf(tmp_file, "%f %f %f %f %f %f %f",
                            &(*interfaces)[ni].altitude[indx], 
                            &(*interfaces)[ni].vp[indx],  &(*interfaces)[ni].vp_grad[indx],
                            &(*interfaces)[ni].vs[indx],  &(*interfaces)[ni].vs_grad[indx],
                            &(*interfaces)[ni].rho[indx], &(*interfaces)[ni].rho_grad[indx]);

                    if (u < 7) {
                        fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                            "please check the interfaces data of #%d layer!\033[0m\n", ni);
                        fflush(stderr);
                        exit(1);
                    }
                    
                }
            }
        }

        // end of the info
        break;
    }    
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
    Interfaces **interfaces)
{
    FILE *file = gfopen(grid_file, "r");
    FILE *tmp_file = tmpfile();

    char line[MAX_BUF_LEN];
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
            fprintf(stderr,"\033[47;31m The MEDIA GRID FILE is wrong, " \
                "please give a number of layers! \033[0m\n");
            fflush(stderr);
            exit(1);
        }

        if (NL < 1) {
            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: No enough layers (minimum is 1)!");
            fflush(stderr);
            exit(1);
        }

        for (int i = 0; i < NL; i++) {
            int ng_i = 0;
            if (fscanf(tmp_file, "%d",&ng_i) < 1) {
                fprintf(stderr,"\033[47;31m The MEDIA GRID FILE is wrong, " \
                    "please give the number of grids in the %d-th layer! \033[0m\n", i);
                fflush(stderr);
                exit(1);
            }
            NI += ng_i;
            NGz.push_back(ng_i);
        }

        /* malloc interface */
        (*interfaces) = new Interfaces [NI];


        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                "please check the given interfaces mesh!\033[0m\n");
            fflush(stderr);
            exit(1);
        }

        if (NX < 2 || NY < 2) {
            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: No enough point (NX >= 2, NY >= 2)!");
            fflush(stderr); 
            exit(1);           
        }


        /* The grid media XOY domain */

//        /* The given media domain must bigger than the calculation domain */
//        if (MAXX < Xmax || MINX > Xmin ) {
//            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: The given media range is smaller than "\
//                "the calculation grid range in x-direction!");
//            fflush(stderr); 
//            exit(1); 
//        }
//        if (MAXY < Ymax || MINY > Ymin) {
//            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: The given media range is smaller than "\
//                "the calculation grid range in y-direction!");
//            fflush(stderr); 
//            exit(1);             
//        }

        /* jlq TODO!!! */
        /* the range need to read */
        size_t ix0 = 0; //(Xmin-MINX)/DX;
        size_t iy0 = 0; //(Ymin-MINY)/DY;
        size_t ix1 = NX-1; //(Xmax-MINX)/DX+1;
        size_t iy1 = NY-1; //(Ymax-MINY)/DY+1;
        size_t nx  = ix1-ix0+1;
        size_t ny  = iy1-iy0+1;
        float minx = ix0 * DX + MINX;
        float miny = iy0 * DY + MINY;


        /* The given media domain must bigger than the calculation domain */
        if (ix0 < 0 || ix1 >= NX ) {
            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: The given media range is smaller than "\
                "the calculation grid range in x-direction!");
            fflush(stderr); 
            exit(1); 
        }
        if (iy0 < 0 || iy0 >= NY) {
            fprintf(stderr,"\033[47;31m%s\033[0m\n", "Error: The given media range is smaller than "\
                "the calculation grid range in y-direction!");
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
                  
            (*interfaces)[ni].altitude = new float [nx*ny];
            (*interfaces)[ni].vp       = new float [nx*ny];
            (*interfaces)[ni].vs       = new float [nx*ny];
            (*interfaces)[ni].rho      = new float [nx*ny];

            for (size_t iy = 0; iy < NY; iy++) {
                for (size_t ix = 0; ix < NX; ix++) {

                    if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {

                        size_t indx = (ix-ix0) + ny * (iy-iy0);

                        int u = fscanf(tmp_file, "%f %f %f %f",
                                       &(*interfaces)[ni].altitude[indx],  // altitude
                                       &(*interfaces)[ni].vp[indx],  
                                       &(*interfaces)[ni].vs[indx],  
                                       &(*interfaces)[ni].rho[indx] );
    
                        if (u < 4) {
                            fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                                "please check the grid data of #%d z-grid!\033[0m\n", ni);
                            fflush(stderr);
                            exit(1);
                        }
                    } else { // read the data no used.
                        float tmp[4];
                        int u = fscanf(tmp_file, "%f %f %f %f", tmp, tmp+1, tmp+2, tmp+3);
                        if (u < 4) {
                            fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                                "please check the grid data of #%d z-grid!\033[0m\n", ni);
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
}
