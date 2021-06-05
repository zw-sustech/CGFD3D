#include <iostream>
#include "media_read_interface_file.hpp"
//using namespace std;

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
            fprintf(stderr,"\033[47;31m%s\033[0m", "Error: No enough interfaes (minimum is 1)!\n");
            fflush(stderr);
            exit(1);
        }

        /* malloc interface */
        (*interfaces) = new Interfaces [NI];

        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"\033[47;31mThe INTERFACES FILE is wrong, " \
                "please check the given interfaces mesh\033[0m\n");
            fflush(stderr);
            exit(1);
        }
        
        if (NX < 2 || NY < 2) {
            fprintf(stderr,"\033[47;31m%s\033[0m", "Error: No enough point (NX >= 2, NY >= 2)!\n");
            fflush(stderr);            
        }

        for (int ni = 0; ni < NI; ni++) {

            (*interfaces)[ni].NX   = NX;
            (*interfaces)[ni].NY   = NY;
            (*interfaces)[ni].DX   = DX;
            (*interfaces)[ni].DY   = DY;
            (*interfaces)[ni].MINX = MINX;
            (*interfaces)[ni].MINY = MINY;
                  
            (*interfaces)[ni].depth    = new float [NX*NY];
            (*interfaces)[ni].vp       = new float [NX*NY];
            (*interfaces)[ni].vs       = new float [NX*NY];
            (*interfaces)[ni].rho      = new float [NX*NY];
            (*interfaces)[ni].vp_grad  = new float [NX*NY];
            (*interfaces)[ni].vs_grad  = new float [NX*NY];
            (*interfaces)[ni].rho_grad = new float [NX*NY];

            for (size_t ix = 0; ix < NX; ix++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    size_t indx = iy*NX + ix;

                    int u = fscanf(tmp_file, "%f %f %f %f %f %f %f",
                            &(*interfaces)[ni].depth[indx], 
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
