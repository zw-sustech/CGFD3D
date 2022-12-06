#include <iostream>
#include <string.h>
#include <math.h>
#include <map>
#include "media_read_file.hpp"
#include "media_utility.hpp"



FILE *gfopen(const char *filename, const char *mode)
{
    FILE *fp;
    if ((fp = fopen(filename,mode)) == NULL) {
        fprintf(stderr, "Error: Cannot open %s, " \
            "please check your file path and run-directory.\n",filename);
        exit(1);
    }
    return fp;
}


void read_interface_file(
    const char *interface_file,
    inter_t *interfaces)
{
    char line[MAX_BUF_LEN];
    FILE *file = gfopen(interface_file, "r");
    FILE *tmp_file = tmpfile();
    char media_type[MAX_BUF_LEN];
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
        if (fscanf(tmp_file, "%s", media_type) < 1) {
            fprintf(stderr,"Error: Please give a media_type in %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        } 

        if (strcmp(media_type, "one_component") == 0) {
            interfaces->media_type = ONE_COMPONENT; 
        } else if (strcmp(media_type, "acoustic_isotropic") == 0) {
            interfaces->media_type = ACOUSTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_isotropic") == 0) {
            interfaces->media_type = ELASTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_vti_prem") == 0) {
            interfaces->media_type = ELASTIC_VTI_PREM; 
        } else if (strcmp(media_type, "elastic_vti_thomsen") == 0) {
            interfaces->media_type = ELASTIC_VTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_vti_cij") == 0) {
            interfaces->media_type = ELASTIC_VTI_CIJ; 
        } else if (strcmp(media_type, "elastic_tti_thomsen") == 0) {
            interfaces->media_type = ELASTIC_TTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_tti_bond") == 0) {
            interfaces->media_type = ELASTIC_TTI_BOND;
        } else if (strcmp(media_type, "elastic_aniso_cij") == 0) {
            interfaces->media_type = ELASTIC_ANISO_CIJ; 
        } else {
            fprintf(stderr,"Error: media_type=%s is not supported, \n"\
                           "       please check %s!\n", media_type, interface_file);
            fflush(stderr);
            exit(1);
        }

        if (fscanf(tmp_file, "%d", &NI) < 1) {
            fprintf(stderr,"Error: please give a number of layers in %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        }

        if (NI < 1) {
            fprintf(stderr, "Error: No enough interfaes (minimum is 1), please check file %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        }

        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"Error: please check the given interfaces mesh in %s!\n", interface_file);
            fflush(stderr);
            exit(1);
        }

        if (NX < 2 || NY < 2) {
            fprintf(stderr, "Error: No enough point (NX >= 2, NY >= 2), please check file %s!\n", interface_file);
            fflush(stderr);   
            exit(1);         
        }

        // share the param of the interface structure
        interfaces->NI   = NI;
        interfaces->NX   = NX;
        interfaces->NY   = NY;
        interfaces->DX   = DX;
        interfaces->DY   = DY;
        interfaces->MINX = MINX;
        interfaces->MINY = MINY;
        

        /* volume and slice of interfaces info */
        size_t inter_line   =  NX;
        size_t inter_slice  =  NX*NY;
        size_t inter_volume =  NI*inter_slice;

        interfaces->elevation = new float[inter_volume];
        switch (interfaces -> media_type) 
        {
        case ONE_COMPONENT:
            interfaces->var      = new float[inter_volume];
            interfaces->var_grad = new float[inter_volume];
            interfaces->var_pow  = new float[inter_volume];
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f",
                                           &(interfaces->elevation[indx]), 
                                           &(interfaces->var[indx]      ),
                                           &(interfaces->var_grad[indx] ),
                                           &(interfaces->var_pow[indx]  ) );
                        if (num_read < 4) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ACOUSTIC_ISOTROPIC:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];
            interfaces->vp      = new float[inter_volume];
            interfaces->vp_grad = new float[inter_volume];
            interfaces->vp_pow  = new float[inter_volume];
       
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f",
                                          &(interfaces->elevation[indx]), 
                                          &(interfaces->rho[indx]      ),
                                          &(interfaces->rho_grad[indx] ),
                                          &(interfaces->rho_pow[indx]  ),
                                          &(interfaces->vp[indx]       ),
                                          &(interfaces->vp_grad[indx]  ),
                                          &(interfaces->vp_pow[indx]   ) );
                        if (num_read < 7) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_ISOTROPIC:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];
            interfaces->vp      = new float[inter_volume];
            interfaces->vp_grad = new float[inter_volume];
            interfaces->vp_pow  = new float[inter_volume];
            interfaces->vs      = new float[inter_volume];
            interfaces->vs_grad = new float[inter_volume];
            interfaces->vs_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f %f %f",
                                          &(interfaces->elevation[indx]), 
                                          &(interfaces->rho[indx]      ),
                                          &(interfaces->rho_grad[indx] ),
                                          &(interfaces->rho_pow[indx]  ),
                                          &(interfaces->vp[indx]       ),
                                          &(interfaces->vp_grad[indx]  ),
                                          &(interfaces->vp_pow[indx]   ),
                                          &(interfaces->vs[indx]       ),
                                          &(interfaces->vs_grad[indx]  ),
                                          &(interfaces->vs_pow[indx]   ) );

                        if (num_read < 10) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_VTI_PREM:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];

            interfaces->vph      = new float[inter_volume];
            interfaces->vph_grad = new float[inter_volume];
            interfaces->vph_pow  = new float[inter_volume];

            interfaces->vpv      = new float[inter_volume];
            interfaces->vpv_grad = new float[inter_volume];
            interfaces->vpv_pow  = new float[inter_volume];

            interfaces->vsh      = new float[inter_volume];
            interfaces->vsh_grad = new float[inter_volume];
            interfaces->vsh_pow  = new float[inter_volume];

            interfaces->vsv      = new float[inter_volume];
            interfaces->vsv_grad = new float[inter_volume];
            interfaces->vsv_pow  = new float[inter_volume];

            interfaces->eta      = new float[inter_volume];
            interfaces->eta_grad = new float[inter_volume];
            interfaces->eta_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf(tmp_file,
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                            &(interfaces->elevation[indx]), 
                            &(interfaces->rho[indx]      ),
                            &(interfaces->rho_grad[indx] ),
                            &(interfaces->rho_pow[indx]  ),
                            &(interfaces->vph[indx]      ),
                            &(interfaces->vph_grad[indx] ),
                            &(interfaces->vph_pow[indx]  ),
                            &(interfaces->vpv[indx]      ),
                            &(interfaces->vpv_grad[indx] ),
                            &(interfaces->vpv_pow[indx]  ),
                            &(interfaces->vsh[indx]      ),
                            &(interfaces->vsh_grad[indx] ),
                            &(interfaces->vsh_pow[indx]  ),
                            &(interfaces->vsv[indx]      ),
                            &(interfaces->vsv_grad[indx] ),
                            &(interfaces->vsv_pow[indx]  ),
                            &(interfaces->eta[indx]      ),
                            &(interfaces->eta_grad[indx] ),
                            &(interfaces->eta_pow[indx]  ) );

                        if (num_read < 19) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_VTI_THOMSEN:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];

            interfaces->vp0      = new float[inter_volume];
            interfaces->vp0_grad = new float[inter_volume];
            interfaces->vp0_pow  = new float[inter_volume];

            interfaces->vs0      = new float[inter_volume];
            interfaces->vs0_grad = new float[inter_volume];
            interfaces->vs0_pow  = new float[inter_volume];

            interfaces->epsilon      = new float[inter_volume];
            interfaces->epsilon_grad = new float[inter_volume];
            interfaces->epsilon_pow  = new float[inter_volume];

            interfaces->delta      = new float[inter_volume];
            interfaces->delta_grad = new float[inter_volume];
            interfaces->delta_pow  = new float[inter_volume];

            interfaces->gamma      = new float[inter_volume];
            interfaces->gamma_grad = new float[inter_volume];
            interfaces->gamma_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf(tmp_file, 
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                            &(interfaces->elevation[indx]), 
                            &(interfaces->rho[indx]      ),
                            &(interfaces->rho_grad[indx] ),
                            &(interfaces->rho_pow[indx]  ),
                            &(interfaces->vp0[indx]      ),
                            &(interfaces->vp0_grad[indx] ),
                            &(interfaces->vp0_pow[indx]  ),
                            &(interfaces->vs0[indx]      ),
                            &(interfaces->vs0_grad[indx] ),
                            &(interfaces->vs0_pow[indx]  ),
                            &(interfaces->epsilon[indx]      ),
                            &(interfaces->epsilon_grad[indx] ),
                            &(interfaces->epsilon_pow[indx]  ),
                            &(interfaces->delta[indx]      ),
                            &(interfaces->delta_grad[indx] ),
                            &(interfaces->delta_pow[indx]  ),
                            &(interfaces->gamma[indx]      ),
                            &(interfaces->gamma_grad[indx] ),
                            &(interfaces->gamma_pow[indx]  ) );

                        if (num_read < 19) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_VTI_CIJ:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];

            interfaces->c11      = new float[inter_volume];
            interfaces->c11_grad = new float[inter_volume];
            interfaces->c11_pow  = new float[inter_volume];

            interfaces->c33      = new float[inter_volume];
            interfaces->c33_grad = new float[inter_volume];
            interfaces->c33_pow  = new float[inter_volume];

            interfaces->c55      = new float[inter_volume];
            interfaces->c55_grad = new float[inter_volume];
            interfaces->c55_pow  = new float[inter_volume];

            interfaces->c66      = new float[inter_volume];
            interfaces->c66_grad = new float[inter_volume];
            interfaces->c66_pow  = new float[inter_volume];

            interfaces->c13      = new float[inter_volume];
            interfaces->c13_grad = new float[inter_volume];
            interfaces->c13_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, 
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                             &(interfaces->elevation[indx]), 
                             &(interfaces->rho[indx]      ),
                             &(interfaces->rho_grad[indx] ),
                             &(interfaces->rho_pow[indx]  ),
                             &(interfaces->c11[indx]      ),
                             &(interfaces->c11_grad[indx] ),
                             &(interfaces->c11_pow[indx]  ),
                             &(interfaces->c33[indx]      ),
                             &(interfaces->c33_grad[indx] ),
                             &(interfaces->c33_pow[indx]  ),
                             &(interfaces->c55[indx]      ),
                             &(interfaces->c55_grad[indx] ),
                             &(interfaces->c55_pow[indx]  ),
                             &(interfaces->c66[indx]      ),
                             &(interfaces->c66_grad[indx] ),
                             &(interfaces->c66_pow[indx]  ),
                             &(interfaces->c13[indx]      ),
                             &(interfaces->c13_grad[indx] ),
                             &(interfaces->c13_pow[indx]  ) );

                        if (num_read < 19) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_TTI_THOMSEN:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];

            interfaces->vp0      = new float[inter_volume];
            interfaces->vp0_grad = new float[inter_volume];
            interfaces->vp0_pow  = new float[inter_volume];

            interfaces->vs0      = new float[inter_volume];
            interfaces->vs0_grad = new float[inter_volume];
            interfaces->vs0_pow  = new float[inter_volume];

            interfaces->epsilon      = new float[inter_volume];
            interfaces->epsilon_grad = new float[inter_volume];
            interfaces->epsilon_pow  = new float[inter_volume];

            interfaces->delta      = new float[inter_volume];
            interfaces->delta_grad = new float[inter_volume];
            interfaces->delta_pow  = new float[inter_volume];

            interfaces->gamma      = new float[inter_volume];
            interfaces->gamma_grad = new float[inter_volume];
            interfaces->gamma_pow  = new float[inter_volume];

            interfaces->azimuth      = new float[inter_volume];
            interfaces->azimuth_grad = new float[inter_volume];
            interfaces->azimuth_pow  = new float[inter_volume];

            interfaces->dip      = new float[inter_volume];
            interfaces->dip_grad = new float[inter_volume];
            interfaces->dip_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, 
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                            &(interfaces->elevation[indx]), 
                            &(interfaces->rho[indx]      ),
                            &(interfaces->rho_grad[indx] ),
                            &(interfaces->rho_pow[indx]  ),
                            &(interfaces->vp0[indx]      ),
                            &(interfaces->vp0_grad[indx] ),
                            &(interfaces->vp0_pow[indx]  ),
                            &(interfaces->vs0[indx]      ),
                            &(interfaces->vs0_grad[indx] ),
                            &(interfaces->vs0_pow[indx]  ),
                            &(interfaces->epsilon[indx]      ),
                            &(interfaces->epsilon_grad[indx] ),
                            &(interfaces->epsilon_pow[indx]  ),
                            &(interfaces->delta[indx]      ),
                            &(interfaces->delta_grad[indx] ),
                            &(interfaces->delta_pow[indx]  ),
                            &(interfaces->gamma[indx]      ),
                            &(interfaces->gamma_grad[indx] ),
                            &(interfaces->gamma_pow[indx]  ),
                            &(interfaces->azimuth[indx]      ),
                            &(interfaces->azimuth_grad[indx] ),
                            &(interfaces->azimuth_pow[indx]  ),
                            &(interfaces->dip[indx]      ),
                            &(interfaces->dip_grad[indx] ),
                            &(interfaces->dip_pow[indx]  ) );

                        if (num_read < 25) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_TTI_BOND:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];

            interfaces->c11      = new float[inter_volume];
            interfaces->c11_grad = new float[inter_volume];
            interfaces->c11_pow  = new float[inter_volume];

            interfaces->c33      = new float[inter_volume];
            interfaces->c33_grad = new float[inter_volume];
            interfaces->c33_pow  = new float[inter_volume];

            interfaces->c55      = new float[inter_volume];
            interfaces->c55_grad = new float[inter_volume];
            interfaces->c55_pow  = new float[inter_volume];

            interfaces->c66      = new float[inter_volume];
            interfaces->c66_grad = new float[inter_volume];
            interfaces->c66_pow  = new float[inter_volume];

            interfaces->c13      = new float[inter_volume];
            interfaces->c13_grad = new float[inter_volume];
            interfaces->c13_pow  = new float[inter_volume];

            interfaces->azimuth      = new float[inter_volume];
            interfaces->azimuth_grad = new float[inter_volume];
            interfaces->azimuth_pow  = new float[inter_volume];

            interfaces->dip      = new float[inter_volume];
            interfaces->dip_grad = new float[inter_volume];
            interfaces->dip_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, 
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                            &(interfaces->elevation[indx]), 
                            &(interfaces->rho[indx]      ),
                            &(interfaces->rho_grad[indx] ),
                            &(interfaces->rho_pow[indx]  ),
                            &(interfaces->c11[indx]      ),
                            &(interfaces->c11_grad[indx] ),
                            &(interfaces->c11_pow[indx]  ),
                            &(interfaces->c33[indx]      ),
                            &(interfaces->c33_grad[indx] ),
                            &(interfaces->c33_pow[indx]  ),
                            &(interfaces->c55[indx]      ),
                            &(interfaces->c55_grad[indx] ),
                            &(interfaces->c55_pow[indx]  ),
                            &(interfaces->c66[indx]      ),
                            &(interfaces->c66_grad[indx] ),
                            &(interfaces->c66_pow[indx]  ),
                            &(interfaces->c13[indx]      ),
                            &(interfaces->c13_grad[indx] ),
                            &(interfaces->c13_pow[indx]  ),
                            &(interfaces->azimuth[indx]      ),
                            &(interfaces->azimuth_grad[indx] ),
                            &(interfaces->azimuth_pow[indx]  ),
                            &(interfaces->dip[indx]      ),
                            &(interfaces->dip_grad[indx] ),
                            &(interfaces->dip_pow[indx]  ));

                        if (num_read < 25) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;

        case ELASTIC_ANISO_CIJ:
            interfaces->rho      = new float[inter_volume];
            interfaces->rho_grad = new float[inter_volume];
            interfaces->rho_pow  = new float[inter_volume];

            interfaces->c11      = new float[inter_volume]; interfaces->c12      = new float[inter_volume];     
            interfaces->c11_grad = new float[inter_volume]; interfaces->c12_grad = new float[inter_volume];     
            interfaces->c11_pow  = new float[inter_volume]; interfaces->c12_pow  = new float[inter_volume]; 

            interfaces->c13      = new float[inter_volume]; interfaces->c14      = new float[inter_volume]; 
            interfaces->c13_grad = new float[inter_volume]; interfaces->c14_grad = new float[inter_volume]; 
            interfaces->c13_pow  = new float[inter_volume]; interfaces->c14_pow  = new float[inter_volume];     

            interfaces->c15      = new float[inter_volume]; interfaces->c16      = new float[inter_volume];    
            interfaces->c15_grad = new float[inter_volume]; interfaces->c16_grad = new float[inter_volume];    
            interfaces->c15_pow  = new float[inter_volume]; interfaces->c16_pow  = new float[inter_volume];    

            interfaces->c22      = new float[inter_volume]; interfaces->c23      = new float[inter_volume];   
            interfaces->c22_grad = new float[inter_volume]; interfaces->c23_grad = new float[inter_volume];   
            interfaces->c22_pow  = new float[inter_volume]; interfaces->c23_pow  = new float[inter_volume];   

            interfaces->c24      = new float[inter_volume]; interfaces->c25      = new float[inter_volume]; 
            interfaces->c24_grad = new float[inter_volume]; interfaces->c25_grad = new float[inter_volume]; 
            interfaces->c24_pow  = new float[inter_volume]; interfaces->c25_pow  = new float[inter_volume]; 

            interfaces->c26      = new float[inter_volume];    
            interfaces->c26_grad = new float[inter_volume];    
            interfaces->c26_pow  = new float[inter_volume];

            interfaces->c33      = new float[inter_volume]; interfaces->c34      = new float[inter_volume];
            interfaces->c33_grad = new float[inter_volume]; interfaces->c34_grad = new float[inter_volume];
            interfaces->c33_pow  = new float[inter_volume]; interfaces->c34_pow  = new float[inter_volume];

            interfaces->c35      = new float[inter_volume]; interfaces->c36      = new float[inter_volume];
            interfaces->c35_grad = new float[inter_volume]; interfaces->c36_grad = new float[inter_volume];
            interfaces->c35_pow  = new float[inter_volume]; interfaces->c36_pow  = new float[inter_volume];

            interfaces->c44      = new float[inter_volume]; interfaces->c45      = new float[inter_volume];
            interfaces->c44_grad = new float[inter_volume]; interfaces->c45_grad = new float[inter_volume];
            interfaces->c44_pow  = new float[inter_volume]; interfaces->c45_pow  = new float[inter_volume];

            interfaces->c46      = new float[inter_volume]; interfaces->c55      = new float[inter_volume];
            interfaces->c46_grad = new float[inter_volume]; interfaces->c55_grad = new float[inter_volume]; 
            interfaces->c46_pow  = new float[inter_volume]; interfaces->c55_pow  = new float[inter_volume];

            interfaces->c56      = new float[inter_volume]; interfaces->c66      = new float[inter_volume]; 
            interfaces->c56_grad = new float[inter_volume]; interfaces->c66_grad = new float[inter_volume];
            interfaces->c56_pow  = new float[inter_volume]; interfaces->c66_pow  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t j = 0; j < NY; j++) {
                    for (size_t i = 0; i < NX; i++) {
                        size_t indx = i + j * inter_line + ni*inter_slice;
                        int num_read = 0;
                        num_read = fscanf( tmp_file, "%f %f %f %f"\
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"\
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f"\
                            "%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",
                            &(interfaces->elevation[indx]),  
                            &(interfaces->rho[indx]), &(interfaces->rho_grad[indx]), &(interfaces->rho_pow[indx]),
                            &(interfaces->c11[indx]), &(interfaces->c11_grad[indx]), &(interfaces->c11_pow[indx]),
                            &(interfaces->c12[indx]), &(interfaces->c12_grad[indx]), &(interfaces->c12_pow[indx]),
                            &(interfaces->c13[indx]), &(interfaces->c13_grad[indx]), &(interfaces->c13_pow[indx]),
                            &(interfaces->c14[indx]), &(interfaces->c14_grad[indx]), &(interfaces->c14_pow[indx]),
                            &(interfaces->c15[indx]), &(interfaces->c15_grad[indx]), &(interfaces->c15_pow[indx]),
                            &(interfaces->c16[indx]), &(interfaces->c16_grad[indx]), &(interfaces->c16_pow[indx]),
                            &(interfaces->c22[indx]), &(interfaces->c22_grad[indx]), &(interfaces->c22_pow[indx]),
                            &(interfaces->c23[indx]), &(interfaces->c23_grad[indx]), &(interfaces->c23_pow[indx]),
                            &(interfaces->c24[indx]), &(interfaces->c24_grad[indx]), &(interfaces->c24_pow[indx]),
                            &(interfaces->c25[indx]), &(interfaces->c25_grad[indx]), &(interfaces->c25_pow[indx]),
                            &(interfaces->c26[indx]), &(interfaces->c26_grad[indx]), &(interfaces->c26_pow[indx]),
                            &(interfaces->c33[indx]), &(interfaces->c33_grad[indx]), &(interfaces->c33_pow[indx]),
                            &(interfaces->c34[indx]), &(interfaces->c34_grad[indx]), &(interfaces->c34_pow[indx]),
                            &(interfaces->c35[indx]), &(interfaces->c35_grad[indx]), &(interfaces->c35_pow[indx]),
                            &(interfaces->c36[indx]), &(interfaces->c36_grad[indx]), &(interfaces->c36_pow[indx]),
                            &(interfaces->c44[indx]), &(interfaces->c44_grad[indx]), &(interfaces->c44_pow[indx]),
                            &(interfaces->c45[indx]), &(interfaces->c45_grad[indx]), &(interfaces->c45_pow[indx]),
                            &(interfaces->c46[indx]), &(interfaces->c46_grad[indx]), &(interfaces->c46_pow[indx]), 
                            &(interfaces->c55[indx]), &(interfaces->c55_grad[indx]), &(interfaces->c55_pow[indx]),
                            &(interfaces->c56[indx]), &(interfaces->c56_grad[indx]), &(interfaces->c56_pow[indx]),                     
                            &(interfaces->c66[indx]), &(interfaces->c66_grad[indx]), &(interfaces->c66_pow[indx]) );

                        if (num_read < 67) {
                            fprintf(stderr, "Error: Insufficient data in %s.\n", interface_file);
                            fflush(stderr);
                            exit(1);
                        }
                        
                    }
                }
            }
        break;
        default: // for self-check
            fprintf(stderr,"Error: Unknow media #%d (code check, please contact Luqian Jiang).\n", interfaces->media_type);
            fflush(stderr);
            exit(1);

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
    inter_t *interfaces)
{
    FILE *tmp_file = gfopen(grid_file, "r");
//    FILE *tmp_file = tmpfile();

    char  line[MAX_BUF_LEN];
    char  media_type[MAX_BUF_LEN];
    int NX = 0, NY = 0; 
    size_t NI = 0;
    float DX = 0, DY = 0;
    float MINX = 0, MINY = 0;

    /* Read every line in the file. */
//    while(fgets(line, MAX_BUF_LEN, file) != NULL)
//    {
//        if (line[0] == '#' || line[0] == '\n')
//            continue;
//        fputs(line, tmp_file);
//    } 
//
//    /* Set the file pointer at the beginning of the stream*/
//    rewind(tmp_file);

    /* Read the temporary data file, Can not use the feof() */
    while(feof(tmp_file) != EOF)
    {

        if (fscanf(tmp_file, "%s", media_type) < 1) {
            fprintf(stderr,"Error: The GRID MEIDA FILE is wrong, please give a media_type!\n");
            fflush(stderr);
            exit(1);
        } 

        if (strcmp(media_type, "one_component") == 0) {
            interfaces->media_type = ONE_COMPONENT; 
        } else if (strcmp(media_type, "acoustic_isotropic") == 0) {
            interfaces->media_type = ACOUSTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_isotropic") == 0) {
            interfaces->media_type = ELASTIC_ISOTROPIC; 
        } else if (strcmp(media_type, "elastic_vti_prem") == 0) {
            interfaces->media_type = ELASTIC_VTI_PREM; 
        } else if (strcmp(media_type, "elastic_vti_thomsen") == 0) {
            interfaces->media_type = ELASTIC_VTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_vti_cij") == 0) {
            interfaces->media_type = ELASTIC_VTI_CIJ; 
        } else if (strcmp(media_type, "elastic_tti_thomsen") == 0) {
            interfaces->media_type = ELASTIC_TTI_THOMSEN; 
        } else if (strcmp(media_type, "elastic_tti_bond") == 0) {
            interfaces->media_type = ELASTIC_TTI_BOND;
        } else if (strcmp(media_type, "elastic_aniso_cij") == 0) {
            interfaces->media_type = ELASTIC_ANISO_CIJ; 
        } else {
            fprintf(stderr,"Error: media_type = %s is not supported, \n"\
                           "       please check %s!\n", 
                           media_type, grid_file);
            fflush(stderr);
            exit(1);
        }

        if (fscanf(tmp_file, "%d", &NL) < 1) {
            fprintf(stderr,"Error: The GRID MEIDA FILE is wrong, " \
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
                fprintf(stderr,"Error: The GRID MEIDA FILE is wrong, " \
                    "please give the number of grids in the %d-th layer! \n", i);
                fflush(stderr);
                exit(1);
            }
            NI += ng_i;
            NGz.push_back(ng_i);
        }

        if (fscanf(tmp_file, "%d %d %f %f %f %f", &NX, &NY, &MINX, &MINY, &DX, &DY) < 6) {
            fprintf(stderr,"Error: The GRID MEIDA FILE is wrong, " \
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
        size_t ix1 = ceil((Xmax-MINX)/DX);
        size_t iy1 = ceil((Ymax-MINY)/DY);
        size_t nx  = ix1-ix0+1;
        size_t ny  = iy1-iy0+1;
        float minx = ix0 * DX + MINX;
        float miny = iy0 * DY + MINY;


        /* The given media domain must bigger than the calculation domain */
        if (ix0 < 0 || ix1 > NX ) {
            fprintf(stderr,"Error: The given media range is smaller than "\
                "the calculation grid range in x-direction! \n");
            fflush(stderr); 
            exit(1); 
        }
        if (iy0 < 0 || iy0 > NY) {
            fprintf(stderr, "Error: The given media range is smaller than "\
                "the calculation grid range in y-direction! \n");
            fflush(stderr); 
            exit(1);             
        }

        interfaces -> NI = NI;
        interfaces -> NX = nx;
        interfaces -> NY = ny;
       // interfaces -> NX = NX;
       // interfaces -> NY = NY;
        interfaces -> DX = DX; 
        interfaces -> DY = DY; 
        interfaces -> MINX = minx;
        interfaces -> MINY = miny;

        // for read interface info
        size_t inter_line = nx;
        size_t inter_slice = nx*ny;
        size_t inter_volume = inter_slice*NI;

        interfaces->elevation = new float[inter_volume];

        switch (interfaces -> media_type) 
        {
        case ONE_COMPONENT:
            interfaces->var = new float[inter_volume];
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {

                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {

                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f",
                                               &(interfaces->elevation[indx]), 
                                               &(interfaces->var[indx]      ));
                            if (num_read < 2) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
                        } else { // read the data no used.
                            float tmp[2];
                            int u = fscanf(tmp_file, "%f %f", tmp, tmp+1);
                            if (u < 2) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }
                        
                    }
                }
            }
        break;

        case ACOUSTIC_ISOTROPIC:
            interfaces->rho = new float[inter_volume];
            interfaces->vp  = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {

                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->vp[indx]       ));
                            if (num_read < 3) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
                        } else { // read the data no used.
                            float tmp[3];
                            int u = fscanf(tmp_file, "%f %f %f", tmp, tmp+1, tmp+2);
                            if (u < 3) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        
                        }
                    }
                }
            }
        break;

        case ELASTIC_ISOTROPIC:
            interfaces->rho = new float[inter_volume];
            interfaces->vp  = new float[inter_volume];
            interfaces->vs  = new float[inter_volume];
            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {

                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f\n",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->vp[indx]       ),
                                              &(interfaces->vs[indx]       ));
                            if (num_read < 4) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
                        } else { // read the data no used.
                            float tmp[4];
                            int u = fscanf(tmp_file, "%f %f %f %f\n", tmp, tmp+1, tmp+2, tmp+3);
                            if (u < 4) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                                
                        }
                    }
                }
            }
        break;

        case ELASTIC_VTI_PREM:
            interfaces->rho = new float[inter_volume];
            interfaces->vph = new float[inter_volume];
            interfaces->vpv = new float[inter_volume];
            interfaces->vsh = new float[inter_volume];
            interfaces->vsv = new float[inter_volume];
            interfaces->eta = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {

                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->vph[indx]      ),
                                              &(interfaces->vpv[indx]      ),
                                              &(interfaces->vsh[indx]      ),
                                              &(interfaces->vsv[indx]      ),
                                              &(interfaces->eta[indx]      ) );
                            if (num_read < 7) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
                        
                        } else { // read the data no used.
                            float tmp[7];
                            int u = fscanf(tmp_file, "%f %f %f %f %f %f %f", 
                                tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6);
                            if (u < 7) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }

                    }
                }
            }
        break;

        case ELASTIC_VTI_THOMSEN:
            interfaces->rho = new float[inter_volume];
            interfaces->vp0 = new float[inter_volume];
            interfaces->vs0 = new float[inter_volume];
            interfaces->epsilon = new float[inter_volume];
            interfaces->delta = new float[inter_volume];
            interfaces->gamma = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {
                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->vp0[indx]      ),
                                              &(interfaces->vs0[indx]      ),
                                              &(interfaces->epsilon[indx]  ),
                                              &(interfaces->delta[indx]    ),
                                              &(interfaces->gamma[indx]    ) );
                            if (num_read < 7) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        } else { // read the data no used.
                            float tmp[7];
                            int u = fscanf(tmp_file, "%f %f %f %f %f %f %f", 
                                tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6);
                            if (u < 7) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }

                    }
                }
            }
        break;

        case ELASTIC_VTI_CIJ:
            interfaces->rho = new float[inter_volume];
            interfaces->c11 = new float[inter_volume];
            interfaces->c33 = new float[inter_volume];
            interfaces->c55 = new float[inter_volume];
            interfaces->c66 = new float[inter_volume];
            interfaces->c13 = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {
                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->c11[indx]      ),
                                              &(interfaces->c33[indx]      ),
                                              &(interfaces->c55[indx]      ),
                                              &(interfaces->c66[indx]      ),
                                              &(interfaces->c13[indx]      ) );
                            if (num_read < 7) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        
                        } else { // read the data no used.
                            float tmp[7];
                            int u = fscanf(tmp_file, "%f %f %f %f %f %f %f", 
                                tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6);
                            if (u < 7) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }
                    }
                }
            }
        break; 

        case ELASTIC_TTI_THOMSEN:
            interfaces->rho = new float[inter_volume];
            interfaces->vp0 = new float[inter_volume];
            interfaces->vs0 = new float[inter_volume];
            interfaces->epsilon = new float[inter_volume];
            interfaces->delta = new float[inter_volume];
            interfaces->gamma = new float[inter_volume];
            interfaces->azimuth = new float[inter_volume];
            interfaces->dip = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {
                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f %f",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->vp0[indx]      ),
                                              &(interfaces->vs0[indx]      ),
                                              &(interfaces->epsilon[indx]  ),
                                              &(interfaces->delta[indx]    ),
                                              &(interfaces->gamma[indx]    ),
                                              &(interfaces->azimuth[indx]  ),
                                              &(interfaces->dip[indx]      ) );
                            if (num_read < 9) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
                        } else { // read the data no used.
                            float tmp[9];
                            int u = fscanf(tmp_file, "%f %f %f %f %f %f %f %f %f", 
                                tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6, tmp+7, tmp+8);
                            if (u < 9) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }
                    }
                }
            }
        break;               

        case ELASTIC_TTI_BOND:
            interfaces->rho = new float[inter_volume];
            interfaces->c11 = new float[inter_volume];
            interfaces->c33 = new float[inter_volume];
            interfaces->c55 = new float[inter_volume];
            interfaces->c66 = new float[inter_volume];
            interfaces->c13 = new float[inter_volume];
            interfaces->azimuth = new float[inter_volume];
            interfaces->dip = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {
                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f %f",
                                              &(interfaces->elevation[indx]), 
                                              &(interfaces->rho[indx]      ),
                                              &(interfaces->c11[indx]      ),
                                              &(interfaces->c33[indx]      ),
                                              &(interfaces->c55[indx]      ),
                                              &(interfaces->c66[indx]      ),
                                              &(interfaces->c13[indx]      ),
                                              &(interfaces->azimuth[indx]  ),
                                              &(interfaces->dip[indx]      ) );
                            if (num_read < 9) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
    
                        } else { // read the data no used.
                            float tmp[9];
                            int u = fscanf(tmp_file, "%f %f %f %f %f %f %f %f %f", 
                                tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6, tmp+7, tmp+8);
                            if (u < 9) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }
                    }
                }
            }
        break;

        case ELASTIC_ANISO_CIJ:
            interfaces->rho = new float[inter_volume];
            interfaces->c11 = new float[inter_volume];
            interfaces->c12 = new float[inter_volume];
            interfaces->c13 = new float[inter_volume];
            interfaces->c14 = new float[inter_volume];
            interfaces->c15 = new float[inter_volume];
            interfaces->c16 = new float[inter_volume];

            interfaces->c22 = new float[inter_volume];
            interfaces->c23 = new float[inter_volume];
            interfaces->c24 = new float[inter_volume];
            interfaces->c25 = new float[inter_volume];
            interfaces->c26 = new float[inter_volume];

            interfaces->c33 = new float[inter_volume];
            interfaces->c34 = new float[inter_volume];
            interfaces->c35 = new float[inter_volume];
            interfaces->c36 = new float[inter_volume];

            interfaces->c44 = new float[inter_volume];
            interfaces->c45 = new float[inter_volume];
            interfaces->c46 = new float[inter_volume];
            interfaces->c55 = new float[inter_volume];
            interfaces->c56 = new float[inter_volume];
            interfaces->c66 = new float[inter_volume];

            for (size_t ni = 0; ni < NI; ni++) {
                for (size_t iy = 0; iy < NY; iy++) {
                    for (size_t ix = 0; ix < NX; ix++) {
                        if (ix >= ix0 && ix <= ix1 && iy >= iy0 && iy <= iy1) {
                            size_t indx = (ix-ix0) + (iy-iy0) * inter_line + ni*inter_slice;
                            int num_read = 0;
                            num_read = fscanf( tmp_file, "%f %f %f %f %f %f %f %f %f %f" \
                                                "%f %f %f %f %f %f %f %f %f %f %f %f %f" ,
                                                &(interfaces->elevation[indx]),
                                                &(interfaces->rho[indx]),
                                                &(interfaces->c11[indx]),
                                                &(interfaces->c12[indx]),
                                                &(interfaces->c13[indx]),
                                                &(interfaces->c14[indx]),
                                                &(interfaces->c15[indx]),
                                                &(interfaces->c16[indx]),
                                                &(interfaces->c22[indx]),
                                                &(interfaces->c23[indx]),
                                                &(interfaces->c24[indx]),
                                                &(interfaces->c25[indx]),
                                                &(interfaces->c26[indx]),
                                                &(interfaces->c33[indx]),
                                                &(interfaces->c34[indx]),
                                                &(interfaces->c35[indx]),
                                                &(interfaces->c36[indx]),
                                                &(interfaces->c44[indx]),
                                                &(interfaces->c45[indx]),
                                                &(interfaces->c46[indx]),
                                                &(interfaces->c55[indx]),
                                                &(interfaces->c56[indx]),
                                                &(interfaces->c66[indx]) );
                            if (num_read < 23) {
                                fprintf(stderr, "Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            } 
                        
                        } else { // read the data no used.
                            float tmp[23];
                            int u = fscanf( tmp_file, "%f %f %f %f %f %f %f %f %f %f" \
                                            "%f %f %f %f %f %f %f %f %f %f %f %f %f" ,
                                            tmp, tmp+1, tmp+2, tmp+3, tmp+4, tmp+5, tmp+6, tmp+7,
                                            tmp+8, tmp+9, tmp+10, tmp+11, tmp+12, tmp+13, tmp+14,
                                            tmp+15, tmp+16, tmp+17, tmp+18, tmp+19, tmp+20, tmp+21, tmp+22);
                            if (u < 23) {
                                fprintf(stderr,"Error: Insufficient data in %s.\n", grid_file);
                                fflush(stderr);
                                exit(1);
                            }
                        }
                    }
                }
            }
        break;

        default:
            fprintf(stderr, "Error: Unknow media_type %s (for code check)\n", media_type);
            fflush(stderr);
            exit(1);
        }
        break;
    }

    checkGridData(NL, NGz, *interfaces, grid_file);

//    fclose(file);   
    fclose(tmp_file);    
}

// check whether the elevation[ng[i]-1] == elevation[ng[i]] 
int checkGridData(int NL, 
    std::vector <int> &NGz,
    inter_t &interfaces, 
    const char *grid_file) 
{   
    if (NL < 2) return 0;

    int inter_slice = interfaces.NX*interfaces.NY;
    int ipoint = 0;
    for (int i = 0; i < NL-1; i++) {
        ipoint += (NGz[i]); 
        for (int indx = 0; indx < inter_slice; indx++) {
            int indx_top = indx + (ipoint-1)*inter_slice;
            int indx_bot = indx + ipoint*inter_slice;
            if (interfaces.elevation[indx_top] != interfaces.elevation[indx_bot]) {
                fprintf(stderr, "Error: The last elevations of #%d layer should equal to the first "\
                                "elevations of #%d layer!, please check %s!", i, i+1, grid_file);        
                fflush(stderr);
                exit(1);
            }
        }      
    }

    return 0;
}


/* 
 * Just read the grid data within the given
 *  [Xmin, Xmax]\times[Ymin, Ymax] domain.
 */
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
    size_t bin_slice)
{
    FILE *fid = gfopen(bin_file, "rb");
    std::vector<int> i(3, 0);

    for (i[2] = 0; i[2] < bin_size[2]; i[2]++) {
      for (i[1] = 0; i[1] < bin_size[1]; i[1]++) {
        for (i[0] = 0; i[0] < bin_size[0]; i[0]++) {
          if (i[0] >= bin_start[0] && i[0] <= bin_end[0] &&
              i[1] >= bin_start[1] && i[1] <= bin_end[1] &&
              i[2] >= bin_start[2] && i[2] <= bin_end[2]) 
          {
            size_t indx = (i[dimx] - bin_start[dimx]) + 
                          (i[dimy] - bin_start[dimy]) * bin_line + 
                           i[dimz] * bin_slice;

            int u = fread(&var[indx], sizeof(float), 1, fid);
            if (u < 1) {
                fprintf(stderr,"Error: Insufficient data in %s.\n",bin_file);
                fflush(stderr);
                exit(1);
            }
          } else {
            float tmp;
            int u = fread(&tmp, sizeof(float), 1, fid);
            if (u < 1) {
                fprintf(stderr,"Error: Insufficient data in %s.\n",bin_file);
                fflush(stderr);
                exit(1);
            }
          }
        }
      }
    }

    fclose(fid);
}
