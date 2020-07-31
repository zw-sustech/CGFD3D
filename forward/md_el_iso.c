/*
********************************************************************************
* Curve grid metric calculation using MacCormack scheme                        *
********************************************************************************
*/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "netcdf.h"
#include "fdlib_mem.h"
#include "medium_elastic_iso.h"

int medium_elastic_iso_init_vars(size_t siz_volume, int *number_of_vars, 
      float **p_m3d, size_t **p_m3d_pos, char ***p_m3d_name)
{
    int ivar;

    *number_of_vars = 3;
    /*
     * 0: rho
     * 1: lambda
     * 2: mu
     */

    // vars
    *p_m3d = (float *) fdlib_mem_calloc_1d_float( 
                 siz_volume * (*number_of_vars), 0.0, "medium_elastic_iso_init_vars");

    if (*p_m3d == NULL) {
        fprintf(stderr,"Error: failed to alloc medium_elastic_iso\n");
        fflush(stderr);
        ierr = -1;
    }

    // position of each var
    *p_m3d_pos = (size_t *) fdlib_mem_calloc_1d_sizet( 
                 *number_of_vars, 0, "medium_elastic_iso_init_vars");

    // name of each var
    *p_m3d_name = (char **) fdlib_mem_malloc_2d_char( 
                 *number_of_vars, FDSYS_MAX_STR_LEN, "medium_elastic_iso_init_vars");

    // init
    ivar = MEDIUM_ELASTIC_ISO_SEQ_RHO;
    (*p_m3d_pos)[ivar] = ivar * siz_volume;
    strcpy((*p_m3d_name)[ivar],"rho");

    ivar = MEDIUM_ELASTIC_ISO_SEQ_LAMBDA;
    (*p_m3d_pos)[ivar] = ivar * siz_volume;
    strcpy((*p_m3d_name)[ivar],"lamda");

    ivar = MEDIUM_ELASTIC_ISO_SEQ_MU;
    (*p_m3d_pos)[ivar] = ivar * siz_volume;
    strcpy((*p_m3d_name)[ivar],"mu");

    return 0;
}

//
//
//
int medium_elastic_iso_import(float *restrict m3d, size_t *restrict m3d_pos, char **restrict m3d_name,
        int number_of_vars, size_t siz_volume, char *in_dir, int *myid3)
{
    int ierr = 0;
    char in_file[FDSYS_MAX_STR_LEN];

    int ncid, varid;

    // construct file name
    sprintf(in_file, "%s/media_mpi%02d%02d%02d.nc", in_dir, myid3[0], myid3[1], myid3[2]);

    // read in nc
    ierr = nc_open(in_file, NC_NOWRITE, &ncid);
    if (ierr != NC_NOERR) errh(ierr);

    for (int ivar=0; ivar<number_of_vars; ivar++) {
        ierr = nc_inq_varid(ncid, m3d_name[ivar], &varid);
        if (ierr != NC_NOERR) errh(ierr);

        ierr = nc_get_vara_float(ncid,varid,m3d+m3d_pos[ivar]);
        if (ierr != NC_NOERR) errh(ierr);
    }
    
    // close file
    ierr = nc_close(ncid);
    if (ierr != NC_NOERR) errh(ierr);

    return ierr;
}

