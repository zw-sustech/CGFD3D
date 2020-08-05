/*********************************************************************
 Main program for seismic wave propagation simulation 

 AUTHOR:
     ZANG Nan
    ZHANG Wei

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <mpi.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "par_funcs.h"
#include "blk_t.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "sv_eliso1st_curv_macdrp.h"

/*******************************************************************************
 * usage function
 ******************************************************************************/

int cgfd3d_elastic_usage()
{
}

/*******************************************************************************
 * main program for seismic wave simulation
 ******************************************************************************/

int main(int argc, char** argv){
{

    char *pname = "cgfd3d-elastic";
    char par_file_name[FD_MAX_STRLEN];

//  Declare variable
    struct fd_t  *fd;
    struct blk_t *blk;  // wavefiled on each block, use pointer to make main/function sytax same
    struct blk_conn_t *blk_conn;  // connection of each blocks
    struct par_t *par;

    int myid,mpi_size;

//-------------------------------------------------------------------------------
// Program Entries:
//-------------------------------------------------------------------------------

//
// init MPI
//
    fprintf(stdout,"init MPI\n"); 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (myid==0) fprintf(stdout,"MPI_COMM_WORD=%d\n", MPI_COMM_WORD); 
    if (myid==0) fprintf(stdout,"comm size=%d\n", mpi_size); 

//
// check commond line
//  todo: check if slave thread can get arguments
//
    if (myid==0) fprintf(stdout,"check commond line arguments ...\n"); 
    if (argc < 1) {
      cgfd3d_usage()
    }

    // get par file name
    if (myid==0) fprintf(stdout,"get input par file name ...\n"); 
    strcpy(par_file_name, argv[1]);
    if (myid==0) fprintf(stdout,"   par file name =  %s\n", par_file_name); 

    if (is_file_not_exist(par_file_name)) {
        exit(-1);
    }

//
// alloc
//
    fd  = (struct fd_t  *) malloc(sizeof(struct fd_t ));
    blk = (struct blk_t *) malloc(sizeof(struct blk_t));
    par = (struct par_t *) malloc(sizeof(struct par_t));

//
// set up fd operator, which determins number of ghost points
//
//  do not support selection scheme by par file right now
    if (myid==0) fprintf(stdout,"set scheme ...\n"); 
    ierr = fd_set_macdrp(fd);
    /*
    if      (par->fd_scheme == PAR_FD_SCHEME_MACDRP)
    {
    }
    else if (par->fd_scheme == FD_SCHEME_CENT_6th)
    {
        //ierr = scheme_macdrp();
    }
    else // read in from par file
    {
    }
    */

// 
// read in pars on master node and exchange to others
// 
    par = (struct par_struct *) malloc(sizeof(par_struct));
    if (myid==0) {
        fprintf(stdout,"readin par file on master node ...\n"); 
        ierr = par_read_file(par_file_name, &par);
    }

    if (myid==0) fprintf(stdout,"exchange par to all nodes ...\n"); 
    ierr = par_exchange(&par);

//
// create mpi topo
//  todo: unqual partition between nodes
//
    if (myid==0) fprintf(stdout,"create MPI topo ...\n"); 
    ierr = fdmpi_topo_creat(&blk_connect, par->num_threads_per_dim);

//
// set parmeters of domain size
//
    if (myid==0) fprintf(stdout,"set slocal/global grid parameters ...\n"); 
    ierr = blk_struct_set_grid_size(&blk, par->nx, par->ny, par->nz);

    // set how many levels per var due to time int scheme
    ierr = blk_struct_set_wave_level(&blk, scheme->number_of_var_levels);

//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

    // allocate grid array
    if (myid==0) fprintf(stdout,"allocate grid vars ...\n"); 
    ierr = gd_curv_init_vars(blk->siz_volume, &(blk->g3d_num_of_vars),
                &(blk->g3d), &(blk->g3d_pos), &(blk->g3d_name));

    // if import coord and metric
    if (par->grid_by_import==1)
    {
        if (myid==0) fprintf(stdout,"import grid vars ...\n"); 

        ierr = gd_curv_import(blk->g3d, blk->g3d_pos, blk->g3d_name,
                blk->g3d_num_of_vars, blk->siz_volume, par->in_metric_dir, myid3);

        ierr = gd_curv_topoall_import(blk->g3d, blk->nx, blk->ny, blk->nz, 
                    myid3[0],myid3[1],myid3[2]);

    }
    //// if cartesian coord
    //else if (par->grid_by_cartesian==1)
    //{
    //    ierr = gd_curv_generate_cartesian(blk->g3d, blk->nx, blk->ny, blk->nz,
    //                par->grid_x0 + myid3[0]*blk->ni - blk->num_of_ghost[0],
    //                par->grid_y0 + myid3[1]*blk->ni - blk->num_of_ghost[1],
    //                par->grid_z0 + myid3[2]*blk->ni - blk->num_of_ghost[2],
    //                par->grid_dx, par->grid_dy, par->grid_dz);
    //}
    // if vmap coord
    else if (par->grid_by_vmap==1)
    {
        ierr = gd_curv_generate_vmap(blk->g3d, blk->nx, blk->ny, blk->nz,
                );
    }
    // if refine coord
    else if (par->grid_by_refine==1)
    {
        grid_coarse = gd_curv_read_coarse(par->grid_in_file, par%grid_dz_taper_length);

        ierr = gd_curv_interp_corse(grid_corse, blk->g3d,
                    );

        gd_curv_coarse_free(grid_coarse);
    }
    else
    {
    }

    // cal metrics and output for QC
    if (par->grid_by_import==0)
    {
        // generate topo over all the domain
        ierr = gd_curv_topoall_generate();

        // output
        if (par->grid_output_to_file==1) {
            ierr = gd_curv_coord_export();
        }

        if (myid==0) fprintf(stdout,"calculate metrics ...\n"); 
        ierr = gd_curv_cal_metric();

        if (myid==0) fprintf(stdout,"exchange metrics ...\n"); 
        ierr = gd_curv_exchange_metric();

        if (par->metric_output_to_file==1) {
            if (myid==0) fprintf(stdout,"export metric to file ...\n"); 
            ierr = gd_curv_metric_export();
        }
    }

//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

    // allocate media vars
    if (myid==0) fprintf(stdout,"allocate media vars ...\n"); 
    ierr = md_el_iso_init_vars(blk->siz_volume, &(blk->m3d_num_of_vars),
                &(blk->m3d_pos), &(blk->m3d_name));

    if (par->medium_by_import==1)
    {
        ierr = md_el_iso_import(blk->m3d, blk->nx, blk->ny, blk->nz, 
                    myid3[0],myid3[1],myid3[2]);
    }
    else
    {

        if (par->medium_output_to_file==1) {
            if (myid==0) fprintf(stdout,"export discrete medium to file ...\n"); 
            ierr = md_el_iso_export();
        }
    }

//-------------------------------------------------------------------------------
//-- estimate/check time step
//-------------------------------------------------------------------------------

    if (par>check_stability==1) {

        //-- estimate time step
        if (myid==0) fprintf(stdout,"   estimate time step ...\n"); 
        ierr = curv_macdrp_elastic_iso_stept_calculate(blk->g3d, blk->m3d,
                blk->nx, blk->ny, blk->nz,
                &dtmax, &dtmaxVp, &dtmaxL, &dtindx)
        
        //-- print for QC
        fprintf(stdout, "-> thisid=%d,%d,%d, dtmax=%f, dtmaxVp=%f, dtmaxL=%f, dtindx=%f\n",
                myid3[0],myid3[1],myid3[2],
                dtmax, dtmaxVp, dtmaxL, dtindx);
        
        dtin[0] = dtmax; dtin[1] = myid;
        ierr = MPI_Reduce(dtin,dtout,1,MPI_2REAL,MPI_MINLOC, 0,MPI_COMM_WORLD);

        dtmax=dtout[1];
        if (myid==0) {
           print *, "Maximum allowed time step is", dtmax, ' in thread', nint(dtout(2))
            fprintf(stdout,"Maximum allowed time step is %f  in thread %d\n", dtmax, (int)dtout[1]);
           //-- auto set stept
           if (stept < 0.0) {
              stept = dtmax;
              fprintf(stdout, "-> Set stept to maximum allowed value\n");
           }
           //-- check value
           if (dtmax<stept) {
              fprint(stdout, "Serious Error: stept %f >dtmax %f\n", stept,dtmax);
              fprint(stdout, " occurs on thread", int(dtout(2)));
              MPI_Abort(MPI_COMM_WORLD,1);
           }
        }
        
        //-- from root to all threads
        ierr = MPI_Bcast(&stept, 1, MPI_REAL, 0, MPI_COMM_WORLD);
    }

//-------------------------------------------------------------------------------
//-- source import or locate on fly
//-------------------------------------------------------------------------------

    if (par->source_by_import)
    {
        if (myid==0) fprintf(stdout,"import source vars ...\n"); 
        ierr = src_import();
    } 
    else
    {
        ierr = src_read_infile(par->in_src_file, myid);
        ierr = src_locate_extend();

        if (par->source_output_to_file==1)
        {
            ierr = src_export();
        }
    }

//-------------------------------------------------------------------------------
//-- setup output
//-------------------------------------------------------------------------------
    if (myid==0) fprintf(stdout,"setup output info ...\n"); 

    // receiver

    // inline

    // snapshot

//-------------------------------------------------------------------------------
//-- allocate main var
//-------------------------------------------------------------------------------

    if (myid==0) fprintf(stdout,"allocate solver vars ...\n"); 
    ierr = wf_el_1st_init_vars(blk->siz_volume, blk->number_of_levels,
          &blk->w3d_num_of_vars, &blk->w3d, &blk->w3d_pos, &blk->w3d_name);

//-------------------------------------------------------------------------------
//-- absorbing boundary etc auxiliary variables
//-------------------------------------------------------------------------------

    if (myid==0) fprintf(stdout,"setup absorbingg boundary ...\n"); 

    if (par->absorbing_boundary_exp==1) {
        ierr = abs_exp_cal_damping();
    }

    if (par->absorbing_boundary_adecfspml==1) {
      // set pml parameters
        ierr = abs_set_cfspml(blk->ni, blk->nj, blk->nk,
            blk->boundary_itype,
            par->abs_num_of_layers,
            par->abs_alpha, par->abs_beta, par->abs_velocity,
            blk->abs_num_of_layers, blk->abs_coefs_dimpos, &blk->abs_coefs
            );

        // alloc pml vars
        ierr = wf_el_1st_init_cfspml_vars(blk->number_of_levels, blk->number_of_vars,
            blk->boundary_itype, blk->abs_blk_indx, blk->abs_blk_vars_siz_volume, &blk->abs_blk_vars
            );
    }

//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------

    if (myid==0) fprintf(stdout,"start time loop ...\n"); 

    t_start = time(NULL);

    ierr = sv_eliso1st_curv_macdrp_allstep(blk->w3d, blk->g3d, blk->m3d,
            blk->nx, blk->ny, blk->nz,
            fd->pair_fdx_all_info, fd->pair_fdx_all_indx, fd->pair_fdx_all_coef,
            fd->pair_fdy_all_info, fd->pair_fdy_all_indx, fd->pair_fdy_all_coef,
            fd->pair_fdz_all_info, fd->pair_fdz_all_indx, fd->pair_fdz_all_coef,
            );

    t_end = time(NULL);

    printf("\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));

//-------------------------------------------------------------------------------
//-- postprocess
//-------------------------------------------------------------------------------

    return 0;
}
