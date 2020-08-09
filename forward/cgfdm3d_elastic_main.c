/*******************************************************************************
 * Curvilinear Grid Finite Difference Seismic Wave Propagation Simulation 
 *
 * Copyright (c) 2020 ZHANG Wei. All rights reserved.
 *
 * Author(s): ZHANG Wei <zhangwei@sustech.edu.cn>
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <mpi.h>

#include "fdlib_mem.h"
#include "fd_t.h"
#include "par_t.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "sv_eliso1st_curv_macdrp.h"

int main(int argc, char** argv)
{
  int verbose = 1; // default fprint
  char par_fname[FD_MAX_STRLEN];

//-------------------------------------------------------------------------------
// start MPI and read par
//-------------------------------------------------------------------------------

  // init MPI

  int myid,mpi_size;
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &myid);
  MPI_Comm_size(comm, &mpi_size);

  // get commond-line argument

  if (myid==0) 
  {
    // argc checking
    if (argc < 2) {
      fprintf(stdout,"usage: cgfdm3d_elastic <par_file> <opt: verbose>\n");
      MPI_Finalize();
      exit(1);
    }

    strncpy(par_fname, argv[1], sizeof(argv[1]));

    if (argc >= 3) {
      verbose = atoi(argv[2]); // verbose number
    }

    // bcast verbose to all nodes
    MPI_Bcast(verbose, 1, MPI_INT, 0, comm);
  }
  else
  {
    // get verbose from id 0
    MPI_Bcast(verbose, 1, MPI_INT, 0, comm);
  }

  if (myid==0 && verbose>0) fprintf(stdout,"comm=%d, size=%d\n", comm, mpi_size); 
  if (myid==0 && verbose>0) fprintf(stdout,"par file =  %s\n", par_fname); 

  // read par

  struct par_t *par = (struct par_t *) malloc(sizeof(struct par_t));

  par_mpi_get(par_fname, myid, comm, par, verbose);

  if (myid==0 && verbose>0) par_print(par);

//-------------------------------------------------------------------------------
// set up based on input info
//-------------------------------------------------------------------------------

  // time
  float   t0 = par->time_start;
  float   dt = par->time_step;
  int     nt_total = (int) ((par->time_end - par->time_start) / dt+0.5);

  // create mpi topo
  //  todo: unqual partition between nodes
  fdmpi = (struct fd_mpi_t *) malloc(sizeof(struct fd_mpi_t));

  fd_mpi_create_topo(fdmpi, par->num_threads_per_dim);

  // alloc
  struct fd_t *fd = (struct fd_t  *) malloc(sizeof(struct fd_t ));

  //  do not support selection scheme by par file right now
  if (myid==0 && verbose>0) fprintf(stdout,"set scheme ...\n"); 

  fd_set_macdrp(fd);

//-------------------------------------------------------------------------------
// init blk
//-------------------------------------------------------------------------------

  struct fd_blk_t *blk = (struct blk_t *) malloc(sizeof(struct blk_t));

  // set parmeters of domain size
  if (myid==0) fprintf(stdout,"set slocal/global grid parameters ...\n"); 

  fd_blk_init(blk,
      par->number_of_x_points,
      par->number_of_y_points,
      par->number_of_z_points,
      par->number_of_x_procs,
      par->number_of_y_procs,
      par->number_of_z_procs,
      par->boundary_type_name,
      par->abs_type_name,
      par->abs_number_of_layers,
      fd->fdx_nghosts,
      fd->fdy_nghosts,
      fd->fdz_nghosts,
      fd->num_rk_stages,
      fdmpi->myid2,
      fdmpi->neighid,
      myid, verbose);

//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

  // allocate grid array
  if (myid==0) fprintf(stdout,"allocate grid vars ...\n"); 

  gd_curv_init_c3d(
      blk->siz_volume,
      &blk->c3d_num_of_vars,
      &blk->c3d,
      &blk->c3d_pos,
      &blk->c3d_name);

  gd_curv_init_g3d(
      blk->siz_volume,
      &blk->g3d_num_of_vars,
      &blk->g3d,
      &blk->g3d_pos,
      &blk->g3d_name);

  // if import coord and metric
  /*
  if (par->coord_by_import==1)
  {
      if (myid==0) fprintf(stdout,"import grid vars ...\n"); 

      gd_curv_import(blk->g3d, blk->g3d_pos, blk->g3d_name,
              blk->g3d_num_of_vars, blk->siz_volume, par->in_metric_dir, myid3);

      gd_curv_topoall_import(blk->g3d, blk->nx, blk->ny, blk->nz, 
                  myid3[0],myid3[1],myid3[2]);

  }
  */
  // if cartesian coord
  if (par->coord_by_cartesian==1)
  {
      gd_curv_gen_cart(blk->c3d, blk->nx, blk->ny, blk->nz, blk->siz_volume,
          blk->nx, par->grid_dx, par->grid_x0 + (blk->gni - fd->fdx_nghosts) * par->grid_dx,
          blk->ny, par->grid_dy, par->grid_y0 + (blk->gnj - fd->fdy_nghosts) * par->grid_dy,
          blk->nz, par->grid_dz, par->grid_z0 + (blk->gnk - fd->fdz_nghosts) * par->grid_dz);
  }
  /*
  // if vmap coord
  if (par->coord_by_vmap==1)
  {
      ierr = gd_curv_gen_vmap(blk->g3d, blk->nx, blk->ny, blk->nz,
              );
  }
  // if refine coord
  if (par->coord_by_refine==1)
  {
      grid_coarse = gd_curv_read_coarse(par->grid_in_file, par%grid_dz_taper_length);

      ierr = gd_curv_interp_corse(grid_corse, blk->g3d,
                  );

      gd_curv_coarse_free(grid_coarse);
  }
  else
  {
  }
  */

  /*
  // generate topo over all the domain
  ierr = gd_curv_topoall_generate();

  // output
  if (par->grid_output_to_file==1) {
      ierr = gd_curv_coord_export();
  }
  */

  // cal metrics and output for QC
  if (par->metric_by_import==0)
  {
    if (myid==0 && verbose>0) fprintf(stdout,"calculate metrics ...\n"); 
    size_t fd_len   = fd->mac_center_all_info[3][1];
    size_t fd_pos   = fd->mac_center_all_info[3][0];
    size_t *fd_indx = fd->mac_center_all_indx+fd_pos;
    float  *fd_coef = fd->mac_center_all_coef+fd_pos;

    gd_curv_cal_metric(
        blk->c3d,
        blk->g3d,
        blk->ni1,
        blk->ni2,
        blk->nj1,
        blk->nj2,
        blk->nk1,
        blk->nk2,
        blk->siz_line,
        blk->siz_slice,
        blk->siz_volume,
        fd_len,
        fd_indx,
        fd_coef
        );

    /*
    if (myid==0) fprintf(stdout,"exchange metrics ...\n"); 
    gd_curv_exchange_metric();

    if (par->metric_output_to_file==1) {
        if (myid==0) fprintf(stdout,"export metric to file ...\n"); 
        gd_curv_metric_export();
    }
    */
  }
  else // import metric
  {
  }

//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

  // allocate media vars
  if (myid==0) fprintf(stdout,"allocate media vars ...\n"); 
  md_el_iso_init_vars(
      blk->siz_volume,
      &blk->m3d_num_of_vars,
      &blk->m3d,
      &blk->m3d_pos,
      &blk->m3d_name);
  
  if (par->medium_by_import==1)
  {
    //md_el_iso_import(blk->m3d, blk->nx, blk->ny, blk->nz, 
    //              myid3[0],myid3[1],myid3[2]);
  }
  else
  {
    md_el_iso_gen_test(
        blk->m3d,
        blk->c3d+GD_CURV_SEQ_X3D * blk->siz_volume,
        blk->c3d+GD_CURV_SEQ_Y3D * blk->siz_volume,
        blk->c3d+GD_CURV_SEQ_Z3D * blk->siz_volume,
        blk->nx,
        blk->ny
        blk->nz
        blk->siz_line,
        blk->siz_slice,
        blk->siz_volume);

    /*
    if (par->medium_output_to_file==1) {
      if (myid==0) fprintf(stdout,"export discrete medium to file ...\n"); 
      md_el_iso_export();
    }
    */
  }
  
  // convert rho to 1 / rho to reduce number of arithmetic cal
  md_el_iso_rho_to_slow(
      blk->m3d,
      blk->nx,
      blk->ny,
      blk->nz,
      blk->siz_line,
      blk->siz_slice,
      blk->siz_volume);

//-------------------------------------------------------------------------------
//-- estimate/check time step
//-------------------------------------------------------------------------------

  /*
  if (par>check_stability==1)
  {
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
  */

//-------------------------------------------------------------------------------
//-- source import or locate on fly
//-------------------------------------------------------------------------------

  //src_read_file(par->src_fname, myid);
  //src_locate_extend();
  src_gen_test(
      t0,
      dt,
      nt_total,
      fd->num_rk_stages,
      blk->num_of_force,
      blk->force_info,
      blk->force_vec_stf,
      blk->num_of_moment,
      blk->moment_loc_point,
      blk->moment_ten_rate,
      verbose);
  
  /*
  if (par->source_output_to_file==1)
  {
      ierr = src_export();
  }
  */

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
  wf_el_1st_init_vars(
      blk->siz_volume,
      blk->number_of_levels,
      &blk->w3d_num_of_vars,
      &blk->w3d,
      &blk->w3d_pos,
      &blk->w3d_name);

//-------------------------------------------------------------------------------
//-- absorbing boundary etc auxiliary variables
//-------------------------------------------------------------------------------

  if (myid==0) fprintf(stdout,"setup absorbingg boundary ...\n"); 
  
  /*
  if (par->absorbing_boundary_exp==1) {
      ierr = abs_exp_cal_damping();
  }
  */
  
  if (blk->abs_itype == FD_BOUNDARY_TYPE_CFSPML) {
    // set pml parameters
    abs_set_cfspml(
        par->cfs_pml_alpha_max,
        par->cfs_pml_beta_max,
        par->cfs_pml_velocity,
        blk->ni1,
        blk->ni2,
        blk->nj1,
        blk->nj2,
        blk->nk1,
        blk->nk2,
        blk->boundary_itype,
        blk->abs_num_of_layers,
        blk->abs_indx,
        blk->abs_coefs_facepos0,
        &(blk->abs_coefs),
        verbose);
  
    // alloc pml vars
    abs_init_vars_cfspml(
        blk->number_of_levels,
        blk->number_of_vars,
        blk->boundary_itype, 
        blk->abs_indx, 
        blk->abs_vars_volsiz, 
        blk->abs_vars_facepos0,
        &blk->abs_vars_size_per_level,
        &blk->abs_blk_vars,
        myid, verbose);
  }

//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------

  // boundary preproc
  if (myid==0 && verbose>0) fprintf(stdout,"cal free surface matrix ...\n"); 

  if (blk->boundary_itype[5] == FD_BOUNDARY_TYPE_FREE)
  {
    sv_eliso1st_curv_macdrp_vel_dxy2dz(
          blk->g3d,
          blk->m3d,
          blk->ni1,
          blk->ni2,
          blk->nj1,
          blk->nj2,
          blk->nk1,
          blk->nk2,
          blk->siz_line,
          blk->siz_volume,
          &blk->matVx2Vz,
          &blk->matVy2Vz,
          myid, verbose);
  }

  if (myid==0) fprintf(stdout,"start time loop ...\n"); 
  
  t_start = time(NULL);
  
  sv_eliso1st_curv_macdrp_allstep(
            blk->w3d,
            blk->g3d,
            blk->m3d,
            blk->ni1,blk->ni2,blk->nj1,blk->nj2,blk->nk1,blk->nk2,
            blk->ni,blk->nj,blk->nk,blk->nx,blk->ny,blk->nz,
            blk->siz_line,
            blk->siz_slice,
            blk->siz_volume,
            blk->boundary_ityp,
            blk->abs_itype,
            blk->abs_num_of_layers,
            blk->abs_indx,
            blk->abs_coefs_facepos0,
            blk->abs_coefs,
            blk->abs_vars_volsiz,
            blk->abs_vars_facepos0,
            blk->abs_vars,
            blk->matVx2Vz,
            blk->matVy2Vz,
            blk->num_of_force,
            blk->force_loc_point,
            blk->force_vec_stf,
            blk->num_of_moment,
            blk->moment_loc_point,
            blk->moment_ten_rate,
            blk->num_of_sta,
            blk->sta_loc_point,
            blk->sta_seismo,
            blk->num_of_snap,
            blk->snap_grid_indx,
            blk->snap_time_indx,
            fd->num_rk_stages,
            fd->rk_a,
            fd->rk_b,
            fd->num_of_pair,
            fd->fdx_max_len,fd->fdy_max_len,fd->fdz_max_len,fd->fdz_num_surf_lay,
            fd->pair_fdx_all_info, fd->pair_fdx_all_indx, fd->pair_fdx_all_coef,
            fd->pair_fdy_all_info, fd->pair_fdy_all_indx, fd->pair_fdy_all_coef,
            fd->pair_fdz_all_info, fd->pair_fdz_all_indx, fd->pair_fdz_all_coef,
            dt,nt_total,t0,
            myid, myid2, comm,
            par->qc_check_nan_num_of_step,
            verbose,
            out_dir);
  
  t_end = time(NULL);
  
  fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));

//-------------------------------------------------------------------------------
//-- postprocess
//-------------------------------------------------------------------------------

  return 0;
}
