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

#include "fd_t.h"
#include "par_t.h"
#include "gd_curv.h"
#include "md_el_iso.h"
#include "wf_el_1st.h"
#include "abs_funcs.h"
#include "src_funcs.h"
#include "io_funcs.h"
#include "sv_eliso1st_curv_macdrp.h"

int main(int argc, char** argv)
{
  int verbose = 1; // default fprint
  char *par_fname;
  char ou_file[FD_MAX_STRLEN];

//-------------------------------------------------------------------------------
// start MPI and read par
//-------------------------------------------------------------------------------

  // init MPI

  int myid, mpi_size;
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

    //strncpy(par_fname, argv[1], sizeof(argv[1]));
    par_fname = argv[1];

    if (argc >= 3) {
      verbose = atoi(argv[2]); // verbose number
      fprintf(stdout,"verbose=%d\n", verbose); fflush(stdout);
    }

    // bcast verbose to all nodes
    MPI_Bcast(&verbose, 1, MPI_INT, 0, comm);
  }
  else
  {
    // get verbose from id 0
    MPI_Bcast(&verbose, 1, MPI_INT, 0, comm);
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
  float   dt = par->size_of_time_step;
  int     nt_total = par->number_of_time_steps+1;

  // alloc
  struct fd_t *fd = (struct fd_t  *) malloc(sizeof(struct fd_t ));

  //  do not support selection scheme by par file right now
  if (myid==0 && verbose>0) fprintf(stdout,"set scheme ...\n"); 

  fd_set_macdrp(fd);

//-------------------------------------------------------------------------------
// init blk
//-------------------------------------------------------------------------------

  struct fd_blk_t *blk = (struct fd_blk_t *) malloc(sizeof(struct fd_blk_t));

  // set parmeters of domain size
  if (myid==0 && verbose>0) fprintf(stdout,"set slocal/global grid parameters ...\n"); 

  fd_blk_init(blk,
              par->grid_name,
              par->number_of_total_grid_points_x,
              par->number_of_total_grid_points_y,
              par->number_of_total_grid_points_z,
              par->number_of_mpiprocs_x,
              par->number_of_mpiprocs_y,
              par->boundary_type_name,
              par->abs_num_of_layers,
              par->output_dir,
              fd->fdx_nghosts,
              fd->fdy_nghosts,
              fd->fdz_nghosts,
              fd->num_rk_stages,
              comm,
              myid, verbose);

//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

  // allocate grid array
  if (myid==0 && verbose>0) fprintf(stdout,"allocate grid vars ...\n"); 

  gd_curv_init_c3d( blk->siz_volume,
                   &blk->c3d_num_of_vars,
                   &blk->c3d,
                   &blk->c3d_pos,
                   &blk->c3d_name,
                   &blk->coord_name);

  gd_curv_init_g3d( blk->siz_volume,
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
      gd_curv_gen_cart(blk->c3d, blk->siz_volume,
                       blk->nx,
                       par->cartesian_grid_dx,
                       par->cartesian_grid_x0 + (blk->gni1 - fd->fdx_nghosts) * par->cartesian_grid_dx,
                       blk->ny,
                       par->cartesian_grid_dy,
                       par->cartesian_grid_y0 + (blk->gnj1 - fd->fdy_nghosts) * par->cartesian_grid_dy,
                       blk->nz,
                       par->cartesian_grid_dz,
                       par->cartesian_grid_z0 + (blk->gnk1 - fd->fdz_nghosts) * par->cartesian_grid_dz);
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

  // generate topo over all the domain
  //ierr = gd_curv_topoall_generate();

  // output
  //if (par->grid_output_to_file==1) {
      io_build_fname(blk->output_dir,"coord",".nc",blk->myid2,ou_file);
      io_var3d_export_nc(ou_file,
                         blk->c3d,
                         blk->c3d_pos,
                         blk->c3d_name,
                         blk->c3d_num_of_vars,
                         blk->coord_name,
                         blk->nx,
                         blk->ny,
                         blk->nz);
  //}

  // cal metrics and output for QC
  if (par->metric_by_import==0)
  {
    if (myid==0 && verbose>0) fprintf(stdout,"calculate metrics ...\n"); 

    gd_curv_cal_metric(blk->c3d,
                       blk->g3d,
                       blk->ni1,
                       blk->ni2,
                       blk->nj1,
                       blk->nj2,
                       blk->nk1,
                       blk->nk2,
                       blk->nx,
                       blk->ny,
                       blk->nz,
                       blk->siz_line,
                       blk->siz_slice,
                       blk->siz_volume,
                       fd->fd_len,
                       fd->fd_indx,
                       fd->fd_coef);

    //if (myid==0) fprintf(stdout,"exchange metrics ...\n"); 
    //gd_curv_exchange_metric();

    //if (par->metric_output_to_file==1) {
        if (myid==0) fprintf(stdout,"export metric to file ...\n"); 
        io_build_fname(blk->output_dir,"metric",".nc",blk->myid2,ou_file);
        io_var3d_export_nc(ou_file,
                           blk->g3d,
                           blk->g3d_pos,
                           blk->g3d_name,
                           blk->g3d_num_of_vars,
                           blk->coord_name,
                           blk->nx,
                           blk->ny,
                           blk->nz);
    //}
  }
  else // import metric
  {
  }

//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

  // allocate media vars
  if (myid==0 && verbose>0) fprintf(stdout,"allocate media vars ...\n"); 
  md_el_iso_init_vars( blk->siz_volume,
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
    md_el_iso_gen_test(blk->m3d,
                       blk->c3d+GD_CURV_SEQ_X3D * blk->siz_volume,
                       blk->c3d+GD_CURV_SEQ_Y3D * blk->siz_volume,
                       blk->c3d+GD_CURV_SEQ_Z3D * blk->siz_volume,
                       blk->nx,
                       blk->ny,
                       blk->nz,
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
  md_el_iso_rho_to_slow(blk->m3d, blk->siz_volume);

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
  //src_gen_test(t0,
  //             dt,
  //             nt_total,
  //             fd->num_rk_stages,
  //             &blk->num_of_force,
  //             &blk->force_info,
  //             &blk->force_vec_stf,
  //             &blk->num_of_moment,
  //             &blk->moment_info,
  //             &blk->moment_ten_rate,
  //             verbose);

  src_gen_test_gauss(blk->siz_line,
                     blk->siz_slice,
                     t0,
                     dt,
                     nt_total,
                     fd->num_rk_stages,
                     &blk->num_of_force,
                     &blk->force_info,
                     &blk->force_vec_stf,
                     &blk->num_of_moment,
                     &blk->moment_info,
                     &blk->moment_ten_rate,
                     &blk->moment_ext_indx,
                     &blk->moment_ext_coef,
                     verbose);

  // manually del source except 0
  //if (myid!=0) {
  //  blk->num_of_moment = 0;
  //}
  
  /*
  if (par->source_output_to_file==1)
  {
      ierr = src_export();
  }
  */

//-------------------------------------------------------------------------------
//-- setup output, may require coord info
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"setup output info ...\n"); 
  
  // receiver
  
  // inline
  
  // snapshot
  fd_blk_set_snapshot(blk,
                      fd->fdx_nghosts,
                      par->number_of_snapshot,
                      par->snapshot_index_start,
                      par->snapshot_index_count,
                      par->snapshot_index_stride,
                      par->snapshot_time_start,
                      par->snapshot_time_count,
                      par->snapshot_time_stride);

//-------------------------------------------------------------------------------
//-- allocate main var
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"allocate solver vars ...\n"); 
  wf_el_1st_init_vars(blk->siz_volume,
                      blk->w3d_num_of_levels,
                      &blk->w3d_num_of_vars,
                      &blk->w3d,
                      &blk->w3d_pos,
                      &blk->w3d_name);

//-------------------------------------------------------------------------------
//-- absorbing boundary etc auxiliary variables
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"setup absorbingg boundary ...\n"); 
  
  /*
  if (par->absorbing_boundary_exp==1) {
      ierr = abs_exp_cal_damping();
  }
  */
  
  if (blk->abs_itype == FD_BOUNDARY_TYPE_CFSPML) {
    // set pml parameters
    abs_set_cfspml(par->cfspml_alpha_max,
                   par->cfspml_beta_max,
                   par->cfspml_velocity,
                   par->cartesian_grid_dx,
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
    abs_init_vars_cfspml(blk->w3d_num_of_levels,
                         blk->w3d_num_of_vars,
                         blk->boundary_itype, 
                         blk->abs_indx, 
                         blk->abs_vars_volsiz, 
                         blk->abs_vars_facepos0,
                         &blk->abs_vars_size_per_level,
                         &blk->abs_vars,
                         myid, verbose);
  }

//-------------------------------------------------------------------------------
//-- setup mesg
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"init mesg ...\n"); 
  fd_blk_init_mpi_mesg(blk,
                       fd->fdx_nghosts,
                       fd->fdy_nghosts);

//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------

  // boundary preproc
  if (myid==0 && verbose>0) fprintf(stdout,"cal free surface matrix ...\n"); 

  if (blk->boundary_itype[5] == FD_BOUNDARY_TYPE_FREE)
  {
    sv_eliso1st_curv_macdrp_vel_dxy2dz(blk->g3d,
                                       blk->m3d,
                                       blk->ni1,
                                       blk->ni2,
                                       blk->nj1,
                                       blk->nj2,
                                       blk->nk1,
                                       blk->nk2,
                                       blk->siz_line,
                                       blk->siz_slice,
                                       blk->siz_volume,
                                       &blk->matVx2Vz,
                                       &blk->matVy2Vz,
                                       myid, verbose);
  }

  if (myid==0 && verbose>0) fprintf(stdout,"start time loop ...\n"); 
  
  time_t t_start = time(NULL);
  
  sv_eliso1st_curv_macdrp_allstep(blk->w3d,blk->w3d_pos,blk->w3d_name,blk->w3d_num_of_vars,
                                  blk->coord_name,
                                  blk->g3d, blk->m3d,
                                  blk->ni1,blk->ni2,blk->nj1,blk->nj2,blk->nk1,blk->nk2,
                                  blk->ni,blk->nj,blk->nk,blk->nx,blk->ny,blk->nz,
                                  blk->siz_line, blk->siz_slice, blk->siz_volume,
                                  blk->boundary_itype, blk->abs_itype,
                                  blk->abs_num_of_layers, blk->abs_indx,
                                  blk->abs_coefs_facepos0, blk->abs_coefs,
                                  blk->abs_vars_size_per_level,
                                  blk->abs_vars_volsiz, blk->abs_vars_facepos0, blk->abs_vars,
                                  blk->matVx2Vz, blk->matVy2Vz,
                                  blk->num_of_force, blk->force_info, blk->force_vec_stf,
                                  blk->num_of_moment, blk->moment_info, blk->moment_ten_rate,
                                  blk->moment_ext_indx,blk->moment_ext_coef,
                                  blk->num_of_sta, blk->sta_loc_point, blk->sta_seismo,
                                  blk->num_of_snap, blk->snap_grid_indx, blk->snap_time_indx,
                                  fd->num_rk_stages, fd->rk_a, fd->rk_b,
                                  fd->num_of_pairs,
                                  fd->fdx_max_half_len,fd->fdy_max_half_len,
                                  fd->fdz_max_len,fd->fdz_num_surf_lay,
                                  fd->pair_fdx_all_info, fd->pair_fdx_all_indx, fd->pair_fdx_all_coef,
                                  fd->pair_fdy_all_info, fd->pair_fdy_all_indx, fd->pair_fdy_all_coef,
                                  fd->pair_fdz_all_info, fd->pair_fdz_all_indx, fd->pair_fdz_all_coef,
                                  dt,nt_total,t0,
                                  myid, blk->myid2, comm,
                                  blk->sbuff, blk->rbuff, blk->s_reqs, blk->r_reqs,
                                  par->check_nan_every_nummber_of_steps,
                                  verbose,
                                  blk->name,
                                  blk->output_dir);
  
  time_t t_end = time(NULL);
  
  fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));

//-------------------------------------------------------------------------------
//-- postprocess
//-------------------------------------------------------------------------------

  return 0;
}
