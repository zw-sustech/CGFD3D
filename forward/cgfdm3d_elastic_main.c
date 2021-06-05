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

// in lib
#include "sacLib.h"

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
  char err_message[FD_MAX_STRLEN];

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
  if (myid==0 && verbose>0) fprintf(stdout,"set local/global grid parameters ...\n"); 

  fd_blk_init(blk,
              par->number_of_total_grid_points_x,
              par->number_of_total_grid_points_y,
              par->number_of_total_grid_points_z,
              par->number_of_mpiprocs_x,
              par->number_of_mpiprocs_y,
              par->boundary_type_name,
              par->abs_num_of_layers,
              par->output_dir,
              par->grid_export_dir,
              par->media_export_dir,
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

  switch (par->grid_generation_itype)
  {
    case PAR_GRID_CARTESIAN :
        if (myid==0) fprintf(stdout,"generate cartesian grid in code ...\n"); 
        float x0 = par->cartesian_grid_origin[0];
        float y0 = par->cartesian_grid_origin[1];
        float z0 = par->cartesian_grid_origin[2];
        float dx = par->cartesian_grid_stepsize[0];
        float dy = par->cartesian_grid_stepsize[1];
        float dz = par->cartesian_grid_stepsize[2];
        gd_curv_gen_cart(blk->c3d, blk->siz_volume,
                         blk->nx,
                         dx,
                         x0 + (blk->gni1 - fd->fdx_nghosts) * dx,
                         blk->ny,
                         dy,
                         y0 + (blk->gnj1 - fd->fdy_nghosts) * dy,
                         blk->nz,
                         dz,
                         z0 + (blk->gnk1 - fd->fdz_nghosts) * dz);
        break;

    case PAR_GRID_IMPORT :
        if (myid==0) fprintf(stdout,"import grid vars ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 
        //gd_curv_import(blk->g3d, blk->g3d_pos, blk->g3d_name,
        //        blk->g3d_num_of_vars, blk->siz_volume, par->in_metric_dir, myid3);
        //gd_curv_topoall_import(blk->g3d, blk->nx, blk->ny, blk->nz, 
        //          myid3[0],myid3[1],myid3[2]);
        break;

    case PAR_GRID_LAYER_INTERP :
        if (myid==0) fprintf(stdout,"gerate grid using layer interp ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 

        int VmapSpacingIsequal[] = { 1, 0, 1 }; // 1:The grid spacing is equal; 0: Is not.
        int NCellPerlay[] = { 25, 15, 20 };
        int nLayers = 3;
        gd_curv_gen_layer(blk->c3d, nLayers,
                          NCellPerlay, VmapSpacingIsequal,
                          blk->nx, blk->ni, blk->gni1, fd->fdx_nghosts, par->number_of_total_grid_points_x,
                          blk->ny, blk->nj, blk->gnj1, fd->fdy_nghosts, par->number_of_total_grid_points_y,
                          blk->nz, blk->nk, blk->gnk1, fd->fdz_nghosts, par->number_of_total_grid_points_z);

        break;

  }

  // generate topo over all the domain
  //ierr = gd_curv_topoall_generate();

  // output
  if (par->is_export_grid==1)
  {
    if (myid==0) fprintf(stdout,"export coord to file ...\n"); 
    gd_curv_coord_export(blk->c3d,
                         blk->c3d_pos,
                         blk->c3d_name,
                         blk->c3d_num_of_vars,
                         blk->nx,
                         blk->ny,
                         blk->nz,
                         blk->ni1,
                         blk->nj1,
                         blk->nk1,
                         blk->ni,
                         blk->nj,
                         blk->nk,
                         blk->gni1,
                         blk->gnj1,
                         blk->gnk1,
                         blk->output_fname_part,
                         blk->grid_export_dir);
  }

  // cal metrics and output for QC
  
  switch (par->metric_method_itype)
  {
    case PAR_METRIC_CALCULATE :
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

        break;

    case PAR_METRIC_IMPORT :
        if (myid==0) fprintf(stdout,"import descreted medium file ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 
        break;
  }

  if (par->is_export_metric==1)
  {
    if (myid==0) fprintf(stdout,"export metric to file ...\n"); 
    gd_curv_metric_export(blk->g3d,
                          blk->g3d_pos,
                          blk->g3d_name,
                          blk->g3d_num_of_vars,
                          blk->nx,
                          blk->ny,
                          blk->nz,
                          blk->ni1,
                          blk->nj1,
                          blk->nk1,
                          blk->ni,
                          blk->nj,
                          blk->nk,
                          blk->gni1,
                          blk->gnj1,
                          blk->gnk1,
                          blk->output_fname_part,
                          blk->grid_export_dir);
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
  
  // read or discrete velocity model
  switch (par->media_input_itype)
  {
    case PAR_MEDIA_CODE :
        if (myid==0) fprintf(stdout,"generate medium in code ...\n"); 
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
        break;

    case PAR_MEDIA_IMPORT :
        if (myid==0) fprintf(stdout,"import descreted medium file ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 
        //md_el_iso_import(blk->m3d, blk->nx, blk->ny, blk->nz, 
        //              myid3[0],myid3[1],myid3[2]);
        break;

    case PAR_MEDIA_3LAY :
        if (myid==0) fprintf(stdout,"read and descrete 3lay medium file ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 

        break;

    case PAR_MEDIA_3GRD :
        if (myid==0) fprintf(stdout,"read and descrete 3grd medium file ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 
        break;

  }

  if (par->is_export_media==1)
  {
    if (myid==0) fprintf(stdout,"export discrete medium to file ...\n"); 

    md_el_iso_export(blk->m3d,
                      blk->m3d_pos,
                      blk->m3d_name,
                      blk->m3d_num_of_vars,
                      blk->nx,
                      blk->ny,
                      blk->nz,
                      blk->ni1,
                      blk->nj1,
                      blk->nk1,
                      blk->ni,
                      blk->nj,
                      blk->nk,
                      blk->gni1,
                      blk->gnj1,
                      blk->gnk1,
                      blk->output_fname_part,
                      blk->media_export_dir);
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
  
  // read or gen source
  if (par->source_input_itype == PAR_SOURCE_FILE)
  {
    if (myid==0) fprintf(stdout,"read source file ...\n"); 
    if (myid==0) fprintf(stdout,"   not implemented yet\n"); 
  }
  else
  {
    if (myid==0) fprintf(stdout,"set single point source term in code ...\n"); 

    // default
    blk->num_of_force=0; blk->num_of_moment=0;
    if (par->source_input_itype == PAR_SOURCE_SINGLE_FORCE)  blk->num_of_force=1;
    if (par->source_input_itype == PAR_SOURCE_SINGLE_MOMENT) blk->num_of_moment=1;

    // set inner vars
    src_gen_single_point_gauss(blk->siz_line,
                               blk->siz_slice,
                               t0,
                               dt,
                               fd->num_rk_stages,
                               fd->rk_rhs_time,
                               blk->gni1,
                               blk->gni2,
                               blk->gnj1,
                               blk->gnj2,
                               blk->gnk1,
                               blk->gnk2,
                               blk->ni1,
                               blk->ni2,
                               blk->nj1,
                               blk->nj2,
                               blk->nk1,
                               blk->nk2,
                               fd->fd_half_len,
                               fd->fd_nghosts,
                               par->source_gridindex,
                               par->source_coords,
                               par->source_force_vector,
                               par->source_moment_tensor,
                               par->wavelet_name,
                               par->wavelet_coefs,
                               par->wavelet_tstart,
                               par->wavelet_tend,
                               &blk->num_of_force,
                               &blk->force_info,
                               &blk->force_vec_stf,
                               &blk->force_ext_indx,
                               &blk->force_ext_coef,
                               &blk->num_of_moment,
                               &blk->moment_info,
                               &blk->moment_ten_rate,
                               &blk->moment_ext_indx,
                               &blk->moment_ext_coef,
                               verbose);
  }

  /*
  if (par->is_export_source==1)
  {
      ierr = src_export();
  }
  */

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
//-- setup output, may require coord info
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"setup output info ...\n"); 

  // receiver: need to do
  io_read_locate_station(par->in_station_file,
                         blk->gni1, blk->gni2,
                         blk->gnj1, blk->gnj2,
                         blk->gnk1, blk->gnk2,
                         blk->ni1, blk->nj1, blk->nk1,
                         blk->siz_line, blk->siz_slice,
                         blk->c3d+blk->c3d_pos[0],
                         blk->c3d+blk->c3d_pos[1],
                         blk->c3d+blk->c3d_pos[2],
                         &(blk->num_of_sta),
                         &(blk->sta_name),
                         &(blk->sta_coord),
                         &(blk->sta_index),
                         &(blk->sta_shift));

  fd_blk_malloc_station(blk, nt_total);

  // inline
  fd_blk_locate_inline(blk,
                    fd->fdx_nghosts,
                    nt_total,
                    par->number_of_receiver_line,
                    par->receiver_line_index_start,
                    par->receiver_line_index_incre,
                    par->receiver_line_count,
                    par->receiver_line_name);
  
  // slice
  fd_blk_set_slice(blk,
                   fd->fdx_nghosts,
                   par->number_of_slice_x,
                   par->number_of_slice_y,
                   par->number_of_slice_z,
                   par->slice_x_index,
                   par->slice_y_index,
                   par->slice_z_index);
  
  // snapshot
  fd_blk_set_snapshot(blk,
                      fd->fdx_nghosts,
                      par->number_of_snapshot,
                      par->snapshot_name,
                      par->snapshot_index_start,
                      par->snapshot_index_count,
                      par->snapshot_index_incre,
                      par->snapshot_time_start,
                      par->snapshot_time_incre,
                      par->snapshot_save_velocity,
                      par->snapshot_save_stress,
                      par->snapshot_save_strain);

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
                   par->cartesian_grid_stepsize[0],
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
//-- qc
//-------------------------------------------------------------------------------
  
  fd_print(fd);
  fd_blk_print(blk);

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
                                  blk->force_ext_indx,blk->force_ext_coef,
                                  blk->num_of_moment, blk->moment_info, blk->moment_ten_rate,
                                  blk->moment_ext_indx,blk->moment_ext_coef,
                                  blk->num_of_sta, blk->sta_index,blk->sta_shift,blk->sta_seismo,
                                  blk->num_of_point, blk->point_loc_indx,blk->point_seismo,
                                  blk->num_of_slice_x, blk->slice_x_indx,blk->slice_x_fname,
                                  blk->num_of_slice_y, blk->slice_y_indx,blk->slice_y_fname,
                                  blk->num_of_slice_z, blk->slice_z_indx,blk->slice_z_fname,
                                  blk->num_of_snap, blk->snap_info, blk->snap_fname,
                                  blk->output_fname_part,
                                  blk->output_dir,
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
                                  par->output_all,
                                  verbose);
  
  time_t t_end = time(NULL);
  
  fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));

//-------------------------------------------------------------------------------
//-- save station and line seismo to sac
//-------------------------------------------------------------------------------

  for (int ir=0; ir<blk->num_of_sta; ir++)
  {
    char *sta_name = blk->sta_name[ir];
    int   num_of_vars = blk->w3d_num_of_vars;
    int   iptr_coord   = ir * FD_NDIM;

    // use fake evt_x etc. since did not implement gather evt_x by mpi
    float evt_x = 0.0;
    float evt_y = 0.0;
    float evt_z = 0.0;
    float evt_d = 0.0;

    //fprintf(stdout,"=== Debug: num_of_vars=%d\n",num_of_vars);fflush(stdout);
    for (int icmp=0; icmp<num_of_vars; icmp++)
    {
      //fprintf(stdout,"=== Debug: icmp=%d\n",icmp);fflush(stdout);

      int iptr_seismo = ir * num_of_vars * nt_total + icmp * nt_total;

      //fprintf(stdout,"=== Debug: icmp=%d,output_dir=%s\n",icmp,blk->output_dir);fflush(stdout);
      //fprintf(stdout,"=== Debug: icmp=%d,source_name=%s\n",icmp,par->source_name);fflush(stdout);
      //fprintf(stdout,"=== Debug: icmp=%d,sta_name=%s\n",icmp,sta_name);fflush(stdout);
      //fprintf(stdout,"=== Debug: icmp=%d,w3d_name=%s\n",icmp,blk->w3d_name[icmp]);fflush(stdout);

      sprintf(ou_file,"%s/%s.%s.%s.sac", blk->output_dir,par->source_name,
                  sta_name,blk->w3d_name[icmp]);

      //fprintf(stdout,"=== Debug: icmp=%d,ou_file=%s\n",icmp,ou_file);fflush(stdout);

      sacExport1C1R(ou_file,
            blk->sta_seismo+iptr_seismo,
            evt_x, evt_y, evt_z, evt_d,
            blk->sta_coord[iptr_coord+0],
            blk->sta_coord[iptr_coord+1],
            blk->sta_coord[iptr_coord+2],
            0.0, dt, nt_total, err_message
            );
    }
  }

  for (int ir=0; ir<blk->num_of_point; ir++)
  {
    int   num_of_vars = blk->w3d_num_of_vars;
    int   iptr_coord  = ir * FD_NDIM;
    int   line_sno  = blk->point_line_sno[ir];
    int   line_offset  = blk->point_line_offset[ir];
    char *line_name = par->receiver_line_name[line_sno];

    //fprintf(stdout,"=== Debug: ir=%d,line_sno=%d,name=%s,nvar=%d\n",
    //    ir,line_sno,line_name,num_of_vars);fflush(stdout);

    // use fake evt_x etc. since did not implement gather evt_x by mpi
    float evt_x = 0.0;
    float evt_y = 0.0;
    float evt_z = 0.0;
    float evt_d = 0.0;

    //fprintf(stdout,"=== Debug: num_of_vars=%d\n",num_of_vars);fflush(stdout);

    for (int icmp=0; icmp<num_of_vars; icmp++)
    {
      int iptr_seismo = ir * num_of_vars * nt_total + icmp * nt_total;

      // evt1.line2.pt2.Vx.sac
      sprintf(ou_file,"%s/%s.%s.no%d.%s.sac", blk->output_dir,par->source_name,
                  line_name,line_offset,blk->w3d_name[icmp]);

      //fprintf(stdout,"=== Debug: icmp=%d,ou_file=%s\n",icmp,ou_file);fflush(stdout);

      sacExport1C1R(ou_file,
            blk->point_seismo+iptr_seismo,
            evt_x, evt_y, evt_z, evt_d,
            blk->point_coord[iptr_coord+0],
            blk->point_coord[iptr_coord+1],
            blk->point_coord[iptr_coord+2],
            0.0, dt, nt_total, err_message
            );
    }
  }

//-------------------------------------------------------------------------------
//-- postprocess
//-------------------------------------------------------------------------------

  MPI_Finalize();

  return 0;
}
