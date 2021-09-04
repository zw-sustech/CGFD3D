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

#include "constants.h"
#include "par_t.h"
// blk_t.h contains most other headers
#include "blk_t.h"

#include "media_discrete_model.h"
#include "sv_eq1st_cart_col.h"

int main(int argc, char** argv)
{
  int verbose = 1; // default fprint
  char *par_fname;
  char err_message[CONST_MAX_STRLEN];

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

  par_t *par = (par_t *) malloc(sizeof(par_t));

  par_mpi_get(par_fname, myid, comm, par, verbose);

  if (myid==0 && verbose>0) par_print(par);

//-------------------------------------------------------------------------------
// init blk_t
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"create blk ...\n"); 

  // malloc blk
  blk_t *blk = (blk_t *) malloc(sizeof(blk_t));

  // malloc inner vars
  blk_init(blk, myid, verbose);

  fd_t            *fd            = blk->fd    ;
  mympi_t         *mympi         = blk->mympi ;
  gdinfo_t        *gdinfo        = blk->gdinfo;
  gd_t            *gdcart            = blk->gd;
  gdcurv_metric_t *gdcurv_metric = blk->gdcurv_metric;
  md_t            *md            = blk->md;
  wav_t           *wav           = blk->wav;
  src_t           *src           = blk->src;
  bdryfree_t      *bdryfree      = blk->bdryfree;
  bdrypml_t       *bdrypml       = blk->bdrypml;
  iorecv_t        *iorecv        = blk->iorecv;
  ioline_t        *ioline        = blk->ioline;
  ioslice_t       *ioslice       = blk->ioslice;
  iosnap_t        *iosnap        = blk->iosnap;

  // set up fd_t
  //    not support selection scheme by par file yet
  if (myid==0 && verbose>0) fprintf(stdout,"set scheme ...\n"); 
  fd_set_macdrp(fd);

  // set mpi
  if (myid==0 && verbose>0) fprintf(stdout,"set mpi topo ...\n"); 
  mympi_set(mympi,
            par->number_of_mpiprocs_x,
            par->number_of_mpiprocs_y,
            comm,
            myid, verbose);

  // set gdinfo
  gd_info_set(gdinfo, mympi,
              par->number_of_total_grid_points_x,
              par->number_of_total_grid_points_y,
              par->number_of_total_grid_points_z,
              par->abs_num_of_layers,
              fd->fdx_nghosts,
              fd->fdy_nghosts,
              fd->fdz_nghosts,
              verbose);

  // set str in blk
  blk_set_output(blk, mympi,
                 par->output_dir,
                 par->grid_export_dir,
                 par->media_export_dir,
                 verbose);

//-------------------------------------------------------------------------------
//-- grid generation or import
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"allocate grid vars ...\n"); 

  if (par->grid_generation_itype != PAR_GRID_CARTESIAN) {
    fprintf(stderr, "ERROR: grid type wrong\n");
  }

  // malloc var in gdcart
  float dx = par->cartesian_grid_stepsize[0];
  float dy = par->cartesian_grid_stepsize[1];
  float dz = par->cartesian_grid_stepsize[2];

  float x0 = par->cartesian_grid_origin[0];
  float y0 = par->cartesian_grid_origin[1];
  float z0 = par->cartesian_grid_origin[2];

  gd_cart_init_set(gdinfo,gdcart,dx,x0,dy,y0,dz,z0);

  fprintf(stdout, " --> done\n"); fflush(stdout);

  // output
  if (par->is_export_grid==1)
  {
    if (myid==0) fprintf(stdout,"export coord to file ...\n"); 
    gd_cart_coord_export(gdinfo, gdcart,
                         blk->output_fname_part,
                         blk->grid_export_dir);
  }
  fprintf(stdout, " --> done\n"); fflush(stdout);

//-------------------------------------------------------------------------------
//-- media generation or import
//-------------------------------------------------------------------------------

  // allocate media vars
  if (myid==0 && verbose>0) fprintf(stdout,"allocate media vars ...\n"); 
  md_init(gdinfo, md, par->media_itype, par->visco_itype);

  // read or discrete velocity model
  switch (par->media_input_itype)
  {
    case PAR_MEDIA_CODE :
        if (myid==0) fprintf(stdout,"generate simple medium in code ...\n"); 

        if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          md_gen_test_el_iso(md);
        }

        if (md->medium_type == CONST_MEDIUM_ELASTIC_ANISO) {
          md_gen_test_el_aniso(md);
        }

        if (md->visco_type == CONST_VISCO_GRAVES_QS) {
          md_gen_test_Qs(md, par->visco_Qs_freq);
        }

        break;

    case PAR_MEDIA_IMPORT :
        if (myid==0) fprintf(stdout,"import discrete medium file ...\n"); 
        if (myid==0) fprintf(stdout,"   not implemented yet\n"); 
        //md_import(blk->m3d, blk->nx, blk->ny, blk->nz, 
        //              myid3[0],myid3[1],myid3[2]);
        break;

    case PAR_MEDIA_3LAY : {
        if (myid==0) fprintf(stdout,"read and discretize 3D layer medium file ...\n"); 

        float *lam3d = md->lambda;
        float  *mu3d = md->mu;
        float *rho3d = md->rho;
        //media_el_iso_layer2model(lam3d, mu3d, rho3d,
        //                         gdcart->x3d,
        //                         gdcart->y3d,
        //                         gdcart->z3d,
        //                         gdcart->nx,
        //                         gdcart->ny,
        //                         gdcart->nz,
        //                         par->media_input_file,
        //                         par->equivalent_medium_method);
        break;
    }

    case PAR_MEDIA_3GRD : {
        if (myid==0) fprintf(stdout,"read and descretize 3D grid medium file ...\n"); 

        float *lam3d = md->lambda;
        float  *mu3d = md->mu;
        float *rho3d = md->rho;

        //media_el_iso_grid2model(lam3d, mu3d, rho3d,
        //                        gdcart->x3d,
        //                        gdcart->y3d,
        //                        gdcart->z3d,
        //                        gdcart->nx,
        //                        gdcart->ny,
        //                        gdcart->nz,
        //                        gdcart->xmin, gdcart->xmax,   //float Xmin, float Xmax,
        //                        gdcart->ymin, gdcart->ymax,   //float Ymin, float Ymax, 
        //                        par->media_input_file,
        //                        par->equivalent_medium_method); 
        break;
    }
  }

  // export grid media
  if (par->is_export_media==1)
  {
    if (myid==0) fprintf(stdout,"export discrete medium to file ...\n"); 

    md_export(gdinfo, md,
              blk->output_fname_part,
              blk->media_export_dir);
  } else {
    if (myid==0) fprintf(stdout,"do not export medium\n"); 
  }

//-------------------------------------------------------------------------------
//-- estimate/check/set time step
//-------------------------------------------------------------------------------

  float   t0 = par->time_start;
  float   dt = par->size_of_time_step;
  int     nt_total = par->number_of_time_steps+1;

  if (par->time_check_stability==1)
  {
    float dt_est[mpi_size];
    float dtmax, dtmaxVp, dtmaxL;
    int   dtmaxi, dtmaxj, dtmaxk;

    //-- estimate time step
    if (myid==0) fprintf(stdout,"   estimate time step ...\n"); 
    blk_dt_esti_cart(gdinfo, gdcart,md,fd->CFL,
            &dtmax, &dtmaxVp, &dtmaxL, &dtmaxi, &dtmaxj, &dtmaxk);
    
    //-- print for QC
    fprintf(stdout, "-> topoid=[%d,%d], dtmax=%f, Vp=%f, L=%f, i=%d, j=%d, k=%d\n",
            mympi->topoid[0],mympi->topoid[1], dtmax, dtmaxVp, dtmaxL, dtmaxi, dtmaxj, dtmaxk);
    
    // receive dtmax from each proc
    MPI_Allgather(&dtmax,1,MPI_REAL,dt_est,1,MPI_REAL,MPI_COMM_WORLD);
  
    if (myid==0)
    {
       int dtmax_mpi_id = 0;
       dtmax = 1e19;
       for (int n=0; n < mpi_size; n++)
       {
        fprintf(stdout,"max allowed dt at each proc: id=%d, dtmax=%f\n", n, dt_est[n]);
        if (dtmax > dt_est[n]) {
          dtmax = dt_est[n];
          dtmax_mpi_id = n;
        }
       }
       fprintf(stdout,"Global maximum allowed time step is %f at thread %d\n", dtmax, dtmax_mpi_id);

       // check valid
       if (dtmax <= 0.0) {
          fprintf(stderr,"ERROR: maximum dt <= 0, stop running\n");
          MPI_Abort(MPI_COMM_WORLD,-1);
       }

       //-- auto set stept
       if (dt < 0.0) {
          dt       = blk_keep_two_digi(dtmax);
          nt_total = (int) (par->time_window_length / dt + 0.5);

          fprintf(stdout, "-> Set dt       = %f according to maximum allowed value\n", dt);
          fprintf(stdout, "-> Set nt_total = %d\n", nt_total);
       }

       //-- if input dt, check value
       if (dtmax < dt) {
          fprintf(stdout, "Serious Error: dt=%f > dtmax=%f, stop!\n", dt, dtmax);
          MPI_Abort(MPI_COMM_WORLD, -1);
       }
    }
    
    //-- from root to all threads
    MPI_Bcast(&dt      , 1, MPI_REAL, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nt_total, 1, MPI_INT , 0, MPI_COMM_WORLD);
  }

//-------------------------------------------------------------------------------
//-- source import or locate on fly
//-------------------------------------------------------------------------------
  
  if (par->source_input_itype == PAR_SOURCE_FILE)
  {
    char fnm_suffix[6];
    int  n = strlen(par->source_input_file);
    strncpy(fnm_suffix,par->source_input_file+n-6,6);
    if(strcmp(fnm_suffix,"anasrc")==0)
    {
      if (myid==0) fprintf(stdout,"***input source type is analysis***\n");
    }
    if(strcmp(fnm_suffix,"valsrc")==0)
    {
      if (myid==0) fprintf(stdout,"***input source type is value sample***\n");
    }
  }
  else
  {
    if (myid==0) fprintf(stdout,"set source using info from par file ...\n"); 

    src_set_by_par(gdinfo, gdcart, src,
                   t0, dt,
                   fd->num_rk_stages, fd->rk_rhs_time,
                   fd->fdx_max_half_len,
                   par->source_name,
                   par->source_number,
                   par->source_index,
                   par->source_inc,
                   par->source_coords,
                   par->source_force_vector,
                   par->source_force_actived,
                   par->source_moment_tensor,
                   par->source_moment_actived,
                   par->wavelet_name,
                   par->wavelet_coefs,
                   par->wavelet_tstart,
                   par->wavelet_tend,
                   comm, myid, verbose);
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
  wav_init(gdinfo, wav, fd->num_rk_stages);

//-------------------------------------------------------------------------------
//-- setup output, may require coord info
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"setup output info ...\n"); 

  // receiver: need to do
  io_recv_read_locate(gdinfo, gdcart, iorecv,
                      nt_total, wav->ncmp, par->in_station_file);

  // line
  io_line_locate(gdinfo, gdcart, ioline,
                 wav->ncmp,
                 nt_total,
                 par->number_of_receiver_line,
                 par->receiver_line_index_start,
                 par->receiver_line_index_incre,
                 par->receiver_line_count,
                 par->receiver_line_name);
  
  // slice
  io_slice_locate(gdinfo, ioslice,
                  par->number_of_slice_x,
                  par->number_of_slice_y,
                  par->number_of_slice_z,
                  par->slice_x_index,
                  par->slice_y_index,
                  par->slice_z_index,
                  blk->output_fname_part,
                  blk->output_dir);
  
  // snapshot
  io_snapshot_locate(gdinfo, iosnap,
                     par->number_of_snapshot,
                     par->snapshot_name,
                     par->snapshot_index_start,
                     par->snapshot_index_count,
                     par->snapshot_index_incre,
                     par->snapshot_time_start,
                     par->snapshot_time_incre,
                     par->snapshot_save_velocity,
                     par->snapshot_save_stress,
                     par->snapshot_save_strain,
                     blk->output_fname_part,
                     blk->output_dir);

//-------------------------------------------------------------------------------
//-- absorbing boundary etc auxiliary variables
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"setup absorbingg boundary ...\n"); 
  
  if (par->bdry_has_cfspml == 1)
  {
    bdry_pml_set(gdinfo, gdcart, wav, bdrypml,
                 mympi->neighid,
                 par->cfspml_is_sides,
                 par->abs_num_of_layers,
                 par->cfspml_alpha_max,
                 par->cfspml_beta_max,
                 par->cfspml_velocity,
                 verbose);
  }

//-------------------------------------------------------------------------------
//-- free surface preproc
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"cal free surface matrix ...\n"); 

  if (par->bdry_has_free == 1)
  {
    bdry_free_set(gdinfo,bdryfree, mympi->neighid, par->free_is_sides, verbose);
  }

//-------------------------------------------------------------------------------
//-- setup mesg
//-------------------------------------------------------------------------------

  if (myid==0 && verbose>0) fprintf(stdout,"init mesg ...\n"); 
  blk_colcent_mesg_init(mympi, gdinfo->ni, gdinfo->nj, gdinfo->nk,
                  fd->fdx_nghosts, fd->fdy_nghosts, wav->ncmp);

//-------------------------------------------------------------------------------
//-- qc
//-------------------------------------------------------------------------------
  
  fd_print(fd);

  blk_print(blk);

  gd_info_print(gdinfo);

  ioslice_print(ioslice);

  iosnap_print(iosnap);

//-------------------------------------------------------------------------------
//-- slover
//-------------------------------------------------------------------------------
  
  // convert rho to 1 / rho to reduce number of arithmetic cal
  md_rho_to_slow(md->rho, md->siz_icmp);

  if (myid==0 && verbose>0) fprintf(stdout,"start solver ...\n"); 
  
  time_t t_start = time(NULL);
  
  sv_eq1st_cart_col_allstep(fd,gdinfo,gdcart,md,
                            src,bdryfree,bdrypml,
                            wav, mympi,
                            iorecv,ioline,ioslice,iosnap,
                            dt,nt_total,t0,
                            blk->output_fname_part,
                            blk->output_dir,
                            par->check_nan_every_nummber_of_steps,
                            par->output_all,
                            verbose);
  
  time_t t_end = time(NULL);
  
  fprintf(stdout,"\n\nRuning Time of time :%f s \n", difftime(t_end,t_start));

//-------------------------------------------------------------------------------
//-- save station and line seismo to sac
//-------------------------------------------------------------------------------

  io_recv_output_sac(iorecv,dt,wav->ncmp,wav->cmp_name,
                      src->evtnm,blk->output_dir,err_message);

  io_line_output_sac(ioline,dt,wav->cmp_name,src->evtnm,blk->output_dir);

//-------------------------------------------------------------------------------
//-- postprocess
//-------------------------------------------------------------------------------

  MPI_Finalize();

  return 0;
}
