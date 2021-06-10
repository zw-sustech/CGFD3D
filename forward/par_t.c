/*
 * 
 */

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "par_t.h"

/*
 * for MPI, master read, broadcast to all procs
 */
void
par_mpi_get(char *par_fname, int myid, MPI_Comm comm, struct par_t *par, int verbose)
{
  char *str;

  if (myid==0)
  {
    FILE *fp=fopen(par_fname,"r");
    if (!fp) {
      fprintf(stderr,"Error: can't open par file: %s\n", par_fname);
      MPI_Finalize();
      exit(1);
    }
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);

    // bcast len to all nodes
    MPI_Bcast(&len, 1, MPI_LONG, 0, comm);

    fseek(fp, 0, SEEK_SET);
    str = (char*)malloc(len+1);
    fread(str, 1, len, fp);
    fclose(fp);

    // bcast str
    MPI_Bcast(str, len+1, MPI_CHAR, 0, comm);
  }
  else
  {
    long len;
    // get len 
    MPI_Bcast(&len, 1, MPI_LONG, 0, comm);

    str = (char*)malloc(len+1);
    // get str
    MPI_Bcast(str, len+1, MPI_CHAR, 0, comm);
  }
    
  par_read_from_str(str, par);

  free(str);
}

/*
 * for non-MPI, read from file
 */
void
par_read_from_file(char *par_fname, int myid, MPI_Comm comm, struct par_t *par, int verbose)
{
  //
  // read whole file inot str
  //
  FILE *fp = fopen(par_fname,"r");

  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);

  fseek(fp, 0, SEEK_SET);
  char *str = (char*)malloc(len+1);
  fread(str, 1, len, fp);
  fclose(fp);

  // read from str
  par_read_from_str(str, par);
}

/*
 * funcs to get par from alread read in str
 */
void 
par_read_from_str(const char *str, struct par_t *par)
{
  // allocate
  par->boundary_type_name = (char **)malloc(FD_NDIM_2 * sizeof(char*));
  for (int i=0; i<FD_NDIM_2; i++) {
    par->boundary_type_name[i] = (char *)malloc(10*sizeof(char));
  }

  // convert str to json
  cJSON *root = cJSON_Parse(str);
  if (NULL == root) {
    printf("Error at parsing json!\n");
    exit(-1);
  }

  cJSON *item;
  cJSON *subitem, *thirditem, *snapitem, *lineitem;

  // no default
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_x")) {
    par->number_of_total_grid_points_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_y")) {
    par->number_of_total_grid_points_y = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_z")) {
    par->number_of_total_grid_points_z = item->valueint;
  }

  // default mpi threads
  par->number_of_mpiprocs_x = 1;
  par->number_of_mpiprocs_y = 1;

  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_x")) {
    par->number_of_mpiprocs_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_y")) {
    par->number_of_mpiprocs_y = item->valueint;
  }

  // set default values to negative
  par->size_of_time_step = -1.0;
  par->number_of_time_steps = -1;
  par->length_of_time_window_in_second = -1.0;

  if (item = cJSON_GetObjectItem(root, "size_of_time_step")) {
    par->size_of_time_step = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_time_steps")) {
    par->number_of_time_steps = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "length_of_time_window_in_second")) {
    par->length_of_time_window_in_second = item->valuedouble;
  }
  if (par->size_of_time_step > 0.0 && par->number_of_time_steps > 0) {
    par->length_of_time_window_in_second = par->size_of_time_step *
                                           par->number_of_time_steps;
  } else {
    fprintf(stderr,"Error: size_of_time_step=%f, number_of_time_steps=%d\n",
           par->size_of_time_step, par->number_of_time_steps);
    fprintf(stderr,"       auto estimate time step has not been implemented\n");
  }

  par->time_start = 0.0;
  par->time_end   = par->time_start + 
          par->number_of_time_steps * par->size_of_time_step;
  //int     nt_total = (int) ((par->time_end - par->time_start) / dt+0.5);

  // cfspml default values
  for (int i = 0; i < FD_NDIM_2; i++) {
    par->abs_num_of_layers[i] = 0;
  }
  // x1 boundary, no default yet
  if (item = cJSON_GetObjectItem(root, "boundary_x_left"))
  {
    if (subitem = cJSON_GetObjectItem(item, "cfspml")) {
       sprintf(par->boundary_type_name[0], "%s", "cfspml");
       par_read_json_cfspml(subitem,
                            par->abs_num_of_layers+0,
                            par->cfspml_alpha_max +0,
                            par->cfspml_beta_max  +0,
                            par->cfspml_velocity  +0);
    }
  }
  // x2 boundary, no default yet
  if (item = cJSON_GetObjectItem(root, "boundary_x_right"))
  {
    if (subitem = cJSON_GetObjectItem(item, "cfspml")) {
       sprintf(par->boundary_type_name[1], "%s", "cfspml");
       par_read_json_cfspml(subitem,
                            par->abs_num_of_layers+1,
                            par->cfspml_alpha_max +1,
                            par->cfspml_beta_max  +1,
                            par->cfspml_velocity  +1);
    }
  }
  // y1 boundary, no default yet
  if (item = cJSON_GetObjectItem(root, "boundary_y_front"))
  {
    if (subitem = cJSON_GetObjectItem(item, "cfspml")) {
       sprintf(par->boundary_type_name[2], "%s", "cfspml");
       par_read_json_cfspml(subitem,
                            par->abs_num_of_layers+2,
                            par->cfspml_alpha_max +2,
                            par->cfspml_beta_max  +2,
                            par->cfspml_velocity  +2);
    }
  }
  // y2 boundary, no default yet
  if (item = cJSON_GetObjectItem(root, "boundary_y_back"))
  {
    if (subitem = cJSON_GetObjectItem(item, "cfspml")) {
       sprintf(par->boundary_type_name[3], "%s", "cfspml");
       par_read_json_cfspml(subitem,
                            par->abs_num_of_layers+3,
                            par->cfspml_alpha_max +3,
                            par->cfspml_beta_max  +3,
                            par->cfspml_velocity  +3);
    }
  }
  // z1 boundary, no default yet
  if (item = cJSON_GetObjectItem(root, "boundary_z_bottom"))
  {
    if (subitem = cJSON_GetObjectItem(item, "cfspml")) {
       sprintf(par->boundary_type_name[4], "%s", "cfspml");
       par_read_json_cfspml(subitem,
                            par->abs_num_of_layers+4,
                            par->cfspml_alpha_max +4,
                            par->cfspml_beta_max  +4,
                            par->cfspml_velocity  +4);
    }
  }
  // z2 boundary, no default yet
  if (item = cJSON_GetObjectItem(root, "boundary_z_top"))
  {
    if (subitem = cJSON_GetObjectItem(item, "cfspml")) {
       sprintf(par->boundary_type_name[5], "%s", "cfspml");
       par_read_json_cfspml(subitem,
                            par->abs_num_of_layers+5,
                            par->cfspml_alpha_max +5,
                            par->cfspml_beta_max  +5,
                            par->cfspml_velocity  +5);
    }
    if (subitem = cJSON_GetObjectItem(item, "free")) {
       sprintf(par->boundary_type_name[5], "%s", "free");
    }
  }

  //
  //-- grid
  //

  // default output grid

  par->grid_generation_itype = PAR_GRID_IMPORT;
  if (item = cJSON_GetObjectItem(root, "grid_generation_method")) {
    // import grid
    if (subitem = cJSON_GetObjectItem(item, "import")) {
       par->grid_generation_itype = PAR_GRID_IMPORT;
       sprintf(par->grid_import_dir, "%s", subitem->valuestring);
    }
    // generate cartesian grid
    if (subitem = cJSON_GetObjectItem(item, "cartesian")) {
       par->grid_generation_itype = PAR_GRID_CARTESIAN;
       if (thirditem = cJSON_GetObjectItem(subitem, "origin")) {
         for (int i = 0; i < FD_NDIM; i++) {
           par->cartesian_grid_origin[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "inteval")) {
         for (int i = 0; i < FD_NDIM; i++) {
           par->cartesian_grid_stepsize[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
       }
    }
    // layer interp
    if (subitem = cJSON_GetObjectItem(item, "layer_interp")) {
       par->grid_generation_itype = PAR_GRID_LAYER_INTERP;
       if (thirditem = cJSON_GetObjectItem(subitem, "in_grid_layer_file")) {
          sprintf(par->in_grid_layer_file, "%s", thirditem->valuestring);
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "refine_factor")) {
         for (int i = 0; i < FD_NDIM; i++) {
           par->grid_layer_interp_factor[i] = cJSON_GetArrayItem(thirditem, i)->valueint;
         }
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "horizontal_start_point")) {
         for (int i = 0; i < FD_NDIM-1; i++) {
           par->grid_layer_startend[i] = cJSON_GetArrayItem(thirditem, i)->valueint;
         }
       }
       if (thirditem = cJSON_GetObjectItem(subitem, "vertical_end_point")) {
         par->grid_layer_startend[FD_NDIM-1] = thirditem->valueint;
       }
    }
  }

  par->is_export_grid = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_grid")) {
     par->is_export_grid = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "grid_export_dir")) {
      sprintf(par->grid_export_dir,"%s",item->valuestring);
  }

  //
  //-- metric
  //

  par->metric_method_itype = PAR_METRIC_CALCULATE;
  if (item = cJSON_GetObjectItem(root, "metric_calculation_method")) {
    if (subitem = cJSON_GetObjectItem(item, "import")) {
        par->metric_method_itype = PAR_METRIC_IMPORT;
        sprintf(par->metric_import_dir, "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "calculate")) {
        par->metric_method_itype = PAR_METRIC_CALCULATE;
    }
  }

  par->is_export_metric = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_metric")) {
     par->is_export_metric = item->valueint;
  }

  //
  //-- medium
  //

  par->media_input_itype = PAR_MEDIA_IMPORT;
  if (item = cJSON_GetObjectItem(root, "media_input")) {
    if (subitem = cJSON_GetObjectItem(item, "import")) {
        par->media_input_itype = PAR_MEDIA_IMPORT;
        sprintf(par->media_import_dir, "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "code_generate")) {
        par->media_input_itype = PAR_MEDIA_CODE;
    }
    if (subitem = cJSON_GetObjectItem(item, "in_3lay_file")) {
        par->media_input_itype = PAR_MEDIA_3LAY;
        sprintf(par->media_input_file, "%s", subitem->valuestring);
        // If input layer model, choose which equivalent medium para method
        if (thirditem = cJSON_GetObjectItem(item, "equivalent_medium_method")) {
          sprintf(par->equivalent_medium_method, "%s", thirditem->valuestring);
        }
    }
    if (subitem = cJSON_GetObjectItem(item, "in_3grd_file")) {
        par->media_input_itype = PAR_MEDIA_3GRD;
        sprintf(par->media_input_file, "%s", subitem->valuestring);
    }
  }

  par->is_export_media = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_media")) {
     par->is_export_media = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "media_export_dir")) {
      sprintf(par->media_export_dir,"%s",item->valuestring);
  }

  //
  //-- source
  //

  // default: grid index to negative for determine grid or loc

  for (int i = 0; i < FD_NDIM; i++) par->source_gridindex[i] = -1;
  
  par->source_input_itype = PAR_SOURCE_SINGLE_FORCE;
  if (item = cJSON_GetObjectItem(root, "source_input")) {
    if (subitem = cJSON_GetObjectItem(item, "single_force"))
    {
       par->source_input_itype = PAR_SOURCE_SINGLE_FORCE;
       if (thirditem = cJSON_GetObjectItem(subitem, "force_vector")) {
         for (int i = 0; i < FD_NDIM; i++) {
           par->source_force_vector[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
       }
       par_read_json_source(subitem,"source_time_functon",
                par->source_name,par->source_coords, par->source_gridindex,
                par->wavelet_name, par->wavelet_coefs,
                &(par->wavelet_tstart), &(par->wavelet_tend));
    }

    if (subitem = cJSON_GetObjectItem(item, "single_moment"))
    {
       par->source_input_itype = PAR_SOURCE_SINGLE_MOMENT;
       if (thirditem = cJSON_GetObjectItem(subitem, "moment_tensor")) {
         for (int i = 0; i < FD_NDIM_2; i++) {
           par->source_moment_tensor[i] = cJSON_GetArrayItem(thirditem, i)->valuedouble;
         }
       }
       par_read_json_source(subitem,"moment_rate_functon",
                par->source_name,par->source_coords, par->source_gridindex,
                par->wavelet_name, par->wavelet_coefs,
                &(par->wavelet_tstart), &(par->wavelet_tend));
    }

    if (subitem = cJSON_GetObjectItem(item, "in_source_file"))
    {
        par->source_input_itype = PAR_SOURCE_FILE;
        sprintf(par->source_input_file, "%s", subitem->valuestring);
    }
  }

  par->is_export_source = 1;
  if (item = cJSON_GetObjectItem(root, "is_export_source")) {
     par->is_export_source = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "source_export_dir")) {
      sprintf(par->source_export_dir,"%s",item->valuestring);
  }

  //-- output dir
  if (item = cJSON_GetObjectItem(root, "output_dir")) {
      sprintf(par->output_dir,"%s",item->valuestring);
  }

  //-- receiver
  if (item = cJSON_GetObjectItem(root, "in_station_file")) {
    sprintf(par->in_station_file, "%s", item->valuestring);
  }

  //-- receiver line
  if (item = cJSON_GetObjectItem(root, "receiver_line"))
  {
    par->number_of_receiver_line = cJSON_GetArraySize(item);
    par->receiver_line_index_start  = (int *)malloc(par->number_of_receiver_line*sizeof(int)*FD_NDIM);
    par->receiver_line_index_incre  = (int *)malloc(par->number_of_receiver_line*sizeof(int)*FD_NDIM);
    par->receiver_line_count  = (int *)malloc(par->number_of_receiver_line*sizeof(int));
    //par->receiver_line_time_interval  = (int *)malloc(par->number_of_receiver_line*sizeof(int));
    par->receiver_line_name = (char **)malloc(par->number_of_receiver_line*sizeof(char*));
    for (int n=0; n<par->number_of_receiver_line; n++) {
      par->receiver_line_name[n] = (char *)malloc(PAR_MAX_STRLEN*sizeof(char));
    }
    // each line
    for (int i=0; i < cJSON_GetArraySize(item) ; i++)
    {
      lineitem = cJSON_GetArrayItem(item, i);

      if (subitem = cJSON_GetObjectItem(lineitem, "name"))
      {
        sprintf(par->receiver_line_name[i],"%s",subitem->valuestring);
      }

      if (subitem = cJSON_GetObjectItem(lineitem, "grid_index_start"))
      {
        for (int j = 0; j < FD_NDIM; j++) {
          par->receiver_line_index_start[i*FD_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(lineitem, "grid_index_incre"))
      {
        for (int j = 0; j < FD_NDIM; j++) {
          par->receiver_line_index_incre[i*FD_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(lineitem, "grid_index_count"))
      {
         par->receiver_line_count[i] = subitem->valueint;
      }

      //if (subitem = cJSON_GetObjectItem(lineitem, "t_index_interval"))
      //{
      //   par->receiver_line_tinterval[i] = cJSON_GetArrayItem(subitem, j)->valueint;
      //}
    }
  }

  // slice
  if (item = cJSON_GetObjectItem(root, "slice"))
  {
    if (subitem = cJSON_GetObjectItem(item, "x_index"))
    {
      par->number_of_slice_x = cJSON_GetArraySize(subitem);
      par->slice_x_index  = (int *)malloc(par->number_of_slice_x*sizeof(int));
      for (int i=0; i < par->number_of_slice_x ; i++)
      {
        par->slice_x_index[i] = cJSON_GetArrayItem(subitem, i)->valueint;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "y_index"))
    {
      par->number_of_slice_y = cJSON_GetArraySize(subitem);
      par->slice_y_index  = (int *)malloc(par->number_of_slice_y*sizeof(int));
      for (int i=0; i < par->number_of_slice_y ; i++)
      {
        par->slice_y_index[i] = cJSON_GetArrayItem(subitem, i)->valueint;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "z_index"))
    {
      par->number_of_slice_z = cJSON_GetArraySize(subitem);
      par->slice_z_index  = (int *)malloc(par->number_of_slice_z*sizeof(int));
      for (int i=0; i < par->number_of_slice_z ; i++)
      {
        par->slice_z_index[i] = cJSON_GetArrayItem(subitem, i)->valueint;
      }
    }
  }

  // snapshot
  if (item = cJSON_GetObjectItem(root, "snapshot"))
  {
    par->number_of_snapshot = cJSON_GetArraySize(item);
    //fprintf(stdout,"size=%d, %d, %d\n", par->number_of_snapshot, sizeof(int), FD_NDIM);
    //fflush(stdout);
    par->snapshot_index_start  = (int *)malloc(par->number_of_snapshot*sizeof(int)*FD_NDIM);
    //if (par->snapshot_index_start == NULL) {
    //  fprintf(stdout,"eror\n");
    //  fflush(stdout);
    //}
    par->snapshot_index_count  = (int *)malloc(par->number_of_snapshot*sizeof(int)*FD_NDIM);
    par->snapshot_index_incre = (int *)malloc(par->number_of_snapshot*sizeof(int)*FD_NDIM);
    par->snapshot_time_start  = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_time_incre = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_save_velocity = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_save_stress  = (int *)malloc(par->number_of_snapshot*sizeof(int));
    par->snapshot_save_strain = (int *)malloc(par->number_of_snapshot*sizeof(int));
    // name of snapshot
    par->snapshot_name = (char **)malloc(par->number_of_snapshot*sizeof(char*));
    for (int n=0; n<par->number_of_snapshot; n++) {
      par->snapshot_name[n] = (char *)malloc(PAR_MAX_STRLEN*sizeof(char));
    }

    // each snapshot
    for (int i=0; i < cJSON_GetArraySize(item) ; i++)
    {
      snapitem = cJSON_GetArrayItem(item, i);

      if (subitem = cJSON_GetObjectItem(snapitem, "name"))
      {
        sprintf(par->snapshot_name[i],"%s",subitem->valuestring);
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "grid_index_start"))
      {
        for (int j = 0; j < FD_NDIM; j++) {
          par->snapshot_index_start[i*FD_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "grid_index_count"))
      {
        for (int j = 0; j < FD_NDIM; j++) {
          par->snapshot_index_count[i*FD_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "grid_index_incre"))
      {
        for (int j = 0; j < FD_NDIM; j++) {
          par->snapshot_index_incre[i*FD_NDIM+j] = cJSON_GetArrayItem(subitem, j)->valueint;
        }
      }

      if (subitem = cJSON_GetObjectItem(snapitem, "time_index_start")) {
        par->snapshot_time_start[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "time_index_incre")) {
        par->snapshot_time_incre[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "save_velocity")) {
        par->snapshot_save_velocity[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "save_stress")) {
        par->snapshot_save_stress[i]  = subitem->valueint;
      }
      if (subitem = cJSON_GetObjectItem(snapitem, "save_strain")) {
        par->snapshot_save_strain[i] = subitem->valueint;
      }
    }
  }

  //-- misc
  if (item = cJSON_GetObjectItem(root, "check_nan_every_nummber_of_steps")) {
      par->check_nan_every_nummber_of_steps = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "output_all")) {
      par->output_all = item->valueint;
  }

  //if (item = cJSON_GetObjectItem(root, "grid_name")) {
  //    sprintf(par->grid_name,"%s",item->valuestring);
  //}

  cJSON_Delete(root);

  // set values to default ones if no input

  return;
}

/*
 * funcs to read cfspml para from json str
 */
void 
par_read_json_cfspml(cJSON *item,
      int *nlay, float *amax, float *bmax, float *vel)
{
  cJSON *subitem;

  if (subitem = cJSON_GetObjectItem(item, "number_of_layers"))
  {
    *nlay = subitem->valueint;
  }
  if (subitem = cJSON_GetObjectItem(item, "alpha_max"))
  {
    *amax = subitem->valuedouble;
  }
  if (subitem = cJSON_GetObjectItem(item, "beta_max"))
  {
    *bmax = subitem->valuedouble;
  }
  if (subitem = cJSON_GetObjectItem(item, "ref_vel"))
  {
    *vel = subitem->valuedouble;
  }
}

/*
 * funcs to read index/wavelet para from json str
 */
void 
par_read_json_source(cJSON *item, char *wavelet_type_name,
      char *src_name, float *src_coord, int *grid_index,
      char *wavelet_name, float *wavelet_coefs, float *t_start, float *t_end)
{
  cJSON *subitem;

  // event name
  if (subitem = cJSON_GetObjectItem(item, "name"))
  {
    sprintf(src_name,"%s",subitem->valuestring);
  }

  // if coordinate
  if (subitem = cJSON_GetObjectItem(item, "location_by_coords")) {
    for (int i = 0; i < FD_NDIM; i++) {
      src_coord[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
    }
  }

  // if grid index
  if (subitem = cJSON_GetObjectItem(item, "location_by_grid_index")) {
    for (int i = 0; i < FD_NDIM; i++) {
      grid_index[i] = cJSON_GetArrayItem(subitem, i)->valueint;
    }
  }

  // stf acts from start time
  if (subitem = cJSON_GetObjectItem(item, "start_time")) {
     *t_start = subitem->valuedouble;
  }
  // stf ends at end time
  if (subitem = cJSON_GetObjectItem(item, "end_time")) {
     *t_end = subitem->valuedouble;
  }

  // wavelet name
  if (subitem = cJSON_GetObjectItem(item, wavelet_type_name))
  {
    sprintf(wavelet_name,"%s",subitem->valuestring);
  }

  // coefs

  // ricker
  if (strcmp(wavelet_name, "ricker")==0) {
    if (subitem = cJSON_GetObjectItem(item, "ricker_center_frequency")) {
      wavelet_coefs[0] = subitem->valuedouble;
    }
    if (subitem = cJSON_GetObjectItem(item, "ricker_peak_time")) {
      wavelet_coefs[1] = subitem->valuedouble;
    }
  }

  // gaussian
  if (strcmp(wavelet_name, "gaussian")==0) {
    if (subitem = cJSON_GetObjectItem(item, "gaussian_rms_width")) {
      wavelet_coefs[0] = subitem->valuedouble;
    }
    if (subitem = cJSON_GetObjectItem(item, "gaussian_peak_time"))
    {
      wavelet_coefs[1] = subitem->valuedouble;
    }
  }
}

void
par_print(struct par_t *par)
{    
  fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> ESTIMATE MEMORY INFO.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "total memory size Byte: %20.5f  B\n", PSV->total_memory_size_Byte);
  //fprintf(stdout, "total memory size KB  : %20.5f KB\n", PSV->total_memory_size_KB  );
  //fprintf(stdout, "total memory size MB  : %20.5f MB\n", PSV->total_memory_size_MB  );
  //fprintf(stdout, "total memory size GB  : %20.5f GB\n", PSV->total_memory_size_GB  );
  //fprintf(stdout, "\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> FOLDER AND FILE INFO.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "   OutFolderName: %s\n", OutFolderName);
  //fprintf(stdout, "       EventName: %s\n", OutPrefix);
  //fprintf(stdout, "     LogFilename: %s\n", LogFilename);
  //fprintf(stdout, " StationFilename: %s\n", StationFilename);
  //fprintf(stdout, "  SourceFilename: %s\n", SourceFilename);
  //fprintf(stdout, "   MediaFilename: %s\n", MediaFilename);
  //fprintf(stdout, "\n");

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> MPI INFO:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " number_of_mpiprocs_x = %-10d\n", par->number_of_mpiprocs_x);
  fprintf(stdout, " number_of_mpiprocs_y = %-10d\n", par->number_of_mpiprocs_y);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> Time Integration INFO:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " size_of_time_step = %10.4e\n", par->size_of_time_step);
  fprintf(stdout, " number_of_time_steps = %-10d\n", par->number_of_time_steps);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> boundary layer information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " boundary_condition  = %10s%10s%10s%10s%10s%10s\n", 
          par->boundary_type_name[0],
          par->boundary_type_name[1],
          par->boundary_type_name[2],
          par->boundary_type_name[3],
          par->boundary_type_name[4],
          par->boundary_type_name[5]);
  fprintf(stdout, " abs_num_of_layers = %10d%10d%10d%10d%10d%10d\n", 
          par->abs_num_of_layers[0],
          par->abs_num_of_layers[1],
          par->abs_num_of_layers[2],
          par->abs_num_of_layers[3],
          par->abs_num_of_layers[4],
          par->abs_num_of_layers[5]);
  fprintf(stdout, " cfspml_velocity = %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n", 
          par->cfspml_velocity[0],
          par->cfspml_velocity[1],
          par->cfspml_velocity[2],
          par->cfspml_velocity[3],
          par->cfspml_velocity[4],
          par->cfspml_velocity[5]);
  fprintf(stdout, " cfspml_alpha_max = %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n", 
          par->cfspml_alpha_max[0],
          par->cfspml_alpha_max[1],
          par->cfspml_alpha_max[2],
          par->cfspml_alpha_max[3],
          par->cfspml_alpha_max[4],
          par->cfspml_alpha_max[5]);
  fprintf(stdout, " cfspml_beta_max = %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n", 
          par->cfspml_beta_max[0],
          par->cfspml_beta_max[1],
          par->cfspml_beta_max[2],
          par->cfspml_beta_max[3],
          par->cfspml_beta_max[4],
          par->cfspml_beta_max[5]);
  fprintf(stdout, "\n");

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> GRID INFO:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " grid_export_dir = %s\n", par->grid_export_dir);
  fprintf(stdout, " number_of_total_grid_points_x = %-10d\n", par->number_of_total_grid_points_x);
  fprintf(stdout, " number_of_total_grid_points_y = %-10d\n", par->number_of_total_grid_points_y);
  fprintf(stdout, " number_of_total_grid_points_z = %-10d\n", par->number_of_total_grid_points_z);

  fprintf(stdout, " grid_generation_itype = %d\n", par->grid_generation_itype);
  if (par->grid_generation_itype==PAR_GRID_CARTESIAN) {
    fprintf(stdout, " cartesian_grid_x0 = %10.4e\n", par->cartesian_grid_origin[0]);
    fprintf(stdout, " cartesian_grid_y0 = %10.4e\n", par->cartesian_grid_origin[1]);
    fprintf(stdout, " cartesian_grid_z0 = %10.4e\n", par->cartesian_grid_origin[2]);
    fprintf(stdout, " cartesian_grid_dx = %10.4e\n", par->cartesian_grid_stepsize[0]);
    fprintf(stdout, " cartesian_grid_dy = %10.4e\n", par->cartesian_grid_stepsize[1]);
    fprintf(stdout, " cartesian_grid_dz = %10.4e\n", par->cartesian_grid_stepsize[2]);
  }

  fprintf(stdout, " metric_method_itype = %d\n", par->metric_method_itype);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> media info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " media_export_dir = %s\n", par->media_export_dir);
  fprintf(stdout, " media_input_itype = %d\n", par->media_input_itype);
  //if (par->media_input_itype == PAR_MEDIA_3LAY) {
  //  fprintf()
  //}
  //fprintf(stdout, "\n --> the media filename is:\n");
  //fprintf(stdout, " velp_file  = %s\n", PSV->fnm_velp);
  //fprintf(stdout, " vels_file  = %s\n", PSV->fnm_vels);
  //fprintf(stdout, " rho_file   = %s\n", PSV->fnm_rho);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> source info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " source_input_itype = %d\n", par->source_input_itype);
  switch (par->source_input_itype)
  {
    case PAR_SOURCE_SINGLE_FORCE :
      fprintf(stdout, " source_input_type = %s\n", "single_force");
      fprintf(stdout, " source_name = %s\n", par->source_name);
      // grid index negative means to use coord
      if (par->source_gridindex[0] < 0) {
        fprintf(stdout, " source location = [%f, %f, %f]\n", 
            par->source_coords[0], par->source_coords[1], par->source_coords[2]);
      } else {
        fprintf(stdout, " source index = [%d, %d, %d]\n", 
            par->source_gridindex[0], par->source_gridindex[1], par->source_gridindex[2]);
      }
      fprintf(stdout, " source_time_functon = %s\n", par->wavelet_name);
      fprintf(stdout, " first two coefs = %f, %f\n", par->wavelet_coefs[0], par->wavelet_coefs[1]);
      fprintf(stdout, " start and end time = %f, %f\n", par->wavelet_tstart,par->wavelet_tend);
      fprintf(stdout, " force vector = [%f, %f, %f]\n", par->source_force_vector[0],
          par->source_force_vector[1], par->source_force_vector[2]);
      break;

    case PAR_SOURCE_SINGLE_MOMENT :
      fprintf(stdout, " source_input_type = %s\n", "single_moment");
      fprintf(stdout, " source_name = %s\n", par->source_name);
      if (par->source_gridindex[0] < 0) {
        fprintf(stdout, " source location = [%f, %f, %f]\n", 
            par->source_coords[0], par->source_coords[1], par->source_coords[2]);
      } else {
        fprintf(stdout, " source index = [%d, %d, %d]\n", 
            par->source_gridindex[0], par->source_gridindex[1], par->source_gridindex[2]);
      }
      fprintf(stdout, " moment_rate_functon = %s\n", par->wavelet_name);
      fprintf(stdout, " first two coefs = %f, %f\n", par->wavelet_coefs[0], par->wavelet_coefs[1]);
      fprintf(stdout, " start and end time = %f, %f\n", par->wavelet_tstart,par->wavelet_tend);
      fprintf(stdout, " moment tensor = [%f, %f, %f, %f, %f, %f]\n",
          par->source_moment_tensor[0],
          par->source_moment_tensor[1],
          par->source_moment_tensor[2],
          par->source_moment_tensor[3],
          par->source_moment_tensor[4],
          par->source_moment_tensor[5]);
      break;

    case PAR_SOURCE_FILE :
      fprintf(stdout, " source_input_type = %s\n", "in_source_file");
      fprintf(stdout, " in_source_file = %s\n", par->source_input_file);
      break;
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output_dir = %s\n", par->output_dir);

  fprintf(stdout, "--> station list file:\n");
  fprintf(stdout, " in_station_file = %s\n", par->in_station_file);

  fprintf(stdout, "--> recivers lines:\n");
  fprintf(stdout, "number_of_receiver_line=%d\n", par->number_of_receiver_line);
  if (par->number_of_receiver_line > 0)
  {
    fprintf(stdout, "#  name  i0  j0  k0   di  dj   dk  count\n");
    for(int n=0; n<par->number_of_receiver_line; n++)
    {
       fprintf(stdout, "%6d %s %6d %6d %6d %6d %6d %6d %6d\n",
           n,
           par->receiver_line_name[n],
           par->receiver_line_index_start[n*3+0],
           par->receiver_line_index_start[n*3+1],
           par->receiver_line_index_start[n*3+2],
           par->receiver_line_index_incre[n*3+0],
           par->receiver_line_index_incre[n*3+1],
           par->receiver_line_index_incre[n*3+2],
           par->receiver_line_count[n]);
    }
  }

  fprintf(stdout, "--> slice:\n");
  fprintf(stdout, "number_of_slice_x=%d: ", par->number_of_slice_x);
  for(int n=0; n<par->number_of_slice_x; n++)
  {
     fprintf(stdout, "%6d", par->slice_x_index[n]);
  }
  fprintf(stdout, "\nnumber_of_slice_y=%d: ", par->number_of_slice_y);
  for(int n=0; n<par->number_of_slice_y; n++)
  {
     fprintf(stdout, "%6d", par->slice_y_index[n]);
  }
  fprintf(stdout, "\nnumber_of_slice_z=%d: ", par->number_of_slice_z);
  for(int n=0; n<par->number_of_slice_z; n++)
  {
     fprintf(stdout, "%6d", par->slice_z_index[n]);
  }
  fprintf(stdout, "\n");

  fprintf(stdout, "--> snapshot information.\n");
  fprintf(stdout, "number_of_snapshot=%d\n", par->number_of_snapshot);
  if (par->number_of_snapshot > 0)
  {
      fprintf(stdout, "#  name  i0  j0  k0  ni  nj  nk  di  dj   dk  it0  nt  dit\n");
      for(int n=0; n<par->number_of_snapshot; n++)
      {
         fprintf(stdout, "%6d %s %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n",
             n,
             par->snapshot_name[n],
             par->snapshot_index_start[n*3+0],
             par->snapshot_index_start[n*3+1],
             par->snapshot_index_start[n*3+2],
             par->snapshot_index_count[n*3+0],
             par->snapshot_index_count[n*3+1],
             par->snapshot_index_count[n*3+2],
             par->snapshot_index_incre[n*3+0],
             par->snapshot_index_incre[n*3+1],
             par->snapshot_index_incre[n*3+2],
             par->snapshot_time_start[n],
             par->snapshot_time_incre[n]);
      }
  }

  fprintf(stdout, "--> qc parameters:\n");
  fprintf(stdout, "check_nan_every_nummber_of_steps=%d\n", par->check_nan_every_nummber_of_steps);
  fprintf(stdout, "output_all=%d\n", par->output_all);

  return;
}
