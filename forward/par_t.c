/*
 * 
 */

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "cJSON.h"
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

  // set non-input default values

  // read parameter
  cJSON *root = cJSON_Parse(str);
  if (NULL == root) {
    printf("Error at parsing json!\n");
    exit(-1);
  }

  cJSON *item;
  cJSON *subitem, *snapitem, *lineitem;

  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_x")) {
    par->number_of_total_grid_points_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_y")) {
    par->number_of_total_grid_points_y = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_total_grid_points_z")) {
    par->number_of_total_grid_points_z = item->valueint;
  }

  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_x")) {
    par->number_of_mpiprocs_x = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_mpiprocs_y")) {
    par->number_of_mpiprocs_y = item->valueint;
  }

  if (item = cJSON_GetObjectItem(root, "size_of_time_step")) {
    par->size_of_time_step = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "number_of_time_steps")) {
    par->number_of_time_steps = item->valueint;
  }

  par->time_start = 0.0;
  par->time_end   = par->time_start + par->number_of_time_steps * par->size_of_time_step;
  //int     nt_total = (int) ((par->time_end - par->time_start) / dt+0.5);

  // boundary
  if (item = cJSON_GetObjectItem(root, "boundary_condition"))
  {
    //int array_size = cJSON_GetArraySize(item);
    //for (int i = 0; i < FD_NDIM_2; i++) {
    //  sprintf(par->boundary_type_name[i], "%s",cJSON_GetArrayItem(item, i)->valuestring);
    //}
    if (subitem = cJSON_GetObjectItem(item, "x_left")) {
       sprintf(par->boundary_type_name[0], "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "x_right")) {
       sprintf(par->boundary_type_name[1], "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "y_front")) {
       sprintf(par->boundary_type_name[2], "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "y_back")) {
       sprintf(par->boundary_type_name[3], "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "z_bottom")) {
       sprintf(par->boundary_type_name[4], "%s", subitem->valuestring);
    }
    if (subitem = cJSON_GetObjectItem(item, "z_top")) {
       sprintf(par->boundary_type_name[5], "%s", subitem->valuestring);
    }
  }

  // cfs-pml
  if (item = cJSON_GetObjectItem(root, "cfspml"))
  {
    if (subitem = cJSON_GetObjectItem(item, "number_of_layers"))
    {
      for (int i = 0; i < FD_NDIM_2; i++) {
        par->abs_num_of_layers[i] = cJSON_GetArrayItem(subitem, i)->valueint;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "alpha_max"))
    {
      for (int i = 0; i < FD_NDIM_2; i++) {
        par->cfspml_alpha_max[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "beta_max"))
    {
      for (int i = 0; i < FD_NDIM_2; i++) {
        par->cfspml_beta_max[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
      }
    }
    if (subitem = cJSON_GetObjectItem(item, "pml_velocity"))
    {
      for (int i = 0; i < FD_NDIM_2; i++) {
        par->cfspml_velocity[i] = cJSON_GetArrayItem(subitem, i)->valuedouble;
      }
    }
  }

  //-- grid
  if (item = cJSON_GetObjectItem(root, "input_grid_type")) {
    sprintf(par->input_grid_type, "%s", item->valuestring);
  }
  //  if cartesian grid
  if (strcmp(par->input_grid_type, "cartesian")==0)
  {
    if (item = cJSON_GetObjectItem(root, "cartesian_grid_origin")) {
      for (int i = 0; i < FD_NDIM; i++) {
        par->cartesian_grid_origin[i] = cJSON_GetArrayItem(item, i)->valuedouble;
      }
    }
    if (item = cJSON_GetObjectItem(root, "cartesian_grid_stepsize")) {
      for (int i = 0; i < FD_NDIM; i++) {
        par->cartesian_grid_stepsize[i] = cJSON_GetArrayItem(item, i)->valuedouble;
      }
    }
  }
  //  if vmap
  if (strcmp(par->input_grid_type, "vmap")==0)
  {
    if (item = cJSON_GetObjectItem(root, "input_vmap_file")) {
        sprintf(par->input_grid_file,"%s",item->valuestring);
    }
  }

  //-- metric
  if (item = cJSON_GetObjectItem(root, "input_metric_type")) {
    sprintf(par->input_metric_type, "%s", item->valuestring);
  }
  //if (strcmp(par->input_metric_type, "calculate")==0)
  //{
  //  par->input_metric_itype = ;
  //}

  //-- medium
  if (item = cJSON_GetObjectItem(root, "input_medium_type")) {
    sprintf(par->input_medium_type, "%s", item->valuestring);
  }
  //if (strcmp(par->input_medium_type, "grid")==0)
  //{
  //  par->input_medium_itype = ;
  //}

  //-- source
  if (item = cJSON_GetObjectItem(root, "input_source_file")) {
    sprintf(par->input_source_file, "%s", item->valuestring);
  }

  //-- receiver
  if (item = cJSON_GetObjectItem(root, "input_receiver_file")) {
    sprintf(par->input_receiver_file, "%s", item->valuestring);
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

  //-- output dir
  if (item = cJSON_GetObjectItem(root, "output_dir")) {
      sprintf(par->output_dir,"%s",item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root, "grid_dir")) {
      sprintf(par->grid_dir,"%s",item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root, "media_dir")) {
      sprintf(par->media_dir,"%s",item->valuestring);
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
  fprintf(stdout, " grid_dir = %s\n", par->grid_dir);
  fprintf(stdout, " number_of_total_grid_points_x = %-10d\n", par->number_of_total_grid_points_x);
  fprintf(stdout, " number_of_total_grid_points_y = %-10d\n", par->number_of_total_grid_points_y);
  fprintf(stdout, " number_of_total_grid_points_z = %-10d\n", par->number_of_total_grid_points_z);

  fprintf(stdout, " input_grid_type = %s\n", par->input_grid_type);
  if (strcmp(par->input_grid_type, "cartesian")==0) {
    fprintf(stdout, " cartesian_grid_x0 = %10.4e\n", par->cartesian_grid_origin[0]);
    fprintf(stdout, " cartesian_grid_y0 = %10.4e\n", par->cartesian_grid_origin[1]);
    fprintf(stdout, " cartesian_grid_z0 = %10.4e\n", par->cartesian_grid_origin[2]);
    fprintf(stdout, " cartesian_grid_dx = %10.4e\n", par->cartesian_grid_stepsize[0]);
    fprintf(stdout, " cartesian_grid_dy = %10.4e\n", par->cartesian_grid_stepsize[1]);
    fprintf(stdout, " cartesian_grid_dz = %10.4e\n", par->cartesian_grid_stepsize[2]);
  }

  fprintf(stdout, " input_metric_type = %s\n", par->input_metric_type);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> media info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " media_dir = %s\n", par->media_dir);
  fprintf(stdout, " input_medium_type = %s\n", par->input_medium_type);

  //fprintf(stdout, "\n --> the media filename is:\n");
  //fprintf(stdout, " velp_file  = %s\n", PSV->fnm_velp);
  //fprintf(stdout, " vels_file  = %s\n", PSV->fnm_vels);
  //fprintf(stdout, " rho_file   = %s\n", PSV->fnm_rho);

  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> source info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " input_source_file = %s\n", par->input_source_file);

  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output_dir = %s\n", par->output_dir);

  fprintf(stdout, "--> individual recivers:\n");
  fprintf(stdout, " input_receiver_file = %s\n", par->input_receiver_file);

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

  /*
  fprintf(stdout, "\n");
  fprintf(stdout, "--> station information.\n");
  fprintf(stdout, " number_of_station  = %4d\n", par->number_of_station);
  fprintf(stdout, " seismo_format_sac  = %4d\n", par->seismo_format_sac );
  fprintf(stdout, " seismo_format_segy = %4d\n", par->seismo_format_segy);
  fprintf(stdout, " SeismoPrefix = %s\n", SeismoPrefix);
  fprintf(stdout, "\n");

  if(par->number_of_station > 0)
  {
      //fprintf(stdout, " station_indx:\n");
      fprintf(stdout, " stations             x           z           i           k:\n");
  }

  for(n=0; n<par->number_of_station; n++)
  {
      indx = 2*n;
      fprintf(stdout, "       %04d  %10.4e  %10.4e  %10d  %10d\n", n+1, 
              par->station_coord[indx], par->station_coord[indx+1],
              par->station_indx [indx], par->station_indx [indx+1]);
  }
  fprintf(stdout, "\n");
  */

  return;
}
