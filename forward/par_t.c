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
  if (myid==0)
  {
    FILE fp=fopen(par_fname,"r");
    if (!fp) {
      fprintf(stderr,"Error: can't open par file: %s\n", par_fname);
      MPI_Finalize();
      exit(1);
    }
    fseek(fp, 0, SEEK_END);
    long len = ftell(fp);

    // bcast len to all nodes
    MPI_Bcast(len, 1, MPI_LONG, 0, comm);

    fseek(fp, 0, SEEK_SET);
    char *str = (char*)malloc(len+1);
    fread(str, 1, len, fp);
    fclose(fp);

    // bcast str
    MPI_Bcast(str, len+1, MPI_CHAR, 0, comm);
  }
  else
  {
    long len;
    // get len 
    MPI_Bcast(len, 1, MPI_LONG, 0, comm);

    char *str = (char*)malloc(len+1);
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
  FILE *fp = fopen(par_file_name,"r");

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
par_read_from_str(char *str, struct par_struct *par)
{
  // allocate
  par->boundary_type_name = (char *)malloc(FD_NDIM_2 * sizeof(char*));
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

  // grid
  if (item = cJSON_GetObjectItem(root, "coord_by_cartesian")) {
    par->coord_by_cartesian = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root, "cartesian_grid_x0")) {
    par->cartesian_grid_x0 = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "cartesian_grid_y0")) {
    par->cartesian_grid_y0 = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "cartesian_grid_z0")) {
    par->cartesian_grid_z0 = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "cartesian_grid_dx")) {
    par->cartesian_grid_dx = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "cartesian_grid_dy")) {
    par->cartesian_grid_dy = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "cartesian_grid_dz")) {
    par->cartesian_grid_dz = item->valuedouble;
  }
  if (item = cJSON_GetObjectItem(root, "metric_by_import")) {
    par->metric_by_import = item->valueint;
  }

  if (item = cJSON_GetObjectItem(root, "medium_by_import")) {
    par->medium_by_import = item->valueint;
  }

  // boundary
  if (item = cJSON_GetObjectItem(root, "boundary_condition"))
  {
    //int array_size = cJSON_GetArraySize(item);
    for (int i = 0; i < FD_NDIM_2; i++) {
      strncpy(par->boundary_type_name[i],
              cJSON_GetArrayItem(item, i)->valuestring,
              sizeof(cJSON_GetArrayItem(item, i)->valuestring));
    }
  }

  // pml
  if (item = cJSON_GetObjectItem(root, "abs_num_of_layers"))
  {
    for (int i = 0; i < FD_NDIM_2; i++) {
      par->abs_num_of_layers[i] = cJSON_GetArrayItem(item, i)->valueint;
    }
  }
  if (item = cJSON_GetObjectItem(root, "cfspml_alpha_max"))
  {
    for (int i = 0; i < FD_NDIM_2; i++) {
      par->cfspml_alpha_max[i] = cJSON_GetArrayItem(item, i)->valuedouble;
    }
  }
  if (item = cJSON_GetObjectItem(root, "cfspml_beta_max"))
  {
    for (int i = 0; i < FD_NDIM_2; i++) {
      par->cfspml_beta_max[i] = cJSON_GetArrayItem(item, i)->valuedouble;
    }
  }
  if (item = cJSON_GetObjectItem(root, "cfspml_velocity"))
  {
    for (int i = 0; i < FD_NDIM_2; i++) {
      par->cfspml_velocity[i] = cJSON_GetArrayItem(item, i)->valuedouble;
    }
  }

  /*
  if (item = cJSON_GetObjectItem(root, "OUT"))
      memcpy(OUT, item->valuestring, strlen(item->valuestring));

  if ( item = cJSON_GetObjectItem(root, "NUM_INTER") ) {
      par->nx = item->valueint;
  } else {
      error_exit("no nx in par file");
  }
  NX    = cJSON_GetObjectItem(root, "NX")->valueint;
  TMAX  = cJSON_GetObjectItem(root, "TMAX")->valuedouble;
  DT    = cJSON_GetObjectItem(root, "DT")->valuedouble;
  if(p = cJSON_GetObjectItem(root, "OUT" )) {
    strcpy(OUT, p->valuestring);
  }
  */

  cJSON_Delete(root);

  //mkpath(OUT, 0700);
  //mkpath(OUT_DIR, 0700);

  //printf("---------------------------------------------\n");
  //printf("TSKP  : %d\n", TSKP);
  //printf("TMAX  : %f\n", TMAX);
  //printf("DT  : %f\n", DT);
  //printf("DH  : %f\n", DH);
  //printf("NX  : %d\n", NX);
  //printf("NY  : %d\n", NY);
  //printf("NZ  : %d\n", NZ);
  //printf("PX  : %d\n", PX);
  //printf("PY  : %d\n", PY);
  //printf("PZ  : %d\n", PZ);
  //printf("---------------------------------------------\n");

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
  fprintf(stdout, "--> GRID INFO:\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " number_of_total_grid_points_x = %-10d\n", par->number_of_total_grid_points_x);
  fprintf(stdout, " number_of_total_grid_points_y = %-10d\n", par->number_of_total_grid_points_y);
  fprintf(stdout, " number_of_total_grid_points_z = %-10d\n", par->number_of_total_grid_points_z);
  fprintf(stdout, " coord_by_cartesian = %-10d\n", par->coord_by_cartesian);

  fprintf(stdout, " cartesian_grid_x0 = %10.4e\n", par->cartesian_grid_x0);
  fprintf(stdout, " cartesian_grid_y0 = %10.4e\n", par->cartesian_grid_y0);
  fprintf(stdout, " cartesian_grid_z0 = %10.4e\n", par->cartesian_grid_z0);
  fprintf(stdout, " cartesian_grid_dx = %10.4e\n", par->cartesian_grid_dx);
  fprintf(stdout, " cartesian_grid_dy = %10.4e\n", par->cartesian_grid_dy);
  fprintf(stdout, " cartesian_grid_dz = %10.4e\n", par->cartesian_grid_dz);


  /*
  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> media info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  if (PSV->media_type == MEDIA_TYPE_LAYER)
  {
      strcpy(str, "layer");
  }
  else if (PSV->media_type == MEDIA_TYPE_GRID)
  {
      strcpy(str, "grid");
  }
  fprintf(stdout, " media_type = %s\n", str);
  if(PSV->media_type == MEDIA_TYPE_GRID)
  {
      fprintf(stdout, "\n --> the media filename is:\n");
      fprintf(stdout, " velp_file  = %s\n", PSV->fnm_velp);
      fprintf(stdout, " vels_file  = %s\n", PSV->fnm_vels);
      fprintf(stdout, " rho_file   = %s\n", PSV->fnm_rho);
  }
  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> source info.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, " number_of_force  = %d\n", PSV->number_of_force);
  if(PSV->number_of_force > 0)
  {
      fprintf(stdout, " force_source           x           z     x_shift     z_shift           i           k:\n");
      for(n=0; n<PSV->number_of_force; n++)
      {
          indx = 2*n;
          fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
                  PSV->force_coord[indx], PSV->force_coord[indx+1],
                  PSV->force_shift[indx], PSV->force_shift[indx+1],
                  PSV->force_indx [indx], PSV->force_indx [indx+1]);
      }
      fprintf(stdout, "\n");
  }

  fprintf(stdout, "\n");
  fprintf(stdout, " number_of_moment = %d\n", PSV->number_of_moment);
  if(PSV->number_of_moment > 0)
  {
      fprintf(stdout, " moment_source          x           z     x_shift     z_shift           i           k:\n");
      for(n=0; n<PSV->number_of_moment; n++)
      {
          indx = 2*n;
          fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
                  PSV->moment_coord[indx], PSV->moment_coord[indx+1],
                  PSV->moment_shift[indx], PSV->moment_shift[indx+1],
                  PSV->moment_indx [indx], PSV->moment_indx [indx+1]);
      }
      fprintf(stdout, "\n");
  }
  */

  fprintf(stdout, "\n");
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
  fprintf(stdout, " abs_velocity = %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n", 
          par->abs_velocity[0],
          par->abs_velocity[1],
          par->abs_velocity[2],
          par->abs_velocity[3],
          par->abs_velocity[4],
          par->abs_velocity[5]);
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
  
  /*
  fprintf(stdout, "\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> output information.\n");
  fprintf(stdout, "-------------------------------------------------------\n");
  fprintf(stdout, "--> snapshot information.\n");
  if (PSV->number_of_snapshot > 0)
  {
      fprintf(stdout, "number_of_snapshot=%d\n", PSV->number_of_snapshot);
      fprintf(stdout, "#   x0    z0    nx    nz    dx    dz    dt     tdim_max    component\n");
      for(n=0; n<PSV->number_of_snapshot; n++)
      {
          indx = 10*n;
          componentV = ' ';
          componentT = ' ';

          if(PSV->snapshot_information[indx+8] == 1) {
              componentV = 'V';
          }

          if(PSV->snapshot_information[indx+9] == 1) {
              componentT = 'T';
          }

          for(i=0; i<7; i++)
          {
              fprintf(stdout, "%6d", PSV->snapshot_information[indx+i]);
          }

          fprintf(stdout, "%12d", PSV->snapshot_information[indx+7]);
          fprintf(stdout, "         %c%c\n", componentV, componentT);
      }
  }

  fprintf(stdout, "\n");
  fprintf(stdout, "--> station information.\n");
  fprintf(stdout, " number_of_station  = %4d\n", PSV->number_of_station);
  fprintf(stdout, " seismo_format_sac  = %4d\n", PSV->seismo_format_sac );
  fprintf(stdout, " seismo_format_segy = %4d\n", PSV->seismo_format_segy);
  fprintf(stdout, " SeismoPrefix = %s\n", SeismoPrefix);
  fprintf(stdout, "\n");

  if(PSV->number_of_station > 0)
  {
      //fprintf(stdout, " station_indx:\n");
      fprintf(stdout, " stations             x           z           i           k:\n");
  }

  for(n=0; n<PSV->number_of_station; n++)
  {
      indx = 2*n;
      fprintf(stdout, "       %04d  %10.4e  %10.4e  %10d  %10d\n", n+1, 
              PSV->station_coord[indx], PSV->station_coord[indx+1],
              PSV->station_indx [indx], PSV->station_indx [indx+1]);
  }
  fprintf(stdout, "\n");
  */

  return;
}
