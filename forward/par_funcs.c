/*
 * ref to set_params.cu
 */

//#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cJSON.h"
#include "par_funcs.h"

int par_read_file(char *par_file_name, struct par_struct *par)
{

  // set non-input default values
  par->number_of_points = ;

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

  //printf("configure json file:\n%s\n", str);

  //
  // read each parameter
  //
  //
  cJSON *root, *item;
  root = cJSON_Parse(str);
  if (NULL == root) {
    printf("Error at parsing json!\n");
    exit(-1);
  }


if (item = cJSON_GetObjectItem(root, "OUT"))
    memcpy(OUT, item->valuestring, strlen(item->valuestring));

if (item = cJSON_GetObjectItem(root, "Fault_grid")){
  //int array_size = cJSON_GetArraySize(item);
    for (int i = 0; i < 4; i++){
      hostParams.Fault_grid[i] = cJSON_GetArrayItem(item, i)->valueint;
    }
  }

  if (item = cJSON_GetObjectItem(root, "TMAX"))
    hostParams.TMAX  = item->valuedouble;
  if (item = cJSON_GetObjectItem(root, "DT"))
    hostParams.DT = item->valuedouble;
  if (item = cJSON_GetObjectItem(root, "DH"))
    hostParams.DH = item->valuedouble;

  if (item = cJSON_GetObjectItem(root, "NX"))
    hostParams.NX = item->valueint;
  if (item = cJSON_GetObjectItem(root, "NY"))
    hostParams.NY = item->valueint;
  if (item = cJSON_GetObjectItem(root, "NZ"))
    hostParams.NZ = item->valueint;
  if (item = cJSON_GetObjectItem(root, "PX"))
    hostParams.PX = item->valueint;
  if (item = cJSON_GetObjectItem(root, "PY"))
    hostParams.PY = item->valueint;
  if (item = cJSON_GetObjectItem(root, "PZ"))
    hostParams.PZ = item->valueint;

  item = cJSON_GetObjectItem(root, "IN_SOURCE");
  memcpy(IN_SOURCE, item->valuestring, strlen(item->valuestring));

  item = cJSON_GetObjectItem(root, "IN_TOPO");
  memcpy(IN_TOPO, item->valuestring, strlen(item->valuestring));

  item = cJSON_GetObjectItem(root, "IN_VELO");
  memcpy(IN_VELO, item->valuestring, strlen(item->valuestring));

  item = cJSON_GetObjectItem(root, "IN_INTER");
  memcpy(IN_INTER, item->valuestring, strlen(item->valuestring));

  item = cJSON_GetObjectItem(root, "IN_REC");
  memcpy(IN_REC, item->valuestring, strlen(item->valuestring));

  item = cJSON_GetObjectItem(root, "OUT_DIR");
  memcpy(OUT_DIR, item->valuestring, strlen(item->valuestring));

  //name_json = cJSON_GetObjectItem(root, "INGRD");
  //if (NULL != name_json)
  //{
  //    INGRD = cJSON_Print(name_json);
  //    printf("INGRD : %s\n", INGRD);
  //    free(name_json);
  //}

  //name_json = cJSON_GetObjectItem(root, "INVEL");
  //if (NULL != name_json)
  //{
  //    INVEL = cJSON_Print(name_json);
  //    printf("INVEL : %s\n", INVEL);
  //    free(name_json);
  //}

  //name_json = cJSON_GetObjectItem(root, "INSRC");
  //if (NULL != name_json)
  //{
  //    INSRC = cJSON_Print(name_json);
  //    printf("INSRC : %s\n", INSRC);
  //    free(name_json);
  //}

  //name_json = cJSON_GetObjectItem(root, "INREC");
  //if (NULL != name_json)
  //{
  //    INREC = cJSON_Print(name_json);
  //    printf("INREC : %s\n", INREC);
  //    free(name_json);
  //}

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

  cJSON_Delete(root);

  mkpath(OUT, 0700);
  mkpath(OUT_DIR, 0700);

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

void print_params() {
  printf( "# Parameters setting:\n"
          "* Tmax           = %g (sec)\n"
          "* DT             = %g (sec)\n"
          "* DH             = %g (m)\n"
          "* NX, NY, NZ, NT = %d %d %d %d\n"
          "* PX, PY, PZ     = %d %d %d\n"
          "* Fault_grid     = %d %d %d %d\n"
          "* Asp_grid       = %d %d %d %d\n"
          "* Dc             = %g (m)\n"
          "* C0             = %g (Pa)\n"
          "* mu_s           = %g\n"
          "* mu_d           = %g\n"
          "* vp             = %g/%g\n"
          "* vs             = %g/%g\n"
          "* rho            = %g/%g\n"
          "* PML_N          = %d\n"
          "* PML_velocity   = %g\n"
          "* PML_fc         = %g\n"
          "* PML_bmax       = %g\n"
          "\n"
          ,
          hostParams.TMAX,
          hostParams.DT,
          hostParams.DH,
          hostParams.NX,
          hostParams.NY,
          hostParams.NZ,
          hostParams.NT,
          hostParams.PX,
          hostParams.PY,
          hostParams.PZ,
          hostParams.Fault_grid[0],
          hostParams.Fault_grid[1],
          hostParams.Fault_grid[2],
          hostParams.Fault_grid[3],
          hostParams.Asp_grid[0],
          hostParams.Asp_grid[1],
          hostParams.Asp_grid[2],
          hostParams.Asp_grid[3],
          hostParams.Dc,
          hostParams.C0,
          hostParams.mu_s,
          hostParams.mu_d,
          hostParams.bi_vp1,
          hostParams.bi_vp2,
          hostParams.bi_vs1,
          hostParams.bi_vs2,
          hostParams.bi_rho1,
          hostParams.bi_rho2,
          hostParams.PML_N,
          hostParams.PML_velocity,
          hostParams.PML_fc,
          hostParams.PML_bmax
          );
  return;
}
