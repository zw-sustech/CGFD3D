/*
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "fdlib_math.h"
#include "fdlib_mem.h"
#include "fd_t.h"
#include "io_funcs.h"

/*
 * export to a single file
 */

void
io_snapshot_export(char *fname,
                   float *restrict var,
                   int nx,
                   int ny,
                   int nz,
                   int *snap_indx,
                   int verbose)
{
  FILE *fp=fopen(fname,"wb");
  if (fp == NULL) {
    fprintf(stderr,"Error: can't create : %s\n", fname);
    exit(1);
  }

  // index triple
  int i1 = snap_indx[0];
  int j1 = snap_indx[1];
  int k1 = snap_indx[2];
  int ic = snap_indx[3];
  int jc = snap_indx[4];
  int kc = snap_indx[5];
  int di = snap_indx[6];
  int dj = snap_indx[7];
  int dk = snap_indx[8];

  for (int n3=0; n3<kc; n3++)
  {
    int k = k1 + n3 * dk;
    for (int n2=0; n2<jc; n2++)
    {
      int j = j1 + n2 * dj;
      for (int n1=0; n1<ic; n1++)
      {
        int i = i1 + n1 * di;
        int iptr = i + j * nx + k * nx * ny;
        fwrite(var+iptr, 1, sizeof(float), fp);
      }
    }
  }

  fclose(fp);
}

/*
 * append to a single file
 */

/*
void
io_snapshot_append(FILE *fp,
                   float *restrict var,
                   size_t ni,
                   size_t nj,
                   size_t nk,
                   int verbose)
{

}
*/
