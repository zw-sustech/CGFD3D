/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "blk_t.h"

//
// malloc inner vars
//

int
blk_init(blk_t *blk,
         const int myid, const int verbose)
{
  int ierr = 0;

  blk->fd            = (fd_t *)malloc(sizeof(fd_t));
  blk->mympi         = (mympi_t *)malloc(sizeof(mympi_t));
  blk->gdinfo        = (gdinfo_t *)malloc(sizeof(gdinfo_t));
  blk->gdcurv        = (gdcurv_t         *)malloc(sizeof(gdcurv_t     ));
  blk->gdcurv_metric = (gdcurv_metric_t *)malloc(sizeof(gdcurv_metric_t));
  blk->mdeliso       = (mdeliso_t      *)malloc(sizeof(mdeliso_t     ));
  blk->wfel1st       = (wfel1st_t      *)malloc(sizeof(wfel1st_t     ));
  blk->src           = (src_t      *)malloc(sizeof(src_t     ));
  blk->bdryfree    = (bdryfree_t  *)malloc(sizeof(bdryfree_t ));
  blk->bdrypml     = (bdrypml_t  *)malloc(sizeof(bdrypml_t ));
  blk->iorecv      = (iorecv_t  *)malloc(sizeof(iorecv_t ));
  blk->ioline      = (ioline_t  *)malloc(sizeof(ioline_t ));
  blk->ioslice     = (ioslice_t  *)malloc(sizeof(ioslice_t ));
  blk->iosnap      = (iosnap_t  *)malloc(sizeof(iosnap_t ));

  sprintf(blk->name, "%s", "single");

  return ierr;
}

// set str
int
blk_set_output(blk_t *blk,
               mympi_t *mympi,
               char *output_dir,
               char *grid_export_dir,
               char *media_export_dir,
               const int verbose)
{
  // set name
  //sprintf(blk->name, "%s", name);

  // output name
  sprintf(blk->output_fname_part,"px%d_py%d", mympi->topoid[0],mympi->topoid[1]);

  // output
  sprintf(blk->output_dir, "%s", output_dir);
  sprintf(blk->grid_export_dir, "%s", grid_export_dir);
  sprintf(blk->media_export_dir, "%s", media_export_dir);

  return 0;
}

/*
 * mpi message
 */

void
blk_mesg_init(mympi_t *mympi,
                int ni,
                int nj,
                int nk,
                int fdx_nghosts,
                int fdy_nghosts,
                int num_of_vars)
{
  // mpi mesg
  int siz_x = (nj * nk * fdx_nghosts) * num_of_vars;
  int siz_y = (ni * nk * fdy_nghosts) * num_of_vars;

  int siz_buff =(  ni * nk * fdy_nghosts * 2
                    + nj * nk * fdx_nghosts * 2) * num_of_vars;
  mympi->sbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc sbuff");
  mympi->rbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc rbuff");
  for (int iptr=0; iptr < siz_buff; iptr++) {
    mympi->rbuff[iptr] = 0.0;
  }

  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + siz_x;
  float *sbuff_y1 = sbuff_x2 + siz_x;
  float *sbuff_y2 = sbuff_y1 + siz_y;

  int tag[4] = { 11, 12, 21, 22 };

  // send
  MPI_Send_init(sbuff_x1, siz_x, MPI_FLOAT, mympi->neighid[0], tag[0], mympi->topocomm, &(mympi->s_reqs[0]));
  MPI_Send_init(sbuff_x2, siz_x, MPI_FLOAT, mympi->neighid[1], tag[1], mympi->topocomm, &(mympi->s_reqs[1]));
  MPI_Send_init(sbuff_y1, siz_y, MPI_FLOAT, mympi->neighid[2], tag[2], mympi->topocomm, &(mympi->s_reqs[2]));
  MPI_Send_init(sbuff_y2, siz_y, MPI_FLOAT, mympi->neighid[3], tag[3], mympi->topocomm, &(mympi->s_reqs[3]));

  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + siz_x;
  float *rbuff_y1 = rbuff_x2 + siz_x;
  float *rbuff_y2 = rbuff_y1 + siz_y;

  // recv
  MPI_Recv_init(rbuff_x1, siz_x, MPI_FLOAT, mympi->neighid[0], tag[1], mympi->topocomm, &(mympi->r_reqs[0]));
  MPI_Recv_init(rbuff_x2, siz_x, MPI_FLOAT, mympi->neighid[1], tag[0], mympi->topocomm, &(mympi->r_reqs[1]));
  MPI_Recv_init(rbuff_y1, siz_y, MPI_FLOAT, mympi->neighid[2], tag[3], mympi->topocomm, &(mympi->r_reqs[2]));
  MPI_Recv_init(rbuff_y2, siz_y, MPI_FLOAT, mympi->neighid[3], tag[2], mympi->topocomm, &(mympi->r_reqs[3]));
}

//- for future: consider different left/right length
//void
//fd_blk_pack_mesg(float *restrict w_cur,
//                 int num_of_vars,
//                 int ni1, int ni2, int nj1, int nj2, int nk1, int nk2,
//                 size_t siz_line, size_t siz_slice, size_t siz_volume,
//                 int   *fdx_info,
//                 int   *fdy_info,
//                 int   *fdz_info,
//                 float *restrict buff)
//{
//  int npt_to_x1 = fdx_info[CONST_INFO_LENGTH_RIGTH];
//  int npt_to_x2 = fdx_info[CONST_INFO_LENGTH_LEFT ];
//
//  size_t iptr_b = 0;
//
//  for (int ivar=0; ivar<num_of_vars; ivar++)
//  {
//    size_t iptr_var = ivar * siz_volume;
//
//    // x1
//    for (int k=nk1; k<=nk2; k++)
//    {
//      size_t iptr_k = iptr_var + k * siz_slice;
//
//      for (int j=nj1; j<=nj2; j++)
//      {
//        size_t iptr_j = iptr_k + j * siz_line;
//
//        for (int i=ni1; i<ni1+npt_to_x1; i++)
//        {
//          buff[iptr_b] = w_cur[iptr_j + i];
//          iptr_b++;
//        }
//      }
//    } 
//
//    // y1y2
//
//  } // ivar
//}

void
blk_pack_mesg(float *restrict w_cur,float *restrict sbuff,
                 int num_of_vars, gdinfo_t *gdinfo,
                 int   fdx_nghosts, int   fdy_nghosts)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;
  size_t siz_volume = gdinfo->siz_volume;

  size_t iptr_b = 0;

  // x1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<ni1+fdx_nghosts; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // x2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni2-fdx_nghosts+1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // y1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<nj1+fdy_nghosts; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // y2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj2-fdy_nghosts+1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  } // ivar

}

void
blk_unpack_mesg(float *restrict rbuff,float *restrict w_cur,
                 int num_of_vars, gdinfo_t *gdinfo,
                 int   fdx_nghosts, int   fdy_nghosts)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;
  size_t siz_volume = gdinfo->siz_volume;

  size_t iptr_b = 0;

  // x1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1-fdx_nghosts; i<ni1; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  }

  // x2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni2+1; i<=ni2+fdx_nghosts; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  }

  // y1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1-fdy_nghosts; j<nj1; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  }

  // y2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj2+1; j<=nj2+fdy_nghosts; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<=ni2; i++)
        {
          w_cur[iptr_j + i] = rbuff[iptr_b];
          iptr_b++;
        }
      }
    }
  } // ivar
}

int
blk_print(blk_t *blk)
{    
  int ierr = 0;

  fprintf(stdout, "\n-------------------------------------------------------\n");
  fprintf(stdout, "print blk %s:\n", blk->name);
  fprintf(stdout, "-------------------------------------------------------\n\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
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

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> media info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //if (blk->media_type == MEDIA_TYPE_LAYER)
  //{
  //    strcpy(str, "layer");
  //}
  //else if (blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    strcpy(str, "grid");
  //}
  //fprintf(stdout, " media_type = %s\n", str);
  //if(blk->media_type == MEDIA_TYPE_GRID)
  //{
  //    fprintf(stdout, "\n --> the media filename is:\n");
  //    fprintf(stdout, " velp_file  = %s\n", blk->fnm_velp);
  //    fprintf(stdout, " vels_file  = %s\n", blk->fnm_vels);
  //    fprintf(stdout, " rho_file   = %s\n", blk->fnm_rho);
  //}
  //fprintf(stdout, "\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> source info.\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, " number_of_force  = %d\n", blk->number_of_force);
  //if(blk->number_of_force > 0)
  //{
  //    fprintf(stdout, " force_source           x           z     x_shift     z_shift           i           k:\n");
  //    for(n=0; n<blk->number_of_force; n++)
  //    {
  //        indx = 2*n;
  //        fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //                blk->force_coord[indx], blk->force_coord[indx+1],
  //                blk->force_shift[indx], blk->force_shift[indx+1],
  //                blk->force_indx [indx], blk->force_indx [indx+1]);
  //    }
  //    fprintf(stdout, "\n");
  //}

  //fprintf(stdout, "\n");
  //fprintf(stdout, " number_of_moment = %d\n", blk->number_of_moment);
  //if(blk->number_of_moment > 0)
  //{
  //    fprintf(stdout, " moment_source          x           z     x_shift     z_shift           i           k:\n");
  //    for(n=0; n<blk->number_of_moment; n++)
  //    {
  //        indx = 2*n;
  //        fprintf(stdout, "         %04d  %10.4e  %10.4e  %10.4e  %10.4e  %10d  %10d\n", n+1, 
  //                blk->moment_coord[indx], blk->moment_coord[indx+1],
  //                blk->moment_shift[indx], blk->moment_shift[indx+1],
  //                blk->moment_indx [indx], blk->moment_indx [indx+1]);
  //    }
  //    fprintf(stdout, "\n");
  //}

  //fprintf(stdout, "-------------------------------------------------------\n");
  //fprintf(stdout, "--> boundary layer information:\n");
  //fprintf(stdout, "-------------------------------------------------------\n");
  //ierr = boundary_id2type(type1, blk->boundary_type[0], errorMsg);
  //ierr = boundary_id2type(type2, blk->boundary_type[1], errorMsg);
  //ierr = boundary_id2type(type3, blk->boundary_type[2], errorMsg);
  //ierr = boundary_id2type(type4, blk->boundary_type[3], errorMsg);
  //fprintf(stdout, " boundary_type         = %10s%10s%10s%10s\n", 
  //        type1, type2, type3, type4);
  //fprintf(stdout, " boundary_layer_number = %10d%10d%10d%10d\n", 
  //        blk->boundary_layer_number[0], blk->boundary_layer_number[1], 
  //        blk->boundary_layer_number[2], blk->boundary_layer_number[3]);
  //fprintf(stdout, "\n");
  //fprintf(stdout, " absorb_velocity       = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->absorb_velocity[0], blk->absorb_velocity[1], blk->absorb_velocity[2], 
  //        blk->absorb_velocity[3]);
  //fprintf(stdout, "\n");
  //fprintf(stdout, " CFS_alpha_max         = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->CFS_alpha_max[0], blk->CFS_alpha_max[1], blk->CFS_alpha_max[2], 
  //        blk->CFS_alpha_max[3]);
  //fprintf(stdout, " CFS_beta_max          = %10.2f%10.2f%10.2f%10.2f\n", 
  //        blk->CFS_beta_max[0], blk->CFS_beta_max[1], blk->CFS_beta_max[2], 
  //        blk->CFS_beta_max[3]);
  
  return ierr;
}
