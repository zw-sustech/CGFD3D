/*********************************************************************
 * setup fd operators
 **********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fdlib_mem.h"
#include "fdlib_math.h"
#include "blk_t.h"

//
// malloc inner vars
//

int
blk_init(blk_t *blk,
         const int myid, const int verbose)
{
  int ierr = 0;

  // alloc struct vars
  blk->fd            = (fd_t *)malloc(sizeof(fd_t));
  blk->fdstg         = (fdstg_t *)malloc(sizeof(fdstg_t));
  blk->mympi         = (mympi_t *)malloc(sizeof(mympi_t));
  blk->gdinfo        = (gdinfo_t *)malloc(sizeof(gdinfo_t));
  blk->gd            = (gd_t        *)malloc(sizeof(gd_t     ));
  blk->gdcurv_metric = (gdcurv_metric_t *)malloc(sizeof(gdcurv_metric_t));
  blk->md            = (md_t      *)malloc(sizeof(md_t     ));
  blk->wav           = (wav_t      *)malloc(sizeof(wav_t     ));
  blk->src           = (src_t      *)malloc(sizeof(src_t     ));
  blk->bdryfree      = (bdryfree_t *)malloc(sizeof(bdryfree_t ));
  blk->bdrypml       = (bdrypml_t  *)malloc(sizeof(bdrypml_t ));
  blk->iorecv        = (iorecv_t   *)malloc(sizeof(iorecv_t ));
  blk->ioline        = (ioline_t   *)malloc(sizeof(ioline_t ));
  blk->ioslice       = (ioslice_t  *)malloc(sizeof(ioslice_t ));
  blk->iosnap        = (iosnap_t   *)malloc(sizeof(iosnap_t ));

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

/*********************************************************************
 * mpi message for collocated grid and center differences
 *********************************************************************/

void
blk_colcent_mesg_init(mympi_t *mympi,
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
  mympi->siz_sbuff = siz_buff;
  mympi->siz_rbuff = siz_buff;

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

void
blk_colcent_pack_mesg(float *restrict w_cur,float *restrict sbuff,
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
blk_colcent_unpack_mesg(float *restrict rbuff,float *restrict w_cur,
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

  return;
}

/*********************************************************************
 * mpi message for macdrp scheme with rk
 *********************************************************************/

void
blk_macdrp_mesg_init(mympi_t *mympi,
                fd_t *fd,
                int ni,
                int nj,
                int nk,
                int num_of_vars)
{
  // alloc
  mympi->pair_siz_sbuff_x1 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_sbuff_x2 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_sbuff_y1 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_sbuff_y2 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_rbuff_x1 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_rbuff_x2 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_rbuff_y1 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_siz_rbuff_y2 = (size_t **)malloc(fd->num_of_pairs * sizeof(size_t *));
  mympi->pair_s_reqs       = (MPI_Request ***)malloc(fd->num_of_pairs * sizeof(MPI_Request **));
  mympi->pair_r_reqs       = (MPI_Request ***)malloc(fd->num_of_pairs * sizeof(MPI_Request **));
  for (int ipair = 0; ipair < fd->num_of_pairs; ipair++)
  {
    mympi->pair_siz_sbuff_x1[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_sbuff_x2[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_sbuff_y1[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_sbuff_y2[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_rbuff_x1[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_rbuff_x2[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_rbuff_y1[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_siz_rbuff_y2[ipair] = (size_t *)malloc(fd->num_rk_stages * sizeof(size_t));
    mympi->pair_s_reqs[ipair] = (MPI_Request **)malloc(fd->num_rk_stages * sizeof(MPI_Request *));
    mympi->pair_r_reqs[ipair] = (MPI_Request **)malloc(fd->num_rk_stages * sizeof(MPI_Request *));

    for (int istage = 0; istage < fd->num_rk_stages; istage++)
    {
      mympi->pair_s_reqs[ipair][istage] = (MPI_Request *)malloc(4 * sizeof(MPI_Request));
      mympi->pair_r_reqs[ipair][istage] = (MPI_Request *)malloc(4 * sizeof(MPI_Request));
    }
  }

  // mpi mesg
  mympi->siz_sbuff = 0;
  mympi->siz_rbuff = 0;
  for (int ipair = 0; ipair < fd->num_of_pairs; ipair++)
  {
    for (int istage = 0; istage < fd->num_rk_stages; istage++)
    {
      fd_op_t *fdx_op = fd->pair_fdx_op[ipair][istage];
      fd_op_t *fdy_op = fd->pair_fdy_op[ipair][istage];

      // to x1, depends on right_len of x1 proc
      mympi->pair_siz_sbuff_x1[ipair][istage] = (nj * nk * fdx_op->right_len) * num_of_vars;
      // to x2, depends on left_len of x2 proc
      mympi->pair_siz_sbuff_x2[ipair][istage] = (nj * nk * fdx_op->left_len ) * num_of_vars;

      mympi->pair_siz_sbuff_y1[ipair][istage] = (ni * nk * fdy_op->right_len) * num_of_vars;
      mympi->pair_siz_sbuff_y2[ipair][istage] = (ni * nk * fdy_op->left_len ) * num_of_vars;

      // from x1, depends on left_len of cur proc
      mympi->pair_siz_rbuff_x1[ipair][istage] = (nj * nk * fdx_op->left_len ) * num_of_vars;
      // from x2, depends on right_len of cur proc
      mympi->pair_siz_rbuff_x2[ipair][istage] = (nj * nk * fdx_op->right_len) * num_of_vars;

      mympi->pair_siz_rbuff_y1[ipair][istage] = (ni * nk * fdy_op->left_len ) * num_of_vars;
      mympi->pair_siz_rbuff_y2[ipair][istage] = (ni * nk * fdy_op->right_len) * num_of_vars;

      size_t siz_s =  mympi->pair_siz_sbuff_x1[ipair][istage]
                    + mympi->pair_siz_sbuff_x2[ipair][istage]
                    + mympi->pair_siz_sbuff_y1[ipair][istage]
                    + mympi->pair_siz_sbuff_y2[ipair][istage];
      size_t siz_r =  mympi->pair_siz_rbuff_x1[ipair][istage]
                    + mympi->pair_siz_rbuff_x2[ipair][istage]
                    + mympi->pair_siz_rbuff_y1[ipair][istage]
                    + mympi->pair_siz_rbuff_y2[ipair][istage];

      if (siz_s > mympi->siz_sbuff) mympi->siz_sbuff = siz_s;
      if (siz_r > mympi->siz_rbuff) mympi->siz_rbuff = siz_r;
    }
  }

  mympi->sbuff = (float *) fdlib_mem_malloc_1d(mympi->siz_sbuff * sizeof(MPI_FLOAT),"alloc sbuff");
  mympi->rbuff = (float *) fdlib_mem_malloc_1d(mympi->siz_rbuff * sizeof(MPI_FLOAT),"alloc rbuff");

  // set up pers communication
  for (int ipair = 0; ipair < fd->num_of_pairs; ipair++)
  {
    for (int istage = 0; istage < fd->num_rk_stages; istage++)
    {
      size_t siz_s_x1 = mympi->pair_siz_sbuff_x1[ipair][istage];
      size_t siz_s_x2 = mympi->pair_siz_sbuff_x2[ipair][istage];
      size_t siz_s_y1 = mympi->pair_siz_sbuff_y1[ipair][istage];
      size_t siz_s_y2 = mympi->pair_siz_sbuff_y2[ipair][istage];

      float *sbuff_x1 = mympi->sbuff;
      float *sbuff_x2 = sbuff_x1 + siz_s_x1;
      float *sbuff_y1 = sbuff_x2 + siz_s_x2;
      float *sbuff_y2 = sbuff_y1 + siz_s_y1;

      // npair: xx, nstage: x, 
      int tag_pair_stage = ipair * 1000 + istage * 100;
      int tag[4] = { tag_pair_stage+11, tag_pair_stage+12, tag_pair_stage+21, tag_pair_stage+22 };

      // send
      MPI_Send_init(sbuff_x1, siz_s_x1, MPI_FLOAT, mympi->neighid[0], tag[0], mympi->topocomm, &(mympi->pair_s_reqs[ipair][istage][0]));
      MPI_Send_init(sbuff_x2, siz_s_x2, MPI_FLOAT, mympi->neighid[1], tag[1], mympi->topocomm, &(mympi->pair_s_reqs[ipair][istage][1]));
      MPI_Send_init(sbuff_y1, siz_s_y1, MPI_FLOAT, mympi->neighid[2], tag[2], mympi->topocomm, &(mympi->pair_s_reqs[ipair][istage][2]));
      MPI_Send_init(sbuff_y2, siz_s_y2, MPI_FLOAT, mympi->neighid[3], tag[3], mympi->topocomm, &(mympi->pair_s_reqs[ipair][istage][3]));

      // recv
      size_t siz_r_x1 = mympi->pair_siz_rbuff_x1[ipair][istage];
      size_t siz_r_x2 = mympi->pair_siz_rbuff_x2[ipair][istage];
      size_t siz_r_y1 = mympi->pair_siz_rbuff_y1[ipair][istage];
      size_t siz_r_y2 = mympi->pair_siz_rbuff_y2[ipair][istage];

      float *rbuff_x1 = mympi->rbuff;
      float *rbuff_x2 = rbuff_x1 + siz_r_x1;
      float *rbuff_y1 = rbuff_x2 + siz_r_x2;
      float *rbuff_y2 = rbuff_y1 + siz_r_y1;

      // recv
      MPI_Recv_init(rbuff_x1, siz_r_x1, MPI_FLOAT, mympi->neighid[0], tag[1], mympi->topocomm, &(mympi->pair_r_reqs[ipair][istage][0]));
      MPI_Recv_init(rbuff_x2, siz_r_x2, MPI_FLOAT, mympi->neighid[1], tag[0], mympi->topocomm, &(mympi->pair_r_reqs[ipair][istage][1]));
      MPI_Recv_init(rbuff_y1, siz_r_y1, MPI_FLOAT, mympi->neighid[2], tag[3], mympi->topocomm, &(mympi->pair_r_reqs[ipair][istage][2]));
      MPI_Recv_init(rbuff_y2, siz_r_y2, MPI_FLOAT, mympi->neighid[3], tag[2], mympi->topocomm, &(mympi->pair_r_reqs[ipair][istage][3]));
    }
  }

  return;
}

void
blk_macdrp_pack_mesg(float *restrict w_cur,float *restrict sbuff,
                 int num_of_vars, gdinfo_t *gdinfo,
                 fd_op_t *fdx_op, fd_op_t *fdy_op)
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

  // to x1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni1; i<ni1+fdx_op->right_len; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // to x2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<=nj2; j++)
      {
        size_t iptr_j = iptr_k + j * siz_line;

        for (int i=ni2- fdx_op->left_len +1; i<=ni2; i++)
        {
          sbuff[iptr_b] = w_cur[iptr_j + i];
          iptr_b++;
        }
      }
    }
  }

  // to y1
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj1; j<nj1+ fdy_op->right_len; j++)
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

  // to y2
  for (int ivar=0; ivar<num_of_vars; ivar++)
  {
    size_t iptr_var = ivar * siz_volume;

    for (int k=nk1; k<=nk2; k++)
    {
      size_t iptr_k = iptr_var + k * siz_slice;

      for (int j=nj2- fdy_op->left_len +1; j<=nj2; j++)
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

  return;
}

void
blk_macdrp_unpack_mesg(float *restrict rbuff,float *restrict w_cur,
                 int num_of_vars, gdinfo_t *gdinfo,
                 fd_op_t *fdx_op, fd_op_t *fdy_op,
                 size_t siz_x1, size_t siz_x2, size_t siz_y1, size_t siz_y2,
                 int *neighid)
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

  // from x1
  float *restrict rbuff_x1 = rbuff;
  if (neighid[0] != MPI_PROC_NULL) {
    iptr_b = 0;
    for (int ivar=0; ivar<num_of_vars; ivar++)
    {
      size_t iptr_var = ivar * siz_volume;

      for (int k=nk1; k<=nk2; k++)
      {
        size_t iptr_k = iptr_var + k * siz_slice;

        for (int j=nj1; j<=nj2; j++)
        {
          size_t iptr_j = iptr_k + j * siz_line;

          for (int i=ni1- fdx_op->left_len; i<ni1; i++)
          {
            w_cur[iptr_j + i] = rbuff_x1[iptr_b];
            iptr_b++;
          }
        }
      }
    }
  }

  // from x2
  float *restrict rbuff_x2 = rbuff_x1 + siz_x1;
  if (neighid[1] != MPI_PROC_NULL) {
    iptr_b = 0;
    for (int ivar=0; ivar<num_of_vars; ivar++)
    {
      size_t iptr_var = ivar * siz_volume;

      for (int k=nk1; k<=nk2; k++)
      {
        size_t iptr_k = iptr_var + k * siz_slice;

        for (int j=nj1; j<=nj2; j++)
        {
          size_t iptr_j = iptr_k + j * siz_line;

          for (int i=ni2+1; i<=ni2+ fdx_op->right_len; i++)
          {
            w_cur[iptr_j + i] = rbuff_x2[iptr_b];
            iptr_b++;
          }
        }
      }
    }
  }

  // from y1
  float *restrict rbuff_y1 = rbuff_x2 + siz_x2;
  if (neighid[2] != MPI_PROC_NULL) {
    iptr_b = 0;
    for (int ivar=0; ivar<num_of_vars; ivar++)
    {
      size_t iptr_var = ivar * siz_volume;

      for (int k=nk1; k<=nk2; k++)
      {
        size_t iptr_k = iptr_var + k * siz_slice;

        for (int j=nj1- fdy_op->left_len; j<nj1; j++)
        {
          size_t iptr_j = iptr_k + j * siz_line;

          for (int i=ni1; i<=ni2; i++)
          {
            w_cur[iptr_j + i] = rbuff_y1[iptr_b];
            iptr_b++;
          }
        }
      }
    }
  }

  // from y2
  float *restrict rbuff_y2 = rbuff_y1 + siz_y1;
  if (neighid[3] != MPI_PROC_NULL) {
    iptr_b = 0;
    for (int ivar=0; ivar<num_of_vars; ivar++)
    {
      size_t iptr_var = ivar * siz_volume;

      for (int k=nk1; k<=nk2; k++)
      {
        size_t iptr_k = iptr_var + k * siz_slice;

        for (int j=nj2+1; j<=nj2+ fdy_op->right_len; j++)
        {
          size_t iptr_j = iptr_k + j * siz_line;

          for (int i=ni1; i<=ni2; i++)
          {
            w_cur[iptr_j + i] = rbuff_y2[iptr_b];
            iptr_b++;
          }
        }
      }
    } // ivar
  }

  return;
}

/*********************************************************************
 * mpi message for el stg scheme
 *********************************************************************/

void
blk_stg_el1st_mesg_init(mympi_t *mympi,
                int ni,
                int nj,
                int nk,
                int fdx_nghosts,
                int fdy_nghosts)
{
  // esti size of mpi mesg 
  int siz_x_per_lay = (nj * nk );
  int siz_y_per_lay = (ni * nk );

  // vel, small than stress
  //int siz_buff = siz_x_per_lay * (2 * fdx_max_half_len - 1) // per var, either side less 1
  //               * CONST_NDIM // ncmp of vel
  //             + siz_y_per_lay * (2 * fdy_max_half_len - 1) // per var, either side less 1
  //               * CONST_NDIM ; // ncmp of vel

  // stress, needs more than vel
  int siz_buff = siz_x_per_lay * (2 * fdx_nghosts - 1) // per var, either side less 1
                 * (CONST_NDIM_2 - 1) // ncmp of stress, one shear is no need
               + siz_y_per_lay * (2 * fdy_nghosts - 1) // per var, either side less 1
                 * (CONST_NDIM_2 - 1) ;

  mympi->sbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc sbuff");
  mympi->rbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc rbuff");
  for (int iptr=0; iptr < siz_buff; iptr++) {
    mympi->rbuff[iptr] = 0.0;
  }
  mympi->siz_sbuff = siz_buff;
  mympi->siz_rbuff = siz_buff;

  // for vel

  int tag_vel[4] = { 11, 12, 21, 22 };

  // send
  int siz_to_x1 = siz_x_per_lay * ( fdx_nghosts  - 1                // Vx
                                   +fdx_nghosts * (CONST_NDIM-1) ); // Vy,Vz

  int siz_to_x2 = siz_x_per_lay * ( fdx_nghosts                         // Vx
                                  +(fdx_nghosts-1) * (CONST_NDIM-1) ); // Vy,Vz

  int siz_to_y1 = siz_y_per_lay * ( fdy_nghosts  - 1                // Vy
                                   +fdy_nghosts * (CONST_NDIM-1) ); // Vx,Vz

  int siz_to_y2 = siz_y_per_lay * ( fdy_nghosts                         // Vy
                                  +(fdy_nghosts-1) * (CONST_NDIM-1) ); // Vx,Vz

  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + siz_to_x1;
  float *sbuff_y1 = sbuff_x2 + siz_to_x2;
  float *sbuff_y2 = sbuff_y1 + siz_to_y1;

  MPI_Send_init(sbuff_x1, siz_to_x1, MPI_FLOAT, mympi->neighid[0], tag_vel[0], mympi->topocomm, &(mympi->s_reqs_vel[0]));
  MPI_Send_init(sbuff_x2, siz_to_x2, MPI_FLOAT, mympi->neighid[1], tag_vel[1], mympi->topocomm, &(mympi->s_reqs_vel[1]));
  MPI_Send_init(sbuff_y1, siz_to_y1, MPI_FLOAT, mympi->neighid[2], tag_vel[2], mympi->topocomm, &(mympi->s_reqs_vel[2]));
  MPI_Send_init(sbuff_y2, siz_to_y2, MPI_FLOAT, mympi->neighid[3], tag_vel[3], mympi->topocomm, &(mympi->s_reqs_vel[3]));

  // recv
  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + siz_to_x2;
  float *rbuff_y1 = rbuff_x2 + siz_to_x1;
  float *rbuff_y2 = rbuff_y1 + siz_to_y2;

  // recv
  MPI_Recv_init(rbuff_x1, siz_to_x2, MPI_FLOAT, mympi->neighid[0], tag_vel[1], mympi->topocomm, &(mympi->r_reqs_vel[0]));
  MPI_Recv_init(rbuff_x2, siz_to_x1, MPI_FLOAT, mympi->neighid[1], tag_vel[0], mympi->topocomm, &(mympi->r_reqs_vel[1]));
  MPI_Recv_init(rbuff_y1, siz_to_y2, MPI_FLOAT, mympi->neighid[2], tag_vel[3], mympi->topocomm, &(mympi->r_reqs_vel[2]));
  MPI_Recv_init(rbuff_y2, siz_to_y1, MPI_FLOAT, mympi->neighid[3], tag_vel[2], mympi->topocomm, &(mympi->r_reqs_vel[3]));

  // for stress

  int tag_stress[4] = { 31, 32, 41, 42 };

  // send
  siz_to_x1 = siz_x_per_lay * ( (fdx_nghosts  ) * (CONST_NDIM  )    // Tii
                               +(fdx_nghosts-1) * (CONST_NDIM-1) ); // Tij

  siz_to_x2 = siz_x_per_lay * ( (fdx_nghosts-1) * (CONST_NDIM  )    // Tii
                               +(fdx_nghosts  ) * (CONST_NDIM-1) ); // Tij

  siz_to_y1 = siz_y_per_lay * ( (fdy_nghosts  ) * (CONST_NDIM  )    // Tii
                               +(fdy_nghosts-1) * (CONST_NDIM-1) ); // Tij

  siz_to_y2 = siz_y_per_lay * ( (fdy_nghosts-1) * (CONST_NDIM  )    // Tii
                               +(fdy_nghosts  ) * (CONST_NDIM-1) ); // Tij

  sbuff_x1 = mympi->sbuff;
  sbuff_x2 = sbuff_x1 + siz_to_x1;
  sbuff_y1 = sbuff_x2 + siz_to_x2;
  sbuff_y2 = sbuff_y1 + siz_to_y1;

  MPI_Send_init(sbuff_x1, siz_to_x1, MPI_FLOAT, mympi->neighid[0], tag_stress[0], mympi->topocomm, &(mympi->s_reqs_stress[0]));
  MPI_Send_init(sbuff_x2, siz_to_x2, MPI_FLOAT, mympi->neighid[1], tag_stress[1], mympi->topocomm, &(mympi->s_reqs_stress[1]));
  MPI_Send_init(sbuff_y1, siz_to_y1, MPI_FLOAT, mympi->neighid[2], tag_stress[2], mympi->topocomm, &(mympi->s_reqs_stress[2]));
  MPI_Send_init(sbuff_y2, siz_to_y2, MPI_FLOAT, mympi->neighid[3], tag_stress[3], mympi->topocomm, &(mympi->s_reqs_stress[3]));

  // recev
  rbuff_x1 = mympi->rbuff;
  rbuff_x2 = rbuff_x1 + siz_to_x2;
  rbuff_y1 = rbuff_x2 + siz_to_x1;
  rbuff_y2 = rbuff_y1 + siz_to_y2;

  MPI_Recv_init(rbuff_x1, siz_to_x2, MPI_FLOAT, mympi->neighid[0], tag_stress[1], mympi->topocomm, &(mympi->r_reqs_stress[0]));
  MPI_Recv_init(rbuff_x2, siz_to_x1, MPI_FLOAT, mympi->neighid[1], tag_stress[0], mympi->topocomm, &(mympi->r_reqs_stress[1]));
  MPI_Recv_init(rbuff_y1, siz_to_y2, MPI_FLOAT, mympi->neighid[2], tag_stress[3], mympi->topocomm, &(mympi->r_reqs_stress[2]));
  MPI_Recv_init(rbuff_y2, siz_to_y1, MPI_FLOAT, mympi->neighid[3], tag_stress[2], mympi->topocomm, &(mympi->r_reqs_stress[3]));

  return;
}

void
blk_stg_el1st_pack_mesg_stress(fdstg_t *fd, gdinfo_t *gdinfo, wav_t *wav,
      float *restrict sbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  // local pointer
  float *restrict Txx   = wav->v5d + wav->Txx_pos;
  float *restrict Tyy   = wav->v5d + wav->Tyy_pos;
  float *restrict Tzz   = wav->v5d + wav->Tzz_pos;
  float *restrict Txz   = wav->v5d + wav->Txz_pos;
  float *restrict Tyz   = wav->v5d + wav->Tyz_pos;
  float *restrict Txy   = wav->v5d + wav->Txy_pos;

  size_t iptr_b = 0;

  //
  // to x 
  //

  // get fd op
  int fdx_max_half_len = fd->fdx_max_half_len;

  // Txx to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<ni1+fdx_max_half_len; i++)
      {
        sbuff[iptr_b] = Txx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tyy to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<ni1+fdx_max_half_len; i++)
      {
        sbuff[iptr_b] = Tyy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tzz to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<ni1+fdx_max_half_len; i++)
      {
        sbuff[iptr_b] = Tzz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txz to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni1; i<ni1+fdx_max_half_len-1; i++)
      {
        sbuff[iptr_b] = Txz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txy to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni1; i<ni1+fdx_max_half_len-1; i++)
      {
        sbuff[iptr_b] = Txy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txx to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni2-fdx_max_half_len+2; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tyy to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni2-fdx_max_half_len+2; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tyy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tzz to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni2-fdx_max_half_len+2; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tzz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txz to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2-fdx_max_half_len+1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txy to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2-fdx_max_half_len+1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  //
  // to y
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // Txx to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j < nj1 + fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tyy to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j < nj1 + fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tyy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tzz to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j < nj1 + fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tzz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tyz to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj1; j < nj1 + fdy_max_half_len -1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tyz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txy to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj1; j < nj1 + fdy_max_half_len -1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txx to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // at j, 1 less
    for (int j=nj2-fdy_max_half_len+2; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tyy to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // at j, 1 less
    for (int j=nj2-fdy_max_half_len+2; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tyy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tzz to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // at j, 1 less
    for (int j=nj2-fdy_max_half_len+2; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tzz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Tyz to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2-fdy_max_half_len+1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Tyz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Txy to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2-fdy_max_half_len+1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Txy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  return;
}

void
blk_stg_el1st_unpack_mesg_stress(fdstg_t *fd,mympi_t *mympi, gdinfo_t *gdinfo, wav_t *wav,
    float *restrict rbuff, size_t siz_rbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;

  // local pointer
  float *restrict Txx   = wav->v5d + wav->Txx_pos;
  float *restrict Tyy   = wav->v5d + wav->Tyy_pos;
  float *restrict Tzz   = wav->v5d + wav->Tzz_pos;
  float *restrict Txz   = wav->v5d + wav->Txz_pos;
  float *restrict Tyz   = wav->v5d + wav->Tyz_pos;
  float *restrict Txy   = wav->v5d + wav->Txy_pos;

  size_t iptr_b = 0;

  //
  // from x 
  //

  // get fd op
  int fdx_max_half_len = fd->fdx_max_half_len;

  // Txx from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni1-fdx_max_half_len+1; i<ni1; i++)
      {
        Txx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tyy from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni1-fdx_max_half_len+1; i<ni1; i++)
      {
        Tyy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tzz from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni1-fdx_max_half_len+1; i<ni1; i++)
      {
        Tzz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txz from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1-fdx_max_half_len; i<ni1; i++)
      {
        Txz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txy from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1-fdx_max_half_len; i<ni1; i++)
      {
        Txy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txx from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
      {
        Txx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tyy from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
      {
        Tyy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tzz from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
      {
        Tzz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txz from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni2+1; i<=ni2+fdx_max_half_len-1; i++)
      {
        Txz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txy from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni2+1; i<=ni2+fdx_max_half_len-1; i++)
      {
        Txy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  //
  // from y 
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // Txx from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj1-fdy_max_half_len+1; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Txx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tyy from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj1-fdy_max_half_len+1; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Tyy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tzz from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj1-fdy_max_half_len+1; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Tzz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tyz from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1-fdy_max_half_len; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Tyz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txy from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1-fdy_max_half_len; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Txy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txx from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Txx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tyy from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Tyy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tzz from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Tzz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Tyz from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj2+1; j<=nj2+fdy_max_half_len-1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Tyz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Txy from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj2+1; j<=nj2+fdy_max_half_len-1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Txy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  if (iptr_b > siz_rbuff) {
    fprintf(stderr, "unpack vel error: iptr_b=%i > siz=%i\n", iptr_b, siz_rbuff);
    exit(-1);
  }

  // reset mpi boundary
  // x1
  if (mympi->neighid[0] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni1-fdx_max_half_len; i<ni1; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Txx[iptr] = 0.0;
          Tyy[iptr] = 0.0;
          Tzz[iptr] = 0.0;
          Tyz[iptr] = 0.0;
          Txz[iptr] = 0.0;
          Txy[iptr] = 0.0;
        }
      }
    }
  }
  // x2
  if (mympi->neighid[1] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Txx[iptr] = 0.0;
          Tyy[iptr] = 0.0;
          Tzz[iptr] = 0.0;
          Tyz[iptr] = 0.0;
          Txz[iptr] = 0.0;
          Txy[iptr] = 0.0;
        }
        // Txy, Txz at i+1/2
        //size_t iptr = ni2 + j * siz_line + k * siz_slice;
        //Txz[iptr] = 0.0;
        //Txy[iptr] = 0.0;
      }
    }
  }
  // y1
  if (mympi->neighid[2] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1-fdy_max_half_len; j<nj1; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Txx[iptr] = 0.0;
          Tyy[iptr] = 0.0;
          Tzz[iptr] = 0.0;
          Tyz[iptr] = 0.0;
          Txz[iptr] = 0.0;
          Txy[iptr] = 0.0;
        }
      }
    }
  }
  // y2
  if (mympi->neighid[3] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Txx[iptr] = 0.0;
          Tyy[iptr] = 0.0;
          Tzz[iptr] = 0.0;
          Tyz[iptr] = 0.0;
          Txz[iptr] = 0.0;
          Txy[iptr] = 0.0;
        } // i
      } // j

      //// Txy, Tyz
      //j = nj2;
      //for (int i=ni1; i<=ni2; i++)
      //{
      //  size_t iptr = i + j * siz_line + k * siz_slice;
      //  Tyz[iptr] = 0.0;
      //  Txy[iptr] = 0.0;
      //} // i
    } // k
  }

  return;
}

void
blk_stg_el1st_pack_mesg_vel(fdstg_t *fd, 
            gdinfo_t *gdinfo, wav_t *wav, float *restrict sbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  size_t iptr_b = 0; 
  //
  // to x 
  //

  // get fd op
  //fd_op_t *this_fd_op = fd->lay_fdx_op + fd->num_of_fdx_op - 1;
  int fdx_max_half_len = fd->fdx_max_half_len;

  // Vx to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vx at i+1/2, one layer already at x neigh
      for (int i=ni1; i<ni1+fdx_max_half_len-1; i++)
      {
        sbuff[iptr_b] = Vx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vy to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<ni1+fdx_max_half_len; i++)
      {
        sbuff[iptr_b] = Vy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vz to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<ni1+fdx_max_half_len; i++)
      {
        sbuff[iptr_b] = Vz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vx to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2-fdx_max_half_len+1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vy to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vy at i, half_len -1 is required
      for (int i=ni2-fdx_max_half_len+2; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vz to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vz at i, half_len -1 is required
      for (int i=ni2-fdx_max_half_len+2; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  //
  // to y
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // Vx to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j < nj1 + fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vy to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vy at j+1/2, 1 less is needed
    for (int j=nj1; j < nj1 + fdy_max_half_len -1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vz to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j < nj1 + fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vx to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vx at j, 1 less
    for (int j=nj2-fdy_max_half_len+2; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vy to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2-fdy_max_half_len+1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vz to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vz at j, 1 less
    for (int j=nj2-fdy_max_half_len+2; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vz[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  return;
}

void
blk_stg_el1st_unpack_mesg_vel(fdstg_t *fd,mympi_t *mympi, gdinfo_t *gdinfo, wav_t *wav,
      float *restrict rbuff, size_t siz_rbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;

  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  size_t iptr_b = 0;

  //
  // from x 
  //

  // get fd op
  int fdx_max_half_len = fd->fdx_max_half_len;

  // Vx from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1-fdx_max_half_len; i<ni1; i++)
      {
        Vx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vy from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vy 1 less
      for (int i=ni1-fdx_max_half_len+1; i<ni1; i++)
      {
        Vy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vz from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1-fdx_max_half_len+1; i<ni1; i++)
      {
        Vz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vx from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vx 1 less
      for (int i=ni2+1; i<=ni2+fdx_max_half_len-1; i++)
      {
        Vx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vy from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
      {
        Vy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vz from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
      {
        Vz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  //
  // from y 
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // Vx from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vx 1 less
    for (int j=nj1-fdy_max_half_len+1; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vy from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1-fdy_max_half_len; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vz from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vz 1 less
    for (int j=nj1-fdy_max_half_len+1; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vx from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vy from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vy 1 less
    for (int j=nj2+1; j<=nj2+fdy_max_half_len-1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vz from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vz[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  if (iptr_b > siz_rbuff) {
    fprintf(stderr, "unpack vel error: iptr_b=%i > siz=%i\n", iptr_b, siz_rbuff);
    exit(-1);
  }

  // reset mpi boundary
  // x1
  if (mympi->neighid[0] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni1-fdx_max_half_len; i<ni1; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vx[iptr] = 0.0;
          Vy[iptr] = 0.0;
          Vz[iptr] = 0.0;
        }
      }
    }
  }
  // x2
  if (mympi->neighid[1] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vx[iptr] = 0.0;
          Vy[iptr] = 0.0;
          Vz[iptr] = 0.0;
        }
        // Vx at i+1/2
        //size_t iptr = ni2 + j * siz_line + k * siz_slice;
        //Vx[iptr] = 0.0;
      }
    }
  }
  // y1
  if (mympi->neighid[2] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1-fdy_max_half_len; j<nj1; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vx[iptr] = 0.0;
          Vy[iptr] = 0.0;
          Vz[iptr] = 0.0;
        }
      }
    }
  }
  // y2
  if (mympi->neighid[3] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vx[iptr] = 0.0;
          Vy[iptr] = 0.0;
          Vz[iptr] = 0.0;
        } // i
      } // j

      //// Vy
      //j = nj2;
      //for (int i=ni1; i<=ni2; i++)
      //{
      //  size_t iptr = i + j * siz_line + k * siz_slice;
      //  Vy[iptr] = 0.0;
      //} // i
    } // k
  }

  return;
}

/*********************************************************************
 * mpi message for ac stg scheme
 *********************************************************************/

void
blk_stg_ac_mesg_init(mympi_t *mympi,
                int ni,
                int nj,
                int nk,
                int fdx_nghosts,
                int fdy_nghosts)
{
  // esti size of mpi mesg 
  int siz_x_per_lay = (nj * nk );
  int siz_y_per_lay = (ni * nk );

  // only one Vi and one P required each side
  int siz_buff = siz_x_per_lay * (2 * fdx_nghosts - 1) // one var, either side less 1
               + siz_y_per_lay * (2 * fdy_nghosts - 1); 

  mympi->sbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc sbuff");
  mympi->rbuff = (float *) fdlib_mem_malloc_1d(siz_buff * sizeof(MPI_FLOAT),"alloc rbuff");
  for (int iptr=0; iptr < siz_buff; iptr++) {
    mympi->rbuff[iptr] = 0.0;
  }
  mympi->siz_sbuff = siz_buff;
  mympi->siz_rbuff = siz_buff;

  // for vel

  int tag_vel[4] = { 11, 12, 21, 22 };

  // send
  int siz_to_x1 = siz_x_per_lay * ( fdx_nghosts  - 1 );                // Vx

  int siz_to_x2 = siz_x_per_lay * ( fdx_nghosts );                         // Vx

  int siz_to_y1 = siz_y_per_lay * ( fdy_nghosts  - 1 );                // Vy

  int siz_to_y2 = siz_y_per_lay * ( fdy_nghosts );                         // Vy

  float *sbuff_x1 = mympi->sbuff;
  float *sbuff_x2 = sbuff_x1 + siz_to_x1;
  float *sbuff_y1 = sbuff_x2 + siz_to_x2;
  float *sbuff_y2 = sbuff_y1 + siz_to_y1;

  MPI_Send_init(sbuff_x1, siz_to_x1, MPI_FLOAT, mympi->neighid[0], tag_vel[0], mympi->topocomm, &(mympi->s_reqs_vel[0]));
  MPI_Send_init(sbuff_x2, siz_to_x2, MPI_FLOAT, mympi->neighid[1], tag_vel[1], mympi->topocomm, &(mympi->s_reqs_vel[1]));
  MPI_Send_init(sbuff_y1, siz_to_y1, MPI_FLOAT, mympi->neighid[2], tag_vel[2], mympi->topocomm, &(mympi->s_reqs_vel[2]));
  MPI_Send_init(sbuff_y2, siz_to_y2, MPI_FLOAT, mympi->neighid[3], tag_vel[3], mympi->topocomm, &(mympi->s_reqs_vel[3]));

  // recv
  float *rbuff_x1 = mympi->rbuff;
  float *rbuff_x2 = rbuff_x1 + siz_to_x2;
  float *rbuff_y1 = rbuff_x2 + siz_to_x1;
  float *rbuff_y2 = rbuff_y1 + siz_to_y2;

  // recv
  MPI_Recv_init(rbuff_x1, siz_to_x2, MPI_FLOAT, mympi->neighid[0], tag_vel[1], mympi->topocomm, &(mympi->r_reqs_vel[0]));
  MPI_Recv_init(rbuff_x2, siz_to_x1, MPI_FLOAT, mympi->neighid[1], tag_vel[0], mympi->topocomm, &(mympi->r_reqs_vel[1]));
  MPI_Recv_init(rbuff_y1, siz_to_y2, MPI_FLOAT, mympi->neighid[2], tag_vel[3], mympi->topocomm, &(mympi->r_reqs_vel[2]));
  MPI_Recv_init(rbuff_y2, siz_to_y1, MPI_FLOAT, mympi->neighid[3], tag_vel[2], mympi->topocomm, &(mympi->r_reqs_vel[3]));

  // for stress

  int tag_stress[4] = { 31, 32, 41, 42 };

  // send
  siz_to_x1 = siz_x_per_lay * (fdx_nghosts  );

  siz_to_x2 = siz_x_per_lay * (fdx_nghosts-1); 

  siz_to_y1 = siz_y_per_lay * (fdy_nghosts  ); 

  siz_to_y2 = siz_y_per_lay * ( fdy_nghosts-1);

  sbuff_x1 = mympi->sbuff;
  sbuff_x2 = sbuff_x1 + siz_to_x1;
  sbuff_y1 = sbuff_x2 + siz_to_x2;
  sbuff_y2 = sbuff_y1 + siz_to_y1;

  MPI_Send_init(sbuff_x1, siz_to_x1, MPI_FLOAT, mympi->neighid[0], tag_stress[0], mympi->topocomm, &(mympi->s_reqs_stress[0]));
  MPI_Send_init(sbuff_x2, siz_to_x2, MPI_FLOAT, mympi->neighid[1], tag_stress[1], mympi->topocomm, &(mympi->s_reqs_stress[1]));
  MPI_Send_init(sbuff_y1, siz_to_y1, MPI_FLOAT, mympi->neighid[2], tag_stress[2], mympi->topocomm, &(mympi->s_reqs_stress[2]));
  MPI_Send_init(sbuff_y2, siz_to_y2, MPI_FLOAT, mympi->neighid[3], tag_stress[3], mympi->topocomm, &(mympi->s_reqs_stress[3]));

  // recev
  rbuff_x1 = mympi->rbuff;
  rbuff_x2 = rbuff_x1 + siz_to_x2;
  rbuff_y1 = rbuff_x2 + siz_to_x1;
  rbuff_y2 = rbuff_y1 + siz_to_y2;

  MPI_Recv_init(rbuff_x1, siz_to_x2, MPI_FLOAT, mympi->neighid[0], tag_stress[1], mympi->topocomm, &(mympi->r_reqs_stress[0]));
  MPI_Recv_init(rbuff_x2, siz_to_x1, MPI_FLOAT, mympi->neighid[1], tag_stress[0], mympi->topocomm, &(mympi->r_reqs_stress[1]));
  MPI_Recv_init(rbuff_y1, siz_to_y2, MPI_FLOAT, mympi->neighid[2], tag_stress[3], mympi->topocomm, &(mympi->r_reqs_stress[2]));
  MPI_Recv_init(rbuff_y2, siz_to_y1, MPI_FLOAT, mympi->neighid[3], tag_stress[2], mympi->topocomm, &(mympi->r_reqs_stress[3]));

  return;
}

void
blk_stg_ac1st_pack_mesg_pressure(fdstg_t *fd, gdinfo_t *gdinfo, wav_t *wav,
      float *restrict sbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  // local pointer
  float *restrict P   = wav->v5d + wav->Txx_pos;

  size_t iptr_b = 0;

  //
  // to x 
  //

  // get fd op
  int fdx_max_half_len = fd->fdx_max_half_len;

  // P to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<ni1+fdx_max_half_len; i++)
      {
        sbuff[iptr_b] = P[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // P to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni2-fdx_max_half_len+2; i<=ni2; i++)
      {
        sbuff[iptr_b] = P[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  //
  // to y
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // P to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j < nj1 + fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = P[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // P to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // at j, 1 less
    for (int j=nj2-fdy_max_half_len+2; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = P[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  return;
}

void
blk_stg_ac1st_unpack_mesg_pressure(fdstg_t *fd,mympi_t *mympi, gdinfo_t *gdinfo, wav_t *wav,
    float *restrict rbuff, size_t siz_rbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;

  // local pointer
  float *restrict P   = wav->v5d + wav->Txx_pos;

  size_t iptr_b = 0;

  //
  // from x 
  //

  // get fd op
  int fdx_max_half_len = fd->fdx_max_half_len;

  // P from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // 1 less
      for (int i=ni1-fdx_max_half_len+1; i<ni1; i++)
      {
        P[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // P from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
      {
        P[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  //
  // from y 
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // P from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // 1 less
    for (int j=nj1-fdy_max_half_len+1; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        P[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // P from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        P[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  if (iptr_b > siz_rbuff) {
    fprintf(stderr, "unpack vel error: iptr_b=%i > siz=%i\n", iptr_b, siz_rbuff);
    exit(-1);
  }

  // reset mpi boundary
  // x1
  if (mympi->neighid[0] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni1-fdx_max_half_len; i<ni1; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          P[iptr] = 0.0;
        }
      }
    }
  }
  // x2
  if (mympi->neighid[1] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          P[iptr] = 0.0;
        }
      }
    }
  }
  // y1
  if (mympi->neighid[2] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1-fdy_max_half_len; j<nj1; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          P[iptr] = 0.0;
        }
      }
    }
  }
  // y2
  if (mympi->neighid[3] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          P[iptr] = 0.0;
        } // i
      } // j
    } // k
  }

  return;
}

void
blk_stg_ac1st_pack_mesg_vel(fdstg_t *fd, 
            gdinfo_t *gdinfo, wav_t *wav, float *restrict sbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_iy;
  size_t siz_slice  = gdinfo->siz_iz;

  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  size_t iptr_b = 0; 

  //
  // to x 
  //

  // get fd op
  //fd_op_t *this_fd_op = fd->lay_fdx_op + fd->num_of_fdx_op - 1;
  int fdx_max_half_len = fd->fdx_max_half_len;

  // Vx to x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vx at i+1/2, one layer already at x neigh
      for (int i=ni1; i<ni1+fdx_max_half_len-1; i++)
      {
        sbuff[iptr_b] = Vx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vx to x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni2-fdx_max_half_len+1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vx[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  //
  // to y
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // Vy to y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vy at j+1/2, 1 less is needed
    for (int j=nj1; j < nj1 + fdy_max_half_len -1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  // Vy to y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj2-fdy_max_half_len+1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        sbuff[iptr_b] = Vy[iptr_j + i];
        iptr_b += 1;
      }
    }
  }

  return;
}

void
blk_stg_ac1st_unpack_mesg_vel(fdstg_t *fd,mympi_t *mympi, gdinfo_t *gdinfo, wav_t *wav,
      float *restrict rbuff, size_t siz_rbuff)
{
  int ni1 = gdinfo->ni1;
  int ni2 = gdinfo->ni2;
  int nj1 = gdinfo->nj1;
  int nj2 = gdinfo->nj2;
  int nk1 = gdinfo->nk1;
  int nk2 = gdinfo->nk2;
  size_t siz_line   = gdinfo->siz_line;
  size_t siz_slice  = gdinfo->siz_slice;

  // local pointer
  float *restrict Vx    = wav->v5d + wav->Vx_pos ;
  float *restrict Vy    = wav->v5d + wav->Vy_pos ;
  float *restrict Vz    = wav->v5d + wav->Vz_pos ;

  size_t iptr_b = 0;

  //
  // from x 
  //

  // get fd op
  int fdx_max_half_len = fd->fdx_max_half_len;

  // Vx from x1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1-fdx_max_half_len; i<ni1; i++)
      {
        Vx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vx from x2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1; j<=nj2; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      // Vx 1 less
      for (int i=ni2+1; i<=ni2+fdx_max_half_len-1; i++)
      {
        Vx[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  //
  // from y 
  //

  // get fd op
  int fdy_max_half_len = fd->fdy_max_half_len;

  // Vy from y1
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    for (int j=nj1-fdy_max_half_len; j<nj1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  // Vy from y2
  for (int k=nk1; k<=nk2; k++)
  {
    size_t iptr_k = k * siz_slice;
    // Vy 1 less
    for (int j=nj2+1; j<=nj2+fdy_max_half_len-1; j++)
    {
      size_t iptr_j = iptr_k + j * siz_line;
      for (int i=ni1; i<=ni2; i++)
      {
        Vy[iptr_j + i] = rbuff[iptr_b];
        iptr_b += 1;
      }
    }
  }

  if (iptr_b > siz_rbuff) {
    fprintf(stderr, "unpack vel error: iptr_b=%i > siz=%i\n", iptr_b, siz_rbuff);
    exit(-1);
  }

  // reset mpi boundary
  // x1
  if (mympi->neighid[0] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni1-fdx_max_half_len; i<ni1; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vx[iptr] = 0.0;
        }
      }
    }
  }
  // x2
  if (mympi->neighid[1] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1; j<=nj2; j++) {
        for (int i=ni2+1; i<=ni2+fdx_max_half_len; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vx[iptr] = 0.0;
        }
      }
    }
  }
  // y1
  if (mympi->neighid[2] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj1-fdy_max_half_len; j<nj1; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vy[iptr] = 0.0;
        }
      }
    }
  }
  // y2
  if (mympi->neighid[3] == MPI_PROC_NULL) {
    for (int k=nk1; k<=nk2; k++) {
      for (int j=nj2+1; j<=nj2+fdy_max_half_len; j++) {
        for (int i=ni1; i<=ni2; i++)
        {
          size_t iptr = i + j * siz_line + k * siz_slice;
          Vy[iptr] = 0.0;
        } // i
      } // j
    } // k
  }

  return;
}

/*********************************************************************
 * estimate dt
 *********************************************************************/

int
blk_dt_esti_curv(gdinfo_t *gdinfo, gd_t *gdcurv, md_t *md,
    float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
    int *dtmaxi, int *dtmaxj, int *dtmaxk)
{
  int ierr = 0;

  float dtmax_local = 1.0e10;
  float Vp;

  float *restrict x3d = gdcurv->x3d;
  float *restrict y3d = gdcurv->y3d;
  float *restrict z3d = gdcurv->z3d;

  for (int k = gdinfo->nk1; k < gdinfo->nk2; k++)
  {
    for (int j = gdinfo->nj1; j < gdinfo->nj2; j++)
    {
      for (int i = gdinfo->ni1; i < gdinfo->ni2; i++)
      {
        size_t iptr = i + j * gdinfo->siz_iy + k * gdinfo->siz_iz;

        if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          Vp = sqrt( (md->lambda[iptr] + 2.0 * md->mu[iptr]) / md->rho[iptr] );
        } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          Vp = sqrt( md->kappa[iptr] / md->rho[iptr] );
        }

        float dtLe = 1.0e20;
        float p0[] = { x3d[iptr], y3d[iptr], z3d[iptr] };

        // min L to 8 adjacent planes
        for (int kk = -1; kk <1; kk++) {
          for (int jj = -1; jj < 1; jj++) {
            for (int ii = -1; ii < 1; ii++) {
              if (ii != 0 && jj !=0 && kk != 0)
              {
                float p1[] = { x3d[iptr-ii], y3d[iptr-ii], z3d[iptr-ii] };
                float p2[] = { x3d[iptr-jj*gdinfo->siz_iy],
                               y3d[iptr-jj*gdinfo->siz_iy],
                               z3d[iptr-jj*gdinfo->siz_iy] };
                float p3[] = { x3d[iptr-kk*gdinfo->siz_iz],
                               y3d[iptr-kk*gdinfo->siz_iz],
                               z3d[iptr-kk*gdinfo->siz_iz] };

                float L = fdlib_math_dist_point2plane(p0, p1, p2, p3);

                if (dtLe > L) dtLe = L;
              }
            }
          }
        }

        // convert to dt
        float dt_point = CFL / Vp * dtLe;

        // if smaller
        if (dt_point < dtmax_local) {
          dtmax_local = dt_point;
          *dtmaxi = i;
          *dtmaxj = j;
          *dtmaxk = k;
          *dtmaxVp = Vp;
          *dtmaxL  = dtLe;
        }

      } // i
    } // i
  } //k

  *dtmax = dtmax_local;

  return ierr;
}

int
blk_dt_esti_cart(gdinfo_t *gdinfo, gd_t *gdcart, md_t *md,
    float CFL, float *dtmax, float *dtmaxVp, float *dtmaxL,
    int *dtmaxi, int *dtmaxj, int *dtmaxk)
{
  int ierr = 0;

  float dtmax_local = 1.0e10;
  float Vp;

  float dx = gdcart->dx;
  float dy = gdcart->dy;
  float dz = gdcart->dz;

  // length to plane
  float p0[] = { 0.0, 0.0, 0.0 };
  float p1[] = {  dx, 0.0, 0.0 };
  float p2[] = { 0.0,  dy, 0.0 };
  float p3[] = { 0.0, 0.0,  dz };

  float dtLe = fdlib_math_dist_point2plane(p0, p1, p2, p3);

  for (int k = gdinfo->nk1; k < gdinfo->nk2; k++)
  {
    for (int j = gdinfo->nj1; j < gdinfo->nj2; j++)
    {
      for (int i = gdinfo->ni1; i < gdinfo->ni2; i++)
      {
        size_t iptr = i + j * gdinfo->siz_iy + k * gdinfo->siz_iz;

        if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          Vp = sqrt( (md->lambda[iptr] + 2.0 * md->mu[iptr]) / md->rho[iptr] );
        } else if (md->medium_type == CONST_MEDIUM_ELASTIC_ISO) {
          Vp = sqrt( md->kappa[iptr] / md->rho[iptr] );
        }

        // convert to dt
        float dt_point = CFL / Vp * dtLe;

        // if smaller
        if (dt_point < dtmax_local) {
          dtmax_local = dt_point;
          *dtmaxi = i;
          *dtmaxj = j;
          *dtmaxk = k;
          *dtmaxVp = Vp;
        }

      } // i
    } // i
  } //k

  *dtmax  = dtmax_local;
  *dtmaxL = dtLe;

  return ierr;
}

float
blk_keep_two_digi(float dt)
{
  char str[40];
  float dt_2;

  sprintf(str, "%4.2e", dt);

  str[3] = '0';

  sscanf(str, "%f", &dt_2);
  
  return dt_2;
}
