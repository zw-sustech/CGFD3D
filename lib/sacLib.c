#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "sacLib.h"

int sacExport1C1R(
        char *fnm, 
        float *seismo, 
        float evla,
        float evlo,
        float evel,
        float evdp,
        float stla,
        float stlo,
        float stel,
        float stdp,
        float delta, 
        int npts, 
        char *message)
{
    /*
     *  Function: sacExport1C1R 
     *
     *  Description: write seismogram to file using sac format.
     *
     *  Author: Wenzhong Cao
     *
     *  Arguments:
     *      char    *fnm     : filename of output seimsogram.
     *      float   *seismo  : seismogram data.
     *      float    evla    : event latitude.
     *      float    evlo    : event longitude.
     *      float    evel    : event elevation.
     *      float    evdp    : event depth.
     *      float    stla    : station latitude.
     *      float    stlo    : station longitude.
     *      float    stel    : station elevation.
     *      float    stdp    : station depth.
     *      float    delta   : time interval.
     *      int      npts    : number of samples.
     *      char    *message : error message if need.
     *
     *  Return:
     *      value == 0 if succeed.
     *      value < 0 if failed.
     */

    FILE *fp= NULL;

    int ierr;
    struct SAC_HEAD *sacHeader;
    float o, b, e, cmpaz, cmpinc;
    enum SAC_HEADER_ENUMS iftype;
    int nvhdr, leven;
    char kstnm[8]="-12345"; 
    char kcmpnm[8]="-12345";

    o       = 0;
    b       = 0;
    e       = delta*(npts-1);
    nvhdr   = 6;
    iftype  = ITIME;
    leven   = 1;

    cmpaz  = 0.0;
    cmpinc = 0.0;

    sacHeader = (struct SAC_HEAD*)calloc(1, sizeof(struct SAC_HEAD));
    memset(sacHeader, 0, sizeof(struct SAC_HEAD));

    fp = fopen(fnm, "wb");
    IF_FILE_NULL_RETURN(fp, fnm)

    /* initialize sac header. */
    sacHeaderInit(sacHeader);

    /* set basic information to sac header. */
    ierr = sacHeaderSet(sacHeader, delta, b, e, nvhdr, npts, iftype, leven, 
            evla, evlo, evel, evdp, stla, stlo, stel, stdp, 
            cmpaz, cmpinc, kstnm, kcmpnm, message);
    IF_ERROR_RETURN(ierr)

    /* write sac file. */
    ierr = sacWrite(fp, sacHeader, seismo, message);
    IF_ERROR_RETURN(ierr)

    fclose(fp);
    
    free(sacHeader);

    return 0;
}

void sacHeaderInit(struct SAC_HEAD *sacHeader)
{
    /*
     *  Function: sacHeaderInit
     *
     *  Description: initialize sac header.
     *
     *  Author: Wenzhong Cao
     *
     *  Arguments:
     *      struct SAC_HEAD *sacHeader: sac header.
     *
     *  Return:
     *      NULL.
     */

    sacHeader->delta    = SAC_HEADER_FLOAT_UNDEFINED;
    sacHeader->depmin   = SAC_HEADER_FLOAT_UNDEFINED;
    sacHeader->depmax   = SAC_HEADER_FLOAT_UNDEFINED;
    sacHeader->scale    = SAC_HEADER_FLOAT_UNDEFINED;
    sacHeader->odelta   = SAC_HEADER_FLOAT_UNDEFINED;
    sacHeader->b        = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->e        = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->o        = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->a        = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->fmt      = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->t0       = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->t1       = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->t2       = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->t3       = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->t4       = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->t5       = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->t6       = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->t7       = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->t8       = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->t9       = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->f        = SAC_HEADER_FLOAT_UNDEFINED;    
    sacHeader->resp0    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp1    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp2    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp3    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp4    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp5    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp6    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp7    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp8    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->resp9    = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->stla     = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->stlo     = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->stel     = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->stdp     = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->evla     = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->evlo     = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->evel     = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->evdp     = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->mag      = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->user0    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user1    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user2    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user3    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user4    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user5    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user6    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user7    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user8    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->user9    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->dist     = SAC_HEADER_FLOAT_UNDEFINED;  
    sacHeader->az       = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->baz      = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->gcarc    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->sb       = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->sdelta   = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->depmen   = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->cmpaz    = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->cmpinc   = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->xminimum = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->xmaximum = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->yminimum = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->ymaximum = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused6  = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused7  = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused8  = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused9  = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused10 = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused11 = SAC_HEADER_FLOAT_UNDEFINED; 
    sacHeader->unused12 = SAC_HEADER_FLOAT_UNDEFINED; 

    sacHeader->nzyear   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nzjday   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nzhour   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nzmin    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nzsec    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nzmsec   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nvhdr    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->norid    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nevid    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->npts     = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nsnpts   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nwfid    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nxsize   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->nysize   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused15 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->iftype   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->idep     = SAC_HEADER_INT_UNDEFINED;
    sacHeader->iztype   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused16 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->iinst    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->istreg   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->ievreg   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->ievtyp   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->iqual    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->isynth   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->imagtyp  = SAC_HEADER_INT_UNDEFINED;
    sacHeader->imagsrc  = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused19 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused20 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused21 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused22 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused23 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused24 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused25 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused26 = SAC_HEADER_INT_UNDEFINED;
    sacHeader->leven    = SAC_HEADER_INT_UNDEFINED;
    sacHeader->lpspol   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->lovrok   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->lcalda   = SAC_HEADER_INT_UNDEFINED;
    sacHeader->unused27 = SAC_HEADER_INT_UNDEFINED;

    strncpy(sacHeader->kstnm , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kevnm , SAC_HEADER_CHAR_UNDEFINED, 16);
    strncpy(sacHeader->khole , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->ko    , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->ka    , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt0   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt1   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt2   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt3   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt4   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt5   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt6   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt7   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt8   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kt9   , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kf    , SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kuser0, SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kuser1, SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kuser2, SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kcmpnm, SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->knetwk, SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kdatrd, SAC_HEADER_CHAR_UNDEFINED, 8);
    strncpy(sacHeader->kinst , SAC_HEADER_CHAR_UNDEFINED, 8);

    return;
}

int sacHeaderSet(
        struct SAC_HEAD *sacHeader, 
        float delta, 
        float b, 
        float e, 
        int nvhdr, 
        int npts, 
        int iftype, 
        int leven, 
        float evla, 
        float evlo,
        float evel,
        float evdp,
        float stla, 
        float stlo, 
        float stel, 
        float stdp, 
        float cmpaz, 
        float cmpinc, 
        char *kstnm, 
        char *kcmpnm, 
        char *message)
{
    /*
     *  Function: sacHeaderSet
     *
     *  Description: set basic information to sac header.
     *
     *  Author: Wenzhong Cao
     *
     *  Arguments:
     *      struct SAC_HEAD *sacHeader: sac header.
     *      float           delta     : time interval.
     *      float           b         : initial time.
     *      float           e         : end time.
     *      int             nvhdr     : header version number.
     *      int             npts      : number of samples.
     *      int             iftype    : type of file.
     *      int             leven     : data-evenly-spaced flag.
     *      float           evla      : event latitude.
     *      float           evlo      : event longitude.
     *      float           evel      : event elevation.
     *      float           evdp      : event depth.
     *      float           stla      : station latitude.
     *      float           stlo      : station longitude.
     *      float           stel      : station elevation.
     *      float           stdp      : station depth.
     *      float           cmpaz     : component azimuth.
     *      float           cmpinc    : component inclination.
     *      char            *kstnm    : station name.
     *      char            *kcmpnm   : component name.
     *      char            *message  : error message if needed.
     *
     *  Return:
     *      value == 0 if succeed.
     *      value < 0 if failed.
     */

    if(sacHeader == NULL)
    {
        sprintf(message, "the struct of sacHeader is NULL @ sacHeaderSet\n");
        return -1;
    }

    /* R: required by SAC */
    sacHeader->delta    = delta;
    sacHeader->b        = b;
    sacHeader->e        = e;
    sacHeader->nvhdr    = nvhdr;
    sacHeader->npts     = npts;
    sacHeader->iftype   = iftype;
    sacHeader->leven    = leven;

    /* T: available in SEED header tables */
    sacHeader->stla     = stla;
    sacHeader->stlo     = stlo;
    sacHeader->stel     = stel;
    sacHeader->stdp     = stdp;
    sacHeader->cmpaz    = cmpaz;
    sacHeader->cmpinc   = cmpinc;

    sacHeader->evla     = evla;
    sacHeader->evlo     = evlo;
    sacHeader->evel     = evel;
    sacHeader->evdp     = evdp;

    strncpy(sacHeader->kstnm , kstnm,  8);
    strncpy(sacHeader->kcmpnm, kcmpnm, 8);
    
    return 0;
}

int sacWrite(
        FILE *fp, 
        struct SAC_HEAD *sacHeader, 
        float *data, 
        char *message)
{
    /*
     *  Function: sacWrite
     *
     *  Description: write sac file.
     *
     *  Author: Wenzhong Cao
     *
     *  Arguments:
     *      FILE            *fp       : FILE pointer.
     *      struct SAC_HEAD *sacHeader: sac header.
     *      float           *data     : seismogram data.
     *      char            *message  : error message if needed.
     *
     *  Return:
     *      value == 0 if succeed.
     *      value < 0 if failed.
     */

    if (fp == NULL) {
        sprintf(message, "The FILE pointer is NULL !!! @ sacWrite");
        return -1;
    }

    fwrite(sacHeader, sizeof(struct SAC_HEAD), 1, fp);
    fwrite(data, sizeof(float), sacHeader->npts, fp);

    return 0;
}

