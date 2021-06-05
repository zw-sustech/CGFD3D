#ifndef __SACLIB_H__
#define __SACLIB_H__ 1

/*---------------------------------------------------------------------------------
 * This lib is write by Wenzhong CAO, caowz@mail.ustc.edu.cn. 
 *  
 * The functions are used for writing seismogram in sac format.
 *---------------------------------------------------------------------------------*/
/*
 * The SAC_HEAD & SAC_HEADER_ENUMS are copy or modify from sac.
 *
 * Warning: This sac header file is to be used with the routines 
 *          in the utils directory of the sac distribution.  For 
 *          the header file associated with the SacIO (libsacio.a)
 *          library, check the include directory of the sac
 *          distribution.
 *
 *          Much of the information provided here will be incorporated 
 *          into the header file within the include directory, making
 *          this header file obsolete.  This also pertains to header
 *          files in the inc directory of the source-only sac distribution.
 *
 *          - Brian Savage 01 August 2008 < savage _AT_ uri _DOT_ edu >
 *
 */

/*************************************************************************
  Name:     sac.h

  Purpose:  structure for header of a SAC (Seismic Analysis Code)
        data file, and prototype for basic SAC I/O

  Notes:    Key to comment flags describing each field:
    Column 1:
        R   required by SAC
          (blank)   optional
    Column 2:
        A = settable from a priori knowledge
        D = available in data
        F = available in or derivable from SEED fixed data header
        T = available in SEED header tables
          (blank) = not directly available from SEED data, header
            tables, or elsewhere

  Problems: none known

  References:   O'Neill, D. (1987).  IRIS Interim Data Distribution Format
                (SAC ASCII), Version 1.0 (12 November 1987).  Incorporated
        Research Institutions for Seismology, 1616 North Fort Myer
        Drive, Suite 1440, Arlington, Virginia 22209.  11 pp.
        Tull, J. (1987).  SAC User's Manual, Version 10.2, October 7,
        1987.  Lawrence Livermore National Laboratory, L-205,
        Livermore, California 94550.  ??? pp.

  Language: C, hopefully ANSI standard

  Author:   Dennis O'Neill

  Revisions:    07/15/88  Dennis O'Neill  Initial preliminary release 0.9
                11/21/88  Dennis O'Neill  Production release 1.0
        01/27/91  Lorraine Hwang  Header number is now version 6
        07/06/93  Xiaoming Ding   structure name sac -> SAC_HEAD
                                  typedef structure to be SACHEAD
        12/06/96  Lupei Zhu   prototype sacio functions
**************************************************************************/

#define SAC_HEADER_FLOAT_UNDEFINED (-12345.0)
#define SAC_HEADER_INT_UNDEFINED   (-12345)
#define SAC_HEADER_CHAR_UNDEFINED  ("-12345  ")
#define SAC_HEADER_UNDEFINED       ("UNDEFINED")

struct SAC_HEAD
{
    float delta;      /* RF time increment, sec    */
    float depmin;     /*    minimum amplitude      */
    float depmax;     /*    maximum amplitude      */
    float scale;      /*    amplitude scale factor */
    float odelta;     /*    observed time inc      */
    float b;          /* RD initial time - wrt nz* */
    float e;          /* RD end time               */
    float o;          /*    event start            */
    float a;          /*    1st arrival time       */
    float fmt;        /*    internal use           */
    float t0;         /*    user-defined time pick */
    float t1;         /*    user-defined time pick */
    float t2;         /*    user-defined time pick */
    float t3;         /*    user-defined time pick */
    float t4;         /*    user-defined time pick */
    float t5;         /*    user-defined time pick */
    float t6;         /*    user-defined time pick */
    float t7;         /*    user-defined time pick */
    float t8;         /*    user-defined time pick */
    float t9;         /*    user-defined time pick */
    float f;          /*    event end, sec > 0     */
    float resp0;      /*    instrument respnse parm*/
    float resp1;      /*    instrument respnse parm*/
    float resp2;      /*    instrument respnse parm*/
    float resp3;      /*    instrument respnse parm*/
    float resp4;      /*    instrument respnse parm*/
    float resp5;      /*    instrument respnse parm*/
    float resp6;      /*    instrument respnse parm*/
    float resp7;      /*    instrument respnse parm*/
    float resp8;      /*    instrument respnse parm*/
    float resp9;      /*    instrument respnse parm*/
    float stla;       /*  T station latitude       */
    float stlo;       /*  T station longitude      */
    float stel;       /*  T station elevation, m   */
    float stdp;       /*  T station depth, m       */
    float evla;       /*    event latitude         */
    float evlo;       /*    event longitude        */
    float evel;       /*    event elevation        */
    float evdp;       /*    event depth            */
    float mag;        /*    magnitude value        */
    float user0;      /*    available to user      */
    float user1;      /*    available to user      */
    float user2;      /*    available to user      */
    float user3;      /*    available to user      */
    float user4;      /*    available to user      */
    float user5;      /*    available to user      */
    float user6;      /*    available to user      */
    float user7;      /*    available to user      */
    float user8;      /*    available to user      */
    float user9;      /*    available to user      */
    float dist;       /*    stn-event distance, km */
    float az;         /*    event-stn azimuth      */
    float baz;        /*    stn-event azimuth      */
    float gcarc;      /*    stn-event dist, degrees*/
    float sb;         /*    saved b value          */
    float sdelta;     /*    saved delta value      */
    float depmen;     /*    mean value, amplitude  */
    float cmpaz;      /*  T component azimuth      */
    float cmpinc;     /*  T component inclination  */
    float xminimum;   /*    XYZ X minimum value    */
    float xmaximum;   /*    XYZ X maximum value    */
    float yminimum;   /*    XYZ Y minimum value    */
    float ymaximum;   /*    XYZ Y maximum value    */
    float unused6;    /*    reserved for future use*/
    float unused7;    /*    reserved for future use*/
    float unused8;    /*    reserved for future use*/
    float unused9;    /*    reserved for future use*/
    float unused10;   /*    reserved for future use*/
    float unused11;   /*    reserved for future use*/
    float unused12;   /*    reserved for future use*/

    int   nzyear;     /*  F zero time of file, yr  */
    int   nzjday;     /*  F zero time of file, day */
    int   nzhour;     /*  F zero time of file, hr  */
    int   nzmin;      /*  F zero time of file, min */
    int   nzsec;      /*  F zero time of file, sec */
    int   nzmsec;     /*  F zero time of file, msec*/
    int   nvhdr;      /* R  header version number  */
    int   norid;      /*    Origin ID              */
    int   nevid;      /*    Event ID               */
    int   npts;       /* RF number of samples      */
    int   nsnpts;     /*    saved npts             */
    int   nwfid;      /*    Waveform ID            */
    int   nxsize;     /*    XYZ X size             */
    int   nysize;     /*    XYZ Y size             */
    int   unused15;   /*    reserved for future use*/
    int   iftype;     /* RA type of file          */
    int   idep;       /*    type of amplitude      */
    int   iztype;     /*    zero time equivalence  */
    int   unused16;   /*    reserved for future use*/
    int   iinst;      /*    recording instrument   */
    int   istreg;     /*    stn geographic region  */
    int   ievreg;     /*    event geographic region*/
    int   ievtyp;     /*    event type             */
    int   iqual;      /*    quality of data        */
    int   isynth;     /*    synthetic data flag    */
    int   imagtyp;    /*    magnitude type         */
    int   imagsrc;    /*    magnitude source       */
    int   unused19;   /*    reserved for future use*/
    int   unused20;   /*    reserved for future use*/
    int   unused21;   /*    reserved for future use*/
    int   unused22;   /*    reserved for future use*/
    int   unused23;   /*    reserved for future use*/
    int   unused24;   /*    reserved for future use*/
    int   unused25;   /*    reserved for future use*/
    int   unused26;   /*    reserved for future use*/
    int   leven;      /* RA data-evenly-spaced flag*/
    int   lpspol;     /*    station polarity flag  */
    int   lovrok;     /*    overwrite permission   */
    int   lcalda;     /*    calc distance, azimuth */
    int   unused27;   /*    reserved for future use*/

    char  kstnm[8];   /*  F station name           */
    char  kevnm[16];  /*    event name             */
    char  khole[8];   /*    man-made event name    */
    char  ko[8];      /*    event origin time id   */
    char  ka[8];      /*    1st arrival time ident */
    char  kt0[8];     /*    time pick 0 ident      */
    char  kt1[8];     /*    time pick 1 ident      */
    char  kt2[8];     /*    time pick 2 ident      */
    char  kt3[8];     /*    time pick 3 ident      */
    char  kt4[8];     /*    time pick 4 ident      */
    char  kt5[8];     /*    time pick 5 ident      */
    char  kt6[8];     /*    time pick 6 ident      */
    char  kt7[8];     /*    time pick 7 ident      */
    char  kt8[8];     /*    time pick 8 ident      */
    char  kt9[8];     /*    time pick 9 ident      */
    char  kf[8];      /*    end of event ident     */
    char  kuser0[8];  /*    available to user      */
    char  kuser1[8];  /*    available to user      */
    char  kuser2[8];  /*    available to user      */
    char  kcmpnm[8];  /*  F component name         */
    char  knetwk[8];  /*    network name           */
    char  kdatrd[8];  /*    date data read         */
    char  kinst[8];   /*    instrument name        */
};

enum SAC_HEADER_ENUMS {
  /* enumerated header values */
  IREAL    = 0,   /* Undocumented                */
  ITIME    = 1,   /* Time series file            */
  IRLIM    = 2,   /* Spectral file-real/imag     */
  IAMPH    = 3,   /* Spectral file-ampl/phase    */
  IXY      = 4,   /* General x vs y file         */
  IUNKN    = 5,   /* Unknown                     */
  IDISP    = 6,   /* Displacement (NM)           */
  IVEL     = 7,   /* Velocity (NM/SEC)           */
  IACC     = 8,   /* Acceleration (NM/SEC/SEC)   */
  IB       = 9,   /* Begin time                  */
  IDAY     = 10,  /* GMT day                     */
  IO       = 11,  /* Event origin time           */
  IA       = 12,  /* First arrival time          */
  IT0      = 13,  /* User defined time pick 0    */
  IT1      = 14,  /* User defined time pick 1    */
  IT2      = 15,  /* User defined time pick 2    */
  IT3      = 16,  /* User defined time pick 3    */
  IT4      = 17,  /* User defined time pick 4    */
  IT5      = 18,  /* User defined time pick 5    */
  IT6      = 19,  /* User defined time pick 6    */
  IT7      = 20,  /* User defined time pick 7    */
  IT8      = 21,  /* User defined time pick 8    */
  IT9      = 22,  /* User defined time pick 9    */
  IRADNV   = 23,  /* Radial (NTS)                */
  ITANNV   = 24,  /* Tangential (NTS)            */
  IRADEV   = 25,  /* Radial (EVENT)              */
  ITANEV   = 26,  /* Tangential (EVENT)          */
  INORTH   = 27,  /* North positive              */
  IEAST    = 28,  /* East positive               */
  IHORZA   = 29,  /* Horizontal (ARB)            */
  IDOWN    = 30,  /* Down positive               */
  IUP      = 31,  /* Up positive                 */
  ILLLBB   = 32,  /* LLL broadband               */
  IWWSN1   = 33,  /* WWSN 15-100                 */
  IWWSN2   = 34,  /* WWSN 30-100                 */
  IHGLP    = 35,  /* High-gain long-period       */
  ISRO     = 36,  /* SRO                         */
  INUCL    = 37,  /* Nuclear event               */
  IPREN    = 38,  /* Nuclear pre-shot event      */
  IPOSTN   = 39,  /* Nuclear post-shot event     */
  IQUAKE   = 40,  /* Earthquake                  */
  IPREQ    = 41,  /* Foreshock                   */
  IPOSTQ   = 42,  /* Aftershock                  */
  ICHEM    = 43,  /* Chemical explosion          */
  IOTHER   = 44,  /* Other                       */
  IGOOD    = 45,  /* Good                        */
  IGLCH    = 46,  /* Gliches                     */
  IDROP    = 47,  /* Dropouts                    */
  ILOWSN   = 48,  /* Low signal to noise ratio   */
  IRLDTA   = 49,  /* Real data                   */
  IVOLTS   = 50,  /* Velocity (volts)            */
  IXYZ     = 51,  /* General XYZ (3-D) file      */
  /* These 18 added to describe magnitude type and source maf 970205 */
  IMB      = 52,  /* Bodywave Magnitude          */
  IMS      = 53,  /* Surface Magnitude           */
  IML      = 54,  /* Local Magnitude             */
  IMW      = 55,  /* Moment Magnitude            */
  IMD      = 56,  /* Duration Magnitude          */
  IMX      = 57,  /* User Defined Magnitude      */
  INEIC    = 58,  /* INEIC                       */
  IPDEQ    = 59,  /* IPDEQ                       */
  IPDEW    = 60,  /* IPDEW                       */
  IPDE     = 61,  /* IPDE                        */
  IISC     = 62,  /* IISC                        */
  IREB     = 63,  /* IREB                        */
  IUSGS    = 64,  /* IUSGS                       */
  IBRK     = 65,  /* IBRK                        */
  ICALTECH = 66,  /* ICALTECH                    */
  ILLNL    = 67,  /* ILLNL                       */
  IEVLOC   = 68,  /* IEVLOC                      */
  IJSOP    = 69,  /* IJSOP                       */
  IUSER    = 70,  /* IUSER                       */
  IUNKNOWN = 71,  /* IUNKNOWN                    */
  /* These  17 added for ievtyp. maf 970325 */
  IQB      = 72,  /* Quarry or mine blast confirmed by quarry */
  IQB1     = 73,  /* Quarry or mine blast with designed shot information-ripple fired */
  IQB2     = 74,  /* Quarry or mine blast with observed shot information-ripple fired */
  IQBX     = 75,  /* Quarry or mine blast - single shot */
  IQMT     = 76,  /* Quarry or mining-induced events: tremors and rockbursts */
  IEQ      = 77,  /* Earthquake                  */
  IEQ1     = 78,  /* Earthquakes in a swarm or aftershock sequence */
  IEQ2     = 79,  /* Felt earthquake             */
  IME      = 80,  /* Marine explosion            */
  IEX      = 81,  /* Other explosion             */
  INU      = 82,  /* Nuclear explosion           */
  INC      = 83,  /* Nuclear cavity collapse     */
  IO_      = 84,  /* Other source of known origin */
  IL       = 85,  /* Local event of unknown origin */
  IR       = 86,  /* Regional event of unknown origin */
  IT       = 87,  /* Teleseismic event of unknown origin */
  IU       = 88,  /* Undetermined or conflicting information  */
  /* These 9 added for ievtype to keep up with database. maf 000530 */
  IEQ3     = 89,  /* Damaging Earthquake         */
  IEQ0     = 90,  /* Probable earthquake         */
  IEX0     = 91,  /* Probable explosion          */
  IQC      = 92,  /* Mine collapse               */
  IQB0     = 93,  /* Probable Mine Blast         */
  IGEY     = 94,  /* Geyser                      */
  ILIT     = 95,  /* Light                       */
  IMET     = 96,  /* Meteroic event              */
  IODOR    = 97   /* Odors                       */
};


#ifndef __ABOUT_ERROR__
#define __ABOUT_ERROR__ 1
  #define STATUS_ERROR -1
  #define STATUS_FOPEN_ERROR -1
  
  /* if FILE pointer is NULL, return STATUS_FOPEN_ERROR. */
  #define IF_FILE_NULL_RETURN(fp, fnm) if (fp==NULL) {sprintf(message, "%s:%d File pointer is NULL, filename is '%s'", __FILE__, __LINE__, fnm); return STATUS_FOPEN_ERROR;}
  /* if FILE pointer is NULL, set ierr=STATUS_FOPEN_ERROR, then goto EXIT. */
  #define IF_FILE_NULL_EXIT(fp, fnm)   if (fp==NULL) {ierr=STATUS_FOPEN_ERROR; sprintf(message, "%s:%d File pointers is NULL, filename is '%s'", __FILE__, __LINE__, fnm); goto EXIT;}
  /* if ierr is less than 0, return ierr. */
  #define IF_ERROR_RETURN(ierr) if ( ierr < 0 ) {return ierr;}
  /* if var is less than 0, goto EXIT. */
  #define IF_ERROR_EXIT(ierr)   if ( ierr < 0 ) {goto EXIT;}
  /* if keyword not found in file, return error. */
  #define IF_ERROR_KEYWORD_IN_FILE_RETURN(ierr, fnm, varnm) if (ierr < 1) {sprintf(message, "%s:%d You must set the [%s] in file [%s]!!!", __FILE__, __LINE__, varnm, fnm); return STATUS_ERROR;}

  #define RETURN_ERROR_MESSAGE(errormsg, status) {snprintf(message, STRLEN, "%s:%d %s", __FILE__, __LINE__, errormsg); return status;}
#endif

/*---------------------------------------------------------------------------------
  This part is for writing sac format seismogram.
  ---------------------------------------------------------------------------------*/
/* write seismogram to file using sac format.*/
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
        float dt, 
        int nt, 
        char *message);

/* initialize sac header. */
void sacHeaderInit(struct SAC_HEAD *sacheader);

/* set basic information to sac header. */
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
        char *message);

/* write sac file. */
int sacWrite(
        FILE *fp, 
        struct SAC_HEAD *sacheader, 
        float *data, 
        char *message);

#endif
