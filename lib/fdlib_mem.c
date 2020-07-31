#include "utility_vector.h"

/*
 * todo:
 *  1 add functions for int, double, char
 *  2 complete header file
 *  3 add required include
 */

/******************************************************************************
 * allocate float array with initial value
 *****************************************************************************/

float *fdlib_mem_calloc_1d_float(size_t n, float v0, char *msg)
{
    int ierr;
    float *buff;
    size_t i;

    ierr = 0;

    if (n <= 0 ) {
        fprintf(stderr, "Error: the length of buff %i is negative!\n", n);
        return NULL;
    }

    // allocate
    buff = (float*) malloc( n * sizeof(float));

    if ( buff == NULL ) {
        fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        ierr = -1;
        exit(ierr);
    }

    // init value
    for (i = 0; i < n; i++ ) {
        buff[i] = v0;
    }

    return buff;
}

/******************************************************************************
  allocate float array without initial value
 *****************************************************************************/

float *fdlib_mem_malloc_1d_float(size_t n, char *msg)
{
    int ierr;
    float *buff;

    ierr = 0;

    if (n <= 0 ) {
        fprintf(stderr, "Error: the length of buff %i is negative!\n", n);
        return NULL;
    }

    // allocate
    buff = (float*) malloc( n * sizeof(float));

    if ( buff == NULL ) {
        fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        ierr = -1;
        exit(ierr);
    }

    return buff;
}

/******************************************************************************
 * free 1d array no matter data type
 *****************************************************************************/

int fdlib_mem_free_1d(void *p)
{
   if (p==NULL) {
       fprintf(stderr, "Pointer is null!\n");
       return -1;
   }

    free(p);
    p = NULL;

    return 0;
}

/******************************************************************************
 * allocate 2d char array
 *****************************************************************************/

char **fdlib_mem_malloc_2d_char(size_t n, int str_len, char *msg)
{
    int ierr;
    float **buff;
    size_t i;

    ierr = 0;

    if (n <= 0 ) {
        fprintf(stderr, "Error: the length of char array %i is negative!\n", n);
        return NULL;
    }

    // allocate
    buff = (char**) malloc( n * sizeof(char *));

    if ( buff == NULL ) {
        fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        ierr = -1;
        exit(ierr);
    }

    // allocate 2
    for (i=0; i<n; i++) {
        buff[i] = (char *)malloc(str_len*sizeof(char));
        // ini buff with empty
        strcpy(buff[i],
    }

    return buff;
}

