#include <stdlib.h>
#include <stdio.h>
#include "fdlib_mem.h"

/*
 * todo:
 *  1 add functions for int, double, char
 *  2 complete header file
 *  3 add required include
 */

/******************************************************************************
 * 1D array
 *****************************************************************************/

//-------------------------------------------------------------------------------
// allocate 1D array without(malloc) initial value, no explicit datatype required
//-------------------------------------------------------------------------------

void *
fdlib_mem_malloc_1d(size_t siz_byte, char *msg)
{
    if (siz_byte <= 0 ) {
        fprintf(stderr, "Error: siz_byte=%i is zero or negative (%s)(!\n", siz_byte, msg);
        fflush(stderr);
        return NULL;
    }

    // allocate
    void *buff = (void *) malloc( siz_byte );
    if ( buff == NULL ) {
        fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    return buff;
}

//-------------------------------------------------------------------------------
// allocate 1D array with(calloc) initial value, explicit datatype is required
//-------------------------------------------------------------------------------

int *
fdlib_mem_calloc_1d_int(size_t n, int v0, char *msg)
{
  if (n <= 0 ) {
    fprintf(stderr, "Error: size=%i is zero or negative (%s)!\n", n, msg);
        fflush(stderr);
    return NULL;
  }
  
  int *buff = (int *) malloc( n * sizeof(int));
  if ( buff == NULL ) {
    fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        fflush(stderr);
    exit(-1);
  }
  
  for (size_t i = 0; i < n; i++ ) buff[i] = v0;
  
  return buff;
}

size_t *
fdlib_mem_calloc_1d_sizet(size_t n, size_t v0, char *msg)
{
  if (n <= 0 ) {
    fprintf(stderr, "Error: size=%i is zero or negative (%s)!\n", n, msg);
        fflush(stderr);
    return NULL;
  }
  
  size_t *buff = (size_t *) malloc( n * sizeof(size_t));
  if ( buff == NULL ) {
    fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        fflush(stderr);
    exit(-1);
  }
  
  for (size_t i = 0; i < n; i++ ) buff[i] = v0;
  
  return buff;
}

float *
fdlib_mem_calloc_1d_float(size_t n, float v0, char *msg)
{
  if (n <= 0 ) {
    fprintf(stderr, "Error: size=%i is zero or negative (%s)!\n", n, msg);
        fflush(stderr);
    return NULL;
  }

  float *var = (float *) malloc( n * sizeof(float));
  if ( var == NULL ) {
    fprintf(stderr, "Error: can't malloc enough mem (%s)!\n", msg);
        fflush(stderr);
    exit(-1);
  }

  for (size_t i = 0; i < n; i++ ) var[i] = v0;

  return var;
}

//-------------------------------------------------------------------------------
// free 1d array no matter data type
//-------------------------------------------------------------------------------

void 
fdlib_mem_free_1d(void *p)
{
  if (p==NULL) {
    fprintf(stderr, "Pointer is null!\n");
        fflush(stderr);
  } else {
    free(p);
  }
}

/******************************************************************************
 * 2nd-level pointers, have to with exlicit datatype?
 *****************************************************************************/

// int type
int **
fdlib_mem_malloc_2l_int(size_t n1, size_t n2, char *msg)
{
  int **buff;
  
  if (n1 <= 0 || n2<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i negative (%s)!\n", n1, n2, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (int**) malloc( n1 * sizeof(int*));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (int*) malloc( n2 * sizeof(int));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }
  }
  
  return buff;
}

int **
fdlib_mem_calloc_2l_int(size_t n1, size_t n2, int v0, char *msg)
{

  int **buff = fdlib_mem_malloc_2l_int(n1, n2, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      buff[i][j] = v0;

  return buff;
}

void
fdlib_mem_free_2l_int(int **p, size_t n1, char *msg)
{
  for (size_t i=0; i < n1; i++ ) {
    free(p[i]);
  }
  free(p);
}

// size_t
size_t **
fdlib_mem_malloc_2l_sizet(size_t n1, size_t n2, char *msg)
{
  size_t **buff;
  
  if (n1 <= 0 || n2<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i negative (%s)!\n", n1, n2, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (size_t**) malloc( n1 * sizeof(size_t*));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (size_t*) malloc( n2 * sizeof(size_t));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }
  }
  
  return buff;
}

size_t **
fdlib_mem_calloc_2l_sizet(size_t n1, size_t n2, size_t v0, char *msg)
{

  size_t **buff = fdlib_mem_malloc_2l_sizet(n1, n2, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      buff[i][j] = v0;

  return buff;
}

void
fdlib_mem_free_2l_sizet(size_t **p, size_t n1, char *msg)
{
  for (size_t i=0; i < n1; i++ ) {
    free(p[i]);
  }
  free(p);
}

// float
float **
fdlib_mem_malloc_2l_float(size_t n1, size_t n2, char *msg)
{
  float **buff;
  
  if (n1 <= 0 || n2<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i negative (%s)!\n", n1, n2, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (float**) malloc( n1 * sizeof(float*));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (float*) malloc( n2 * sizeof(float));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }
  }
  
  return buff;
}

float **
fdlib_mem_calloc_2l_float(size_t n1, size_t n2, float v0, char *msg)
{

  float **buff = fdlib_mem_malloc_2l_float(n1, n2, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      buff[i][j] = v0;

  return buff;
}

void
fdlib_mem_free_2l_float(float **p, size_t n1, char *msg)
{
  for (size_t i=0; i < n1; i++ ) {
    free(p[i]);
  }
  free(p);
}

// char
char **
fdlib_mem_malloc_2l_char(size_t n1, size_t n2, char *msg)
{
  char **buff;
  
  if (n1 <= 0 || n2<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i negative (%s)!\n", n1, n2, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (char**) malloc( n1 * sizeof(char*));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (char*) malloc( n2 * sizeof(char));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }
  }
  
  return buff;
}

void
fdlib_mem_free_2l_char(char **p, size_t n1, char *msg)
{
  for (size_t i=0; i < n1; i++ ) {
    free(p[i]);
  }
  free(p);
}

/******************************************************************************
 * 3rd-level pointers, have to with exlicit datatype?
 *****************************************************************************/

// int
int ***
fdlib_mem_malloc_3l_int(size_t n1, size_t n2, size_t n3, char *msg)
{
  int ***buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i negative (%s)!\n", n1, n2, n3, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (int***) malloc( n1 * sizeof(int**));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (int**) malloc( n2 * sizeof(int*));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (int*) malloc( n3 * sizeof(int));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }
    }
  }
  
  return buff;
}

int ***
fdlib_mem_calloc_3l_int(size_t n1, size_t n2, size_t n3, int v0, char *msg)
{

  int ***buff = fdlib_mem_malloc_3l_int(n1, n2, n3, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      for (size_t k=0; k < n3; k++)
       buff[i][j][k] = v0;

  return buff;
}

void
fdlib_mem_free_3l_int(int ***p, size_t n1, size_t n2, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) free(p[i][j]);

    free(p[i]);
  }
  free(p);
}

// size_t
size_t ***
fdlib_mem_malloc_3l_sizet(size_t n1, size_t n2, size_t n3, char *msg)
{
  size_t ***buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i negative (%s)!\n", n1, n2, n3, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (size_t***) malloc( n1 * sizeof(size_t**));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (size_t**) malloc( n2 * sizeof(size_t*));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (size_t*) malloc( n3 * sizeof(size_t));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }
    }
  }
  
  return buff;
}

size_t ***
fdlib_mem_calloc_3l_sizet(size_t n1, size_t n2, size_t n3, size_t v0, char *msg)
{

  size_t ***buff = fdlib_mem_malloc_3l_sizet(n1, n2, n3, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      for (size_t k=0; k < n3; k++)
       buff[i][j][k] = v0;

  return buff;
}

void
fdlib_mem_free_3l_sizet(size_t ***p, size_t n1, size_t n2, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) free(p[i][j]);

    free(p[i]);
  }
  free(p);
}

// float
float ***
fdlib_mem_malloc_3l_float(size_t n1, size_t n2, size_t n3, char *msg)
{
  float ***buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i negative (%s)!\n", n1, n2, n3, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (float***) malloc( n1 * sizeof(float**));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (float**) malloc( n2 * sizeof(float*));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (float*) malloc( n3 * sizeof(float));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }
    }
  }
  
  return buff;
}

float ***
fdlib_mem_calloc_3l_float(size_t n1, size_t n2, size_t n3, float v0, char *msg)
{

  float ***buff = fdlib_mem_malloc_3l_float(n1, n2, n3, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      for (size_t k=0; k < n3; k++)
       buff[i][j][k] = v0;

  return buff;
}

void
fdlib_mem_free_3l_float(float ***p, size_t n1, size_t n2, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) free(p[i][j]);

    free(p[i]);
  }
  free(p);
}

// char
char ***
fdlib_mem_malloc_3l_char(size_t n1, size_t n2, size_t n3, char *msg)
{
  char ***buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i negative (%s)!\n", n1, n2, n3, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (char***) malloc( n1 * sizeof(char**));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (char**) malloc( n2 * sizeof(char*));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (char*) malloc( n3 * sizeof(char));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }
    }
  }
  
  return buff;
}

void
fdlib_mem_free_3l_char(char ***p, size_t n1, size_t n2, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) free(p[i][j]);

    free(p[i]);
  }
  free(p);
}

/******************************************************************************
 * 4th-level pointers, have to with exlicit datatype?
 *****************************************************************************/

// int
int ****
fdlib_mem_malloc_4l_int(size_t n1, size_t n2, size_t n3, size_t n4, char *msg)
{
  int ****buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0 || n4<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i n4=%i negative (%s)!\n", n1, n2, n3, n4, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (int****) malloc( n1 * sizeof(int***));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (int***) malloc( n2 * sizeof(int**));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (int**) malloc( n3 * sizeof(int*));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }

      // alloc 4th-level
      for (size_t k=0; k<n3; k++)
      {
        buff[i][j][k] = (int*) malloc( n4 * sizeof(int));
        if ( buff[i][j][k] == NULL ) {
            fprintf(stderr, "Error: allocate 4 fails (%s)!\n", msg);
        fflush(stderr);
            exit(-1);
        }
      }
    }
  }
  
  return buff;
}

int ****
fdlib_mem_calloc_4l_int(size_t n1, size_t n2, size_t n3, size_t n4, int v0, char *msg)
{

  int ****buff = fdlib_mem_malloc_4l_int(n1, n2, n3, n4, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      for (size_t k=0; k < n3; k++)
        for (size_t m=0; m < n4; m++)
       buff[i][j][k][m] = v0;

  return buff;
}

void
fdlib_mem_free_4l_int(int ****p, size_t n1, size_t n2, size_t n3, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) 
    {
      for (size_t k=0; k < n3; k++) 
      {
        free(p[i][j][k]);
      }

      free(p[i][j]);
    }

    free(p[i]);
  }
  free(p);
}

// size_t
size_t ****
fdlib_mem_malloc_4l_sizet(size_t n1, size_t n2, size_t n3, size_t n4, char *msg)
{
  size_t ****buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0 || n4<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i n4=%i negative (%s)!\n", n1, n2, n3, n4, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (size_t****) malloc( n1 * sizeof(size_t***));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (size_t***) malloc( n2 * sizeof(size_t**));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (size_t**) malloc( n3 * sizeof(size_t*));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }

      // alloc 4th-level
      for (size_t k=0; k<n3; k++)
      {
        buff[i][j][k] = (size_t*) malloc( n4 * sizeof(size_t));
        if ( buff[i][j][k] == NULL ) {
            fprintf(stderr, "Error: allocate 4 fails (%s)!\n", msg);
        fflush(stderr);
            exit(-1);
        }
      }
    }
  }
  
  return buff;
}

size_t ****
fdlib_mem_calloc_4l_sizet(size_t n1, size_t n2, size_t n3, size_t n4, size_t v0, char *msg)
{

  size_t ****buff = fdlib_mem_malloc_4l_sizet(n1, n2, n3, n4, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      for (size_t k=0; k < n3; k++)
        for (size_t m=0; m < n4; m++)
       buff[i][j][k][m] = v0;

  return buff;
}

void
fdlib_mem_free_4l_sizet(size_t ****p, size_t n1, size_t n2, size_t n3, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) 
    {
      for (size_t k=0; k < n3; k++) 
      {
        free(p[i][j][k]);
      }

      free(p[i][j]);
    }

    free(p[i]);
  }
  free(p);
}

// float
float ****
fdlib_mem_malloc_4l_float(size_t n1, size_t n2, size_t n3, size_t n4, char *msg)
{
  float ****buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0 || n4<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i n4=%i negative (%s)!\n", n1, n2, n3, n4, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (float****) malloc( n1 * sizeof(float***));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (float***) malloc( n2 * sizeof(float**));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (float**) malloc( n3 * sizeof(float*));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }

      // alloc 4th-level
      for (size_t k=0; k<n3; k++)
      {
        buff[i][j][k] = (float*) malloc( n4 * sizeof(float));
        if ( buff[i][j][k] == NULL ) {
            fprintf(stderr, "Error: allocate 4 fails (%s)!\n", msg);
        fflush(stderr);
            exit(-1);
        }
      }
    }
  }
  
  return buff;
}

float ****
fdlib_mem_calloc_4l_float(size_t n1, size_t n2, size_t n3, size_t n4, float v0, char *msg)
{

  float ****buff = fdlib_mem_malloc_4l_float(n1, n2, n3, n4, msg);

  for (size_t i=0; i < n1; i++ ) 
    for (size_t j=0; j < n2; j++)
      for (size_t k=0; k < n3; k++)
        for (size_t m=0; m < n4; m++)
       buff[i][j][k][m] = v0;

  return buff;
}

void
fdlib_mem_free_4l_float(float ****p, size_t n1, size_t n2, size_t n3, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) 
    {
      for (size_t k=0; k < n3; k++) 
      {
        free(p[i][j][k]);
      }

      free(p[i][j]);
    }

    free(p[i]);
  }
  free(p);
}

// char
char ****
fdlib_mem_malloc_4l_char(size_t n1, size_t n2, size_t n3, size_t n4, char *msg)
{
  char ****buff;
  
  if (n1 <= 0 || n2<=0 || n3<=0 || n4<=0) {
    fprintf(stderr, "Error: n1=%i n2=%i n3=%i n4=%i negative (%s)!\n", n1, n2, n3, n4, msg);
        fflush(stderr);
    return NULL;
  }
  
  // alloc 1st-level
  buff = (char****) malloc( n1 * sizeof(char***));
  if ( buff == NULL ) {
      fprintf(stderr, "Error: allocate 1 fails (%s)!\n", msg);
        fflush(stderr);
      exit(-1);
  }
  
  // alloc 2nd-level
  for (size_t i=0; i<n1; i++)
  {
    buff[i] = (char***) malloc( n2 * sizeof(char**));
    if ( buff[i] == NULL ) {
        fprintf(stderr, "Error: allocate 2 fails (%s)!\n", msg);
        fflush(stderr);
        exit(-1);
    }

    // alloc 3rd-level
    for (size_t j=0; j<n2; j++)
    {
      buff[i][j] = (char**) malloc( n3 * sizeof(char*));
      if ( buff[i][j] == NULL ) {
          fprintf(stderr, "Error: allocate 3 fails (%s)!\n", msg);
        fflush(stderr);
          exit(-1);
      }
      // alloc 4th-level
      for (size_t k=0; k<n3; k++)
      {
        buff[i][j][k] = (char*) malloc( n4 * sizeof(char));
        if ( buff[i][j][k] == NULL ) {
            fprintf(stderr, "Error: allocate 4 fails (%s)!\n", msg);
        fflush(stderr);
            exit(-1);
        }
      }
    }
  }
  
  return buff;
}

void
fdlib_mem_free_4l_char(char ****p, size_t n1, size_t n2, size_t n3, char *msg)
{
  for (size_t i=0; i < n1; i++ )
  {
    for (size_t j=0; j < n2; j++) 
    {
      for (size_t k=0; k < n3; k++) 
      {
        free(p[i][j][k]);
      }

      free(p[i][j]);
    }

    free(p[i]);
  }
  free(p);
}

