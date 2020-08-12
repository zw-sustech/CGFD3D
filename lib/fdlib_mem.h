#ifndef FDLIB_MEM_H
#define FDLIB_MEM_H

void *
fdlib_mem_malloc_1d(size_t siz_byte, char *msg);

int *
fdlib_mem_calloc_1d_int(size_t n, int v0, char *msg);

size_t *
fdlib_mem_calloc_1d_sizet(size_t n, size_t v0, char *msg);

float *
fdlib_mem_calloc_1d_float(size_t n, float v0, char *msg);

void 
fdlib_mem_free_1d(void *p);

int **
fdlib_mem_malloc_2l_int(size_t n1, size_t n2, char *msg);

int **
fdlib_mem_calloc_2l_int(size_t n1, size_t n2, int v0, char *msg);

void
fdlib_mem_free_2l_int(int **p, size_t n1, char *msg);

size_t **
fdlib_mem_malloc_2l_sizet(size_t n1, size_t n2, char *msg);

size_t **
fdlib_mem_calloc_2l_sizet(size_t n1, size_t n2, size_t v0, char *msg);

void
fdlib_mem_free_2l_sizet(size_t **p, size_t n1, char *msg);

float **
fdlib_mem_malloc_2l_float(size_t n1, size_t n2, char *msg);

float **
fdlib_mem_calloc_2l_float(size_t n1, size_t n2, float v0, char *msg);

void
fdlib_mem_free_2l_float(float **p, size_t n1, char *msg);

char **
fdlib_mem_malloc_2l_char(size_t n1, size_t n2, char *msg);

void
fdlib_mem_free_2l_char(char **p, size_t n1, char *msg);

int ***
fdlib_mem_malloc_3l_int(size_t n1, size_t n2, size_t n3, char *msg);

int ***
fdlib_mem_calloc_3l_int(size_t n1, size_t n2, size_t n3, int v0, char *msg);

void
fdlib_mem_free_3l_int(int ***p, size_t n1, size_t n2, char *msg);

size_t ***
fdlib_mem_malloc_3l_sizet(size_t n1, size_t n2, size_t n3, char *msg);

size_t ***
fdlib_mem_calloc_3l_sizet(size_t n1, size_t n2, size_t n3, size_t v0, char *msg);

void
fdlib_mem_free_3l_sizet(size_t ***p, size_t n1, size_t n2, char *msg);

float ***
fdlib_mem_malloc_3l_float(size_t n1, size_t n2, size_t n3, char *msg);

float ***
fdlib_mem_calloc_3l_float(size_t n1, size_t n2, size_t n3, float v0, char *msg);

void
fdlib_mem_free_3l_float(float ***p, size_t n1, size_t n2, char *msg);

char ***
fdlib_mem_malloc_3l_char(size_t n1, size_t n2, size_t n3, char *msg);

void
fdlib_mem_free_3l_char(char ***p, size_t n1, size_t n2, char *msg);

int ****
fdlib_mem_malloc_4l_int(size_t n1, size_t n2, size_t n3, size_t n4, char *msg);

int ****
fdlib_mem_calloc_4l_int(size_t n1, size_t n2, size_t n3, size_t n4, int v0, char *msg);

void
fdlib_mem_free_4l_int(int ****p, size_t n1, size_t n2, size_t n3, char *msg);

size_t ****
fdlib_mem_malloc_4l_sizet(size_t n1, size_t n2, size_t n3, size_t n4, char *msg);

size_t ****
fdlib_mem_calloc_4l_sizet(size_t n1, size_t n2, size_t n3, size_t n4, size_t v0, char *msg);

void
fdlib_mem_free_4l_sizet(size_t ****p, size_t n1, size_t n2, size_t n3, char *msg);

float ****
fdlib_mem_malloc_4l_float(size_t n1, size_t n2, size_t n3, size_t n4, char *msg);

float ****
fdlib_mem_calloc_4l_float(size_t n1, size_t n2, size_t n3, size_t n4, float v0, char *msg);

void
fdlib_mem_free_4l_float(float ****p, size_t n1, size_t n2, size_t n3, char *msg);

char ****
fdlib_mem_malloc_4l_char(size_t n1, size_t n2, size_t n3, size_t n4, char *msg);

void
fdlib_mem_free_4l_char(char ****p, size_t n1, size_t n2, size_t n3, char *msg);

#endif
