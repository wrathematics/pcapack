#ifndef __PCAPACK_MATRIX_H__
#define __PCAPACK_MATRIX_H__


#include <stdlib.h>

typedef struct dblvec_t
{
  int length;
  double *restrict data;
} dblvec_t;

typedef struct intvec_t
{
  int length;
  int *restrict data;
} intvec_t;

typedef struct matrix_t
{
  int nrows;
  int ncols;
  double *restrict data;
} matrix_t;


#define newmat(m,n,x) x=malloc(m*n*sizeof(*x))
#define new0mat(m,n,x) x=calloc(m*n, sizeof(*x))
#define copymat(m, n, x, y) memcpy(y, x, m*n*sizeof(*x))


#endif
