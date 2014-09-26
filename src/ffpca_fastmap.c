#include "ffpca.h"

void cma_(int *n, int *p, double *x, int *k, int *info);

// 
SEXP R_cma(SEXP X, SEXP K)
{
  int m = nrows(X), n = ncols(X);
  int info = 0;
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  cma_(&m, &n, REAL(CPX), INTEGER(K), &info);
  
  if (info != 0) error("TODO\n");
  
  UNPROTECT(1);
  return CPX;
}


