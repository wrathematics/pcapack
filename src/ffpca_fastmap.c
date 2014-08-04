#include "ffpca.h"

// 
SEXP R_cma(SEXP M, SEXP N, SEXP X, SEXP K)
{
  const int m = INTEGER(M)[0], n = INTEGER(N)[0];
  SEXP CPX;
  PROTECT(CPX = allocMatrix(REALSXP, m, n));
  
  memcpy(REAL(CPX), REAL(X), m*n*sizeof(double));
  
  cma_(&m, &n, REAL(CPX), INTEGER(K));
  
  UNPROTECT(1);
  return CPX;
}


