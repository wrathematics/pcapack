#include "ffpca.h"

// 
SEXP R_DTRAN(SEXP M, SEXP N, SEXP A)
{
  const int m = INTEGER(M)[0], n = INTEGER(N)[0];
  const double alpha = 1.0;
  const double beta = 0.0;
  
  SEXP C;
  PROTECT(C = allocMatrix(REALSXP, n, m));
  
  dtran_(&m, &n, &alpha, REAL(A), &beta, REAL(C));
  
  UNPROTECT(1);
  return C;
} 

