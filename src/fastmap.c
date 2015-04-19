#include "pcapack.h"

int pcapack_cma(int n, int p, double *x, int k);

SEXP R_cma(SEXP X, SEXP K)
{
  R_INIT;
  const int n = nrows(X), p = ncols(X);
  const int k = INT(K);
  int info;
  SEXP CPX;
  newRmat(CPX, n, p, "dbl");
  
  memcpy(DBL(CPX), DBL(X), n*p*sizeof(double));
  
  info = pcapack_cma(n, p, DBL(CPX), k;
  
  if (info != 0) error("TODO\n");
  
  R_END;
  return CPX;
}


