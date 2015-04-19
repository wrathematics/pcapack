#include "pcapack.h"


SEXP R_pcapack_cma(SEXP X, SEXP K)
{
  R_INIT;
  const int n = nrows(X), p = ncols(X);
  const int k = INT(K);
  int info;
  SEXP CPX;
  newRmat(CPX, n, p, "dbl");
  
  memcpy(DBLP(CPX), DBLP(X), n*p*sizeof(double));
  
  info = pcapack_cma(n, p, DBLP(CPX), k);
  
  if (info != 0) error("TODO\n");
  
  R_END;
  return CPX;
}

