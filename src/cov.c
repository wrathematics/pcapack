#include "pcapack.h"


SEXP R_cov(SEXP x)
{
  R_INIT;
  double *x_cp;
  const int m = nrows(x), n = ncols(x);
  int info = 0;
  SEXP cov;
  
  
  newRmat(cov, n, n, "double");
  
  x_cp = (double *) R_alloc(m*n, sizeof(double));
  memcpy(x_cp, REAL(x), m*n*sizeof(double));
  
  info = pcapack_cov(COR_PEARSON, m, n, x_cp, REAL(cov));
  chkinfo(info);
  
  R_END;
  return cov;
}



SEXP R_cov_naive(SEXP x)
{
  R_INIT;
  double *x_cp;
  const int m = nrows(x), n = ncols(x);
  int info = 0;
  SEXP cov;
  
  
  newRmat(cov, n, n, "double");
  
  x_cp = (double *) R_alloc(m*n, sizeof(double));
  memcpy(x_cp, REAL(x), m*n*sizeof(double));
  
  info = pcapack_cov_naive(m, n, x_cp, REAL(cov));
  chkinfo(info);
  
  R_END;
  return cov;
}
