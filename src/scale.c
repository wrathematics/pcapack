#include "pcapack.h"

SEXP R_scale(SEXP centerx, SEXP scalex, SEXP x)
{
  R_INIT;
  const int m = nrows(x), n = ncols(x);
  int info = 0;
  SEXP scaled;
  
  newRmat(scaled, m, n, "double");
  memcpy(REAL(scaled), REAL(x), m*n*sizeof(double));
  
  info = pcapack_scale(INT(centerx), INT(scalex), m, n, REAL(scaled));
  chkinfo(info);
  
  R_END;
  return scaled;
}
