#include "pcapack.h"

SEXP R_pca_predict(SEXP centerx, SEXP scalex, SEXP x)
{
  R_INIT;
  const int m = nrows(x), n = ncols(x);
  int info = 0;
  SEXP scaled;
  
  newRmat(scaled, m, n, "double");
  memcpy(REAL(scaled), REAL(x), m*n*sizeof(double));
  
  info = pcapack_scale(INT(centerx), INT(scalex), m, n, REAL(scaled));
  
  R_END;
  return scaled;
}




scale(newdata, object$center, object$scale) %*% object$rotation

