#include "ffpca.h"

// 
SEXP R_pca(SEXP M, SEXP N, SEXP K, SEXP X, SEXP CENTER, SEXP SCALE, SEXP RETROT)
{
  const int m = INTEGER(M)[0], n = INTEGER(N)[0], k = INTEGER(K)[0];
  
  SEXP RET, RET_NAMES, SDEV, TROT, INFO;
  
  /* Protect R objects. */
  PROTECT(RET = allocVector(VECSXP, 3));
  PROTECT(RET_NAMES = allocVector(STRSXP, 3));
  PROTECT(INFO = allocVector(INTSXP, 1));
  PROTECT(SDEV = allocVector(REALSXP, k));
  PROTECT(TROT = allocMatrix(REALSXP, k, n));
  
  SET_VECTOR_ELT(RET, 0, INFO);
  SET_VECTOR_ELT(RET, 1, SDEV);
  SET_VECTOR_ELT(RET, 2, TROT);
  SET_STRING_ELT(RET_NAMES, 0, mkChar("info")); 
  SET_STRING_ELT(RET_NAMES, 1, mkChar("sdev")); 
  SET_STRING_ELT(RET_NAMES, 2, mkChar("trot"));
  setAttrib(RET, R_NamesSymbol, RET_NAMES);
  
  INTEGER(INFO)[0] = 0;
  dpca_(&m, &n, &k, REAL(X), REAL(SDEV), REAL(TROT), CHARPT(RETROT, 0), 
    CHARPT(CENTER, 0), CHARPT(SCALE, 0), INTEGER(INFO));
  
  UNPROTECT(5);
  return RET;
} 


