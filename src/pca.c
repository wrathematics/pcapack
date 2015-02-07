#include "pcapack.h"

#define setDimNames(X, Y, Z) \
  newRlist(X, 2); \
  SET_VECTOR_ELT(X, 0, RNULL); \
  SET_VECTOR_ELT(X, 1, Y); \
  setAttrib(Z, R_DimNamesSymbol, X);


SEXP R_pca_svd(SEXP M, SEXP N, SEXP X, SEXP CENTER, SEXP SCALE, SEXP RETROT)
{
  R_INIT;
  const int m = INT(M), n = INT(N);
  bool retrot = (bool) INTEGER(RETROT)[0];
  bool center = (bool) INTEGER(CENTER)[0];
  bool scale = (bool) INTEGER(SCALE)[0];
  int info = 0;
  
  const int k = m<n?m:n;
  
  SEXP RET, RET_NAMES, SDEV, TROT;
  SEXP pcnames, dimnames;
  
  
  newRvec(SDEV, k, "double");
  newRmat(TROT, k, n, "double");
  
  prcomp_svd_(&m, &n, DBLP(X), DBLP(SDEV), DBLP(TROT), &retrot, &center, &scale, &info);
  
/*  if (info != 0)*/
/*    error(_("info=%d from Lapack routine '%s'"), info, "dgesdd");*/
  
  pcnames = make_pca_default_colnames(n);
  setDimNames(dimnames, pcnames, TROT);
  
  RET_NAMES = make_list_names(2, "sdev", "rotation");
  RET = make_list(RET_NAMES, 2, SDEV, TROT);
  
  R_END;
  return RET;
}




SEXP R_pca_eigcov(SEXP M, SEXP N, SEXP X, SEXP RETROT)
{
  R_INIT;
  const int m = INT(M), n = INT(N);
  bool retrot = (bool) INTEGER(RETROT)[0];
  int info = 0;
  
  const int k = m<n?m:n;
  
  SEXP RET, RET_NAMES, SDEV, TROT;
  SEXP pcnames, dimnames;
  
  
  newRvec(SDEV, k, "double");
  newRmat(TROT, k, n, "double");
  
  prcomp_eigcov_(&m, &n, DBLP(X), DBLP(SDEV), DBLP(TROT), &retrot, &info);
  
/*  if (info != 0)*/
/*    error(_("info=%d from Lapack routine '%s'"), info, "dgesdd");*/
  
  pcnames = make_pca_default_colnames(n);
  setDimNames(dimnames, pcnames, TROT);
  
  RET_NAMES = make_list_names(2, "sdev", "rotation");
  RET = make_list(RET_NAMES, 2, SDEV, TROT);
  
  R_END;
  return RET;
} 


