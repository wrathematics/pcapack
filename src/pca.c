#include "pcapack.h"

#define setDimNames(X, Y, Z) \
  newRlist(X, 2); \
  SET_VECTOR_ELT(X, 0, RNULL); \
  SET_VECTOR_ELT(X, 1, Y); \
  setAttrib(Z, R_DimNamesSymbol, X);


SEXP R_pcapack_prcomp_svd(SEXP X, SEXP CENTER, SEXP SCALE, SEXP RETX)
{
  R_INIT;
  int info = 0;
  const int m = nrows(X), n = ncols(X);
  const bool retx = (bool) INT(RETX);
  const bool centerx = (bool) INT(CENTER);
  const bool scalex = (bool) INT(SCALE);
  const int k = m<n?m:n;
  double *x;
  SEXP RET, RET_NAMES, SDEV, TROT, XRET;
  SEXP pcnames, dimnames, xnames;
  
  
  newRvec(SDEV, k, "double");
  newRmat(TROT, k, n, "double");
  
  if (retx)
  {
    newRmat(XRET, m, n, "double"); // FIXME need to pass in mxn, return mxk
    memcpy(DBLP(XRET), DBLP(X), m*n*sizeof(double));
    x = DBLP(XRET);
  }
  else
    x = DBLP(X);
  
  pcapack_prcomp_svd(centerx, scalex, retx, m, n, x, DBLP(SDEV), DBLP(TROT));
  
/*  if (info != 0)*/
/*    error(_("info=%d from Lapack routine '%s'"), info, "dgesdd");*/
  
  pcnames = make_pca_default_colnames(n);
  setDimNames(dimnames, pcnames, TROT);
  
  if (retx)
  {
    xnames = make_pca_default_colnames(n);
    setDimNames(dimnames, xnames, XRET);
    
    RET_NAMES = make_list_names(3, "sdev", "rotation", "x");
    RET = make_list(RET_NAMES, 3, SDEV, TROT, XRET);
  }
  else
  {
    RET_NAMES = make_list_names(2, "sdev", "rotation");
    RET = make_list(RET_NAMES, 2, SDEV, TROT);
  }
  
  
  R_END;
  return RET;
}



SEXP R_pcapack_prcomp_eigcov(SEXP X, SEXP RETX)
{
  R_INIT;
  const int m = nrows(X), n = ncols(X);
  bool retx = (bool) INT(RETX);
  double *x;
  int info = 0;
  
  const int k = m<n?m:n;
  
  SEXP RET, RET_NAMES, SDEV, TROT, XRET;
  SEXP pcnames, dimnames, xnames;
  
  
  newRvec(SDEV, k, "double");
  newRmat(TROT, k, n, "double");
  
  if (retx)
  {
    newRmat(XRET, m, n, "double"); // FIXME need to pass in mxn, return mxk
    memcpy(DBLP(XRET), DBLP(X), m*n*sizeof(double));
    x = DBLP(XRET);
  }
  else
    x = DBLP(X);
  
  info = pcapack_prcomp_eigcov(&retx, m, n, x, DBLP(SDEV), DBLP(TROT));
  
  pcnames = make_pca_default_colnames(n);
  setDimNames(dimnames, pcnames, TROT);
  

  if (retx)
  {
    xnames = make_pca_default_colnames(n);
    setDimNames(dimnames, xnames, XRET);
    
    RET_NAMES = make_list_names(3, "sdev", "rotation", "x");
    RET = make_list(RET_NAMES, 3, SDEV, TROT, XRET);
  }
  else
  {
    RET_NAMES = make_list_names(2, "sdev", "rotation");
    RET = make_list(RET_NAMES, 2, SDEV, TROT);
  }
  

  R_END;
  return RET;
}

