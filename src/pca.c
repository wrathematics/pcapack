#include "ffpca.h"
#include <stdbool.h>
#include <RNACI.h>

// 
SEXP R_pca(SEXP M, SEXP N, SEXP K, SEXP X, SEXP CENTER, SEXP SCALE, SEXP RETROT)
{
  R_INIT;
  const int m = INT(M), n = INT(N), k = INT(K);
  bool retrot = (bool) INTEGER(RETROT)[0];
  bool center = (bool) INTEGER(CENTER)[0];
  bool scale = (bool) INTEGER(SCALE)[0];
  int info = 0;
  
  SEXP RET, RET_NAMES, SDEV, TROT;
  
  
  newRvec(SDEV, k, "double");
  newRmat(TROT, k, n, "double");
  
  prcomp_svd_(&m, &n, &k, DBLP(X), DBLP(SDEV), DBLP(TROT), &retrot, 
    &center, &scale, &info);
  
/*  if (info < 0)*/
/*    warning()*/
  
  RET_NAMES = make_list_names(2, "sdev", "rotation");
  RET = make_list(RET_NAMES, 2, SDEV, TROT);
  
  R_END;
  return RET;
} 


