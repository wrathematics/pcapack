#include "pcapack.h"


SEXP R_pcapack_svd(SEXP NU, SEXP NV, SEXP M, SEXP N, SEXP X)
{
  R_INIT;
  int m = INT(M), n = INT(N);
  int info = 0;
  const int minmn = MIN(m, n);
  
  const bool retu = INT(NU) ? true : false;
  const bool retvt = INT(NV) ? true : false;
  
  SEXP RET, RET_NAMES;
  SEXP S, U, VT;
  double *u, *vt;
  
  newRvec(S, minmn, "dbl");
  
  if (INT(NU) == 0 && INT(NV) == 0)
  {
    u = NULL;
    vt = NULL;
  }
  else if (INT(NU) == 0 && m >= n)
  {
    u = NULL;
    newRmat(VT, minmn, n, "dbl");
    vt = REAL(VT);
  }
  else if (INT(NV) == 0 && m < n)
  {
    vt = NULL;
    newRmat(U, m, minmn, "dbl");
    u = REAL(U);
  }
  else if (INT(NU) <= minmn && INT(NV) <= minmn)
  {
    newRmat(U, m, minmn, "dbl");
    newRmat(VT, minmn, n, "dbl");
    u = REAL(U);
    vt = REAL(VT);
  }
  else
  {
    newRmat(U, m, m, "dbl");
    newRmat(VT, n, n, "dbl");
    u = REAL(U);
    vt = REAL(VT);
  }
  
  
  info = pcapack_svd(false, INT(NU), INT(NV), m, n, DBLP(X), DBLP(S), u, vt);
  
/*  if (info != 0)*/
/*    error(_("info=%d from Lapack routine '%s'"), info, "dgesdd");*/
  
  if (retu && retvt)
  {
    RET_NAMES = make_list_names(3, "d", "u", "vt");
    RET = make_list(RET_NAMES, 3, S, U, VT);
  }
  else if (retu)
  {
    RET_NAMES = make_list_names(2, "d", "u");
    RET = make_list(RET_NAMES, 2, S, U);
  }
  else if (retvt)
  {
    RET_NAMES = make_list_names(2, "d", "vt");
    RET = make_list(RET_NAMES, 2, S, VT);
  }
  else
  {
    RET_NAMES = make_list_names(1, "d");
    RET = make_list(RET_NAMES, 1, S);
  }
  
  R_END;
  return RET;
} 


