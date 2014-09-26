#include "ffpca.h"
#include <stdbool.h>
#include <RNACI.h>


void LA_svd_(int *nu, int *nv, int *m, int *n, double *x, double *s, 
  double *u, double *vt, int *info);

SEXP R_pcapack_svd(SEXP NU, SEXP NV, SEXP M, SEXP N, SEXP X)
{
  R_INIT;
  int m = INT(M), n = INT(N);
  int info = 0;
  const int minmn = MIN(m, n);
  
  const bool retu = INT(NU)?true:false;
  const bool retvt = INT(NV)?true:false;
  
  SEXP RET, RET_NAMES;
  SEXP S, U, VT;
  
  newRvec(S, minmn, "dbl");
  newRmat(U, m, minmn, "dbl");
  newRmat(VT, minmn, n, "dbl");
  
  LA_svd_(INTP(NU), INTP(NV), &m, &n, DBLP(X), DBLP(S), DBLP(U), DBLP(VT), &info);
  
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


