#include "pcapack.h"


SEXP R_pcapack_randsvd(SEXP retu, SEXP NV, SEXP X)
{
  R_INIT;
  int m = nrows(X), n = ncols(X);
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
  
  
  info = pcapack_randsvd(INT(retu), INT(retvt), INT(k), const int niter, const int method, const int m, const int n, double *restrict x, double *restrict s, double *restrict u, double *restrict vt);
  info = pcapack_randsvd(INT(NU), INT(NV), INT(METHOD), m, n, *x, int k, int niter, double *s, bool retu, double *u, bool retvt, double *vt)
  
  if (retu && retvt)
  {
    RET_NAMES = make_list_names(3, "d", "u", "vt");
    RET = make_list(RET_NAMES, 3, S, U, VT);
  }
  
  R_END;
  return RET;
} 




