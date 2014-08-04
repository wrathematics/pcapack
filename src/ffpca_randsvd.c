#include "ffpca.h"

#define MIN(a,b) (a<b?a:b)
#define PTCT(x) (PROTECT(x);ptct++)

// R wrapper for drandsvd_
SEXP R_randsvd(SEXP METHOD, SEXP M, SEXP N, SEXP A, SEXP K, SEXP Q, SEXP JOBU, SEXP JOBVT)
{
  const char jobu = RCHAR(JOBU)[0], jobvt = RCHAR(JOBVT)[0];
  
  const int m = INTEGER(M)[0], n = INTEGER(N)[0], k = INTEGER(K)[0];
  const int rc_min = MIN(m, n);
  
  int info = 0;
  int numret = 0;
  
  double *a;
  
  SEXP S, U, VT;
  SEXP RET, RET_NAMES;
  
  
  // Allocate
  PROTECT(S = allocVector(REALSXP, k));
  if (jobu == 'V')
  {
    numret++;
    PROTECT(U = allocMatrix(REALSXP, m, k));
  }
  else
    PROTECT(U = allocVector(REALSXP, 1));
  
  if (jobvt == 'V')
  {
    numret++;
    PROTECT(VT = allocMatrix(REALSXP, n, k));
  }
  else
    PROTECT(U = allocVector(REALSXP, 1));
  
  a = R_alloc(m*n, sizeof(double));
  memcpy(a, REAL(A), m*n*sizeof(double));
  
  
  // Compute
  drandsvd_(RCHAR(METHOD), RCHAR(JOBU), RCHAR(JOBVT), &m, &n, a, &k, INTEGER(Q), REAL(S), REAL(U), REAL(VT), &info);
  
  // Return wrangling --- herein I do a very silly thing to try to avoid redundant copying
  PROTECT(RET = allocVector(VECSXP, numret));
  PROTECT(RET_NAMES = allocVector(STRSXP, numret));
  setAttrib(RET, R_NamesSymbol, RET_NAMES);
  
  SET_VECTOR_ELT(RET, 0, S);
  SET_STRING_ELT(RET_NAMES, 0, mkChar("d")); 
  
  if (jobu == 'V')
  {
    SET_VECTOR_ELT(RET, 1, U);
    SET_STRING_ELT(RET_NAMES, 1, mkChar("u")); 
    
    if (jobvt == 'V')
    {
      SET_VECTOR_ELT(RET, 2, VT);
      SET_STRING_ELT(RET_NAMES, 2, mkChar("vt")); 
    }
  }
  else
  {
    if (jobvt == 'V')
    {
      SET_VECTOR_ELT(RET, 1, VT);
      SET_STRING_ELT(RET_NAMES, 1, mkChar("vt")); 
    }
  }
  
  
  UNPROTECT(5);
  return RET;
} 


