#ifndef R_FFPCA_H
#define R_FFPCA_H


#include <R.h>
#include <Rinternals.h>

//#include "pcapack/include/rand_svd.h"
//#include "pcapack/include/fastmap.h"

// Produce a copy of a real SEXP matrix
#define COPYMAT(M, N, X, CPX) (memcpy(REAL(CPX), REAL(X), M*N*sizeof(double)))

#define MIN(a,b) a<b?a:b

#endif
