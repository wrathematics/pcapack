#ifndef R_FFPCA_H
#define R_FFPCA_H


#include <R.h>
#include <Rinternals.h>

#include "pcapack/include/rand_svd.h"
#include "pcapack/include/fastmap.h"

// Produce a copy of a real SEXP matrix
#define COPYMAT(M, N, X, CPX) (memcpy(REAL(CPX), REAL(X), M*N*sizeof(double)))

// Character pointers
#define CHARPT(x,i) ((char*)CHAR(STRING_ELT(x,i)))
#define RCHAR(x) ((char*)CHAR(STRING_ELT(x,0)))


#endif
