#ifndef R_PCAPACK_H
#define R_PCAPACK_H


#include "pcapack/src/pcapack.h"
#include <stdbool.h>
#include <RNACI.h>

//#include "pcapack/include/rand_svd.h"
//#include "pcapack/include/fastmap.h"

// Produce a copy of a real SEXP matrix
#define COPYMAT(M, N, X, CPX) (memcpy(REAL(CPX), REAL(X), M*N*sizeof(double)))

#define MIN(m,n) m<n?m:n
#define MAX(m,n) m<n?n:m

#define chkinfo(info) if(info!=0)error("returned non-zero info\n")


// scale.c
SEXP R_scale(SEXP centerx, SEXP scalex, SEXP x);

// svd.c
SEXP R_pcapack_svd(SEXP NU, SEXP NV, SEXP X);

// utils.c
SEXP make_pca_default_colnames(const int n);


#endif
