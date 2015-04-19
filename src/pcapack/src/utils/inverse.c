// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include "../lapack.h"


/**
 * @file
 * @brief Matrix inverse.
 *
 * @details
 * Matrix inverse via LU decomposition.  Uses dgetrf and dgetri.
 * 
 * @param n
 * Inputs.  Problem size (rows/cols of x)
 * @param x
 * Input/Output.  The data matrix.
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_inverse(int n, double *restrict x)
{
  int info = 0;
  int *ipiv;
  int lwork;
  double tmp;
  double *work;
  
  
  // Factor x = LU
  ipiv = malloc(n * sizeof(*ipiv));
  dgetrf_(&n, &n, x, &n, ipiv, &info);
  if (info != 0) goto cleanup;
  
  
  // Invert
  lwork = -1;
  dgetri_(&n, x, &n, ipiv, &tmp, &lwork, &info);
  if (info != 0) goto cleanup;
  
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  dgetri_(&n, x, &n, ipiv, work, &lwork, &info);
  
  
  free(work);
  cleanup:
  free(ipiv);
  
  return info;
}

