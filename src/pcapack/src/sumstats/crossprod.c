// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include "sumstats.h"
#include "../misc.h"


// make symmetric via copying from one triangle to the other.
int pcapack_symmetrize(const int triang, const int m, const int n, double *x)
{
  int i, j;
  const int k = MIN(m, n);
  
  if (m == 0 || n == 0) return 0;
  if (triang != UPPER && triang != LOWER) return -1;
  
  // NOTE keep these serial, as all the cache misses inherent to these make the performance tank
  
  // Copy upper ONTO lower
  if (likely(triang == UPPER))
  {
    for (j=0; j<k; j++)
    {
      for (i=0; i<j; i++)
        x[j + m*i] = x[i + m*j];
    }
  }
  // Copy lower ONTO upper
  else if (unlikely(triang == LOWER))
  {
    for (j=0; j<k; j++)
    {
      for (i=j+1; i<k; i++)
        x[j + m*i] = x[i + m*j];
    }
  }
  
  return 0;
}



// t(x) * x
int pcapack_crossprod(int m, int n, double *x, double alpha, double *c)
{
  int info = 0;
  
  dsyrk_(&(char){'u'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
  info = pcapack_symmetrize(UPPER, n, n, c);
  
  return info;
}


// x * t(x)
int pcapack_tcrossprod(int m, int n, double *x, double alpha, double *c)
{
  int info = 0;
  
  dsyrk_(&(char){'u'}, &(char){'n'}, &m, &n, &alpha, x, &m, &(double){0.0}, c, &m);
  info = pcapack_symmetrize(UPPER, m, m, c);
  
  return info;
}



int pcapack_inverse(int n, double *x)
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
  dgetri_(&n, x, &n, ipiv, work);
  
  
  free(work);
  cleanup:
  free(ipiv);
  
  return info;
}

