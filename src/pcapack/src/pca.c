// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pcapack.h"
#include "lapack.h"
#include "misc.h"


int pcapack_prcomp_svd(bool centerx, bool scalex, bool retx, int m, int n, double *restrict x, double *restrict sdev, double *restrict rotation)
{
  char trans = 'n';
  int info;
  int minmn = MIN(m, n);
  double *x_cp;
  double *u;
  double tmp;
  
  
  if (centerx || scalex || retx)
  {
    x_cp = malloc(m*n * sizeof(*x));
    memcpy(x_cp, x, m*n*sizeof(*x));
  }
  else
    x_cp = x;
  
  u = malloc(m*minmn * sizeof(*x));
  
  
  pcapack_scale(centerx, scalex, m, n, x_cp);
  
  info = pcapack_svd(false, n, n, m, n, x_cp, sdev, u, rotation);
  if (info != 0) goto cleanup;
  
  pcapack_xpose(minmn, n, rotation);
  
  if (retx)
    dgemm_(&trans, &trans, &m, &minmn, &n, &(double){1.0}, x_cp, &m, rotation, &n, &(double){0.0}, x, &m);
  
  tmp = 1. / (MAX(1., sqrt((double) m-1)));
  dscal_(&minmn, &tmp, sdev, &(int){1});
  
  cleanup:
    if (centerx || scalex || retx) free(x_cp);
    free(u);
  
  return info;
}



static inline void sqrt_rev(int n, double *sdev)
{
  int i;
  double tmp;
  
  for (i=0; i<n/2; i++)
  {
    tmp = sqrt(sdev[i]);
    sdev[i] = sqrt(sdev[n-i-1]);
    sdev[n-i-1] = tmp;
  }
  
  if (n % 2 == 1) sdev[n/2] = sqrt(sdev[n/2]);
  
  return;
}

int pcapack_prcomp_eig(bool retx, int m, int n, double *x, double *sdev, double *rotation)
{
  int info = 0;
  int i;
  double tmp;
  double *cov;
  double *x_cp;
  char trans = 'N';
  
  
  if (retx)
  {
    x_cp = malloc(m*n * sizeof(*x));
    memcpy(x_cp, x, m*n*sizeof(*x));
  }
  else
    x_cp = x;
  
  cov = malloc(n*n * sizeof(*cov));
  
  info = pcapack_cov(COR_PEARSON, m, n, x, cov);
  if (info != 0) goto cleanup;
  
  tmp = 1. - 1./((double) m);
  dscal_(&(int){n*n}, &tmp, cov, &(int){1});
  
  // Take eigen decomposition
  info = pcapack_eig(true, false, true, n, cov, sdev, rotation);

  if (info != 0) goto cleanup;

  // sdev = rev(sqrt(sdev))
  sqrt_rev(n, sdev);
  
  pcapack_xpose(n, n, rotation);
  
  if (retx)
    dgemm_(&trans, &trans, &m, &n, &n, &(double){1.}, x, &m, rotation, &n, &(double){0.}, x, &m);
  
  cleanup:
    free(cov);
    if (retx) free(x_cp);
  
  return info;
}

