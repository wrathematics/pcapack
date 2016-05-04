// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <string.h>

#include "sumstats.h"
#include "../utils/rank.h"
#include "../matlib/matlib.h"


// centering is done in-place

/**
 * @file
 * @brief Covariance.
 *
 * @details
 * Computes the variance-covariance matrix.  Centering is done in-place.
 * 
 * @param method
 * Input.  The form the covariance matrix takes (pearson, kendall, 
 * spearman).  Currently only pearson works.
 * @param m,n
 * Inputs.  Problem size (dims of x)
 * @param x
 * Input.  The data matrix.
 * @param coc
 * Output.  The covariance matrix.
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_cov(const int method, int m, int n, double *restrict x, double *restrict cov)
{
  int info = 0;
  const double alpha = 1. / ((double) m-1);
  
  pcapack_scale(true, false, m, n, x);
  
  info = pcapack_crossprod(true, m, n, x, alpha, cov);
  
  return info;
}



static inline double covar(const int n, const double *x, const double *y)
{
  const double recip_n = (double) 1. / (n-1);
  double sum_xy = 0., sum_x = 0., sum_y = 0.;
  double tx, ty;
  
  #pragma omp simd reduction(+: sum_xy, sum_x, sum_y)
  for (int i=0; i<n; i++)
  {
    tx = x[i];
    ty = y[i];
    
    sum_xy += tx*ty;
    sum_x += tx;
    sum_y += ty;
  }
  
  return (sum_xy - (sum_x*sum_y*((double) 1./n))) * recip_n;
}

int pcapack_cov_naive(const int m, const int n, const double *restrict x, double *restrict cov)
{
  int i, j;
  
  #pragma omp parallel for
  for (j=0; j<n; j++)
  {
    for (i=j; i<n; i++)
      cov[i + n*j] = covar(m, x+m*j, x+m*i);
  }
  
  int info = pcapack_symmetrize(LOWER, n, n, cov);
  
  return info;
}

