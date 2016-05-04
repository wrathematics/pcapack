// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "../misc.h"
#include "lapack.h"
#include "matlib.h"


bool pcapack_is_symmetric(const int m, const int n, const double *restrict x)
{
  int i, j;
  const int k = MIN(m, n);
  
  if (m == 0 || n == 0) return 0;
  
  // NOTE keep this serial (cache misses)
  for (j=0; j<k; j++)
  {
    for (i=0; i<j; i++)
    {
      if (x[j + m*i] != x[i + m*j])
        return false;
    }
  }
  
  return true;
}



/**
 * @file
 * @brief Symmetrize a matrix.
 *
 * @details
 * Makes x symmetric via copying from one triangle to the other.
 * 
 * @param triang
 * The triangle to copy FROM.  So if trang is UPPER, then the upper
 * triangle is copied ONTO the lower triangle.
 * @param m,n
 * Inputs.  Problem size (dims of x)
 * @param x
 * In/Output.  The data matrix.
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_symmetrize(const int triang, const int m, const int n, double *x)
{
  int i, j;
  const int k = MIN(m, n);
  
  if (m == 0 || n == 0) return 0;
  if (triang != UPPER && triang != LOWER) return -1;
  
  // NOTE keep this serial (cache misses)
  // FIXME investigate further???
  
  // Copy upper ONTO lower
  //#pragma omp paralell for simd
  if (likely(triang == UPPER))
  {
    for (j=0; j<k; j++)
    {
      for (i=0; i<j; i++)
        x[j + m*i] = x[i + m*j];
    }
  }
  // Copy lower ONTO upper
  //#pragma omp paralell for simd
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

