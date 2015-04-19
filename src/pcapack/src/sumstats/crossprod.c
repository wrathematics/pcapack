// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <string.h>
#include "sumstats.h"
#include "../misc.h"
#include "../lapack.h"


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



/**
 * @file
 * @brief Matrix crossproduct.
 *
 * @details
 * Computes t(x) * x using a rank-k symmetric update.  Uses dsyrk.
 * 
 * @param m,n
 * Inputs.  Problem size (dims of x)
 * @param x
 * Input.  The data matrix.
 * @param alpha
 * Value to multiply against output C.
 * @param c
 * Output (the crossproduct).
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_crossprod(int m, int n, const double *restrict x, double alpha, double *restrict c)
{
  int info = 0;
  
  dsyrk_(&(char){'u'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
  info = pcapack_symmetrize(UPPER, n, n, c);
  
  return info;
}



/**
 * @file
 * @brief Transpose of matrix crossproduct.
 *
 * @details
 * Computes x * t(x) using a rank-k symmetric update.  Uses dsyrk.
 * 
 * @param m,n
 * Inputs.  Problem size (dims of x)
 * @param x
 * Input.  The data matrix.
 * @param alpha
 * Value to multiply against output C.
 * @param c
 * Output (the tcrossproduct).
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_tcrossprod(int m, int n, const double *x, double alpha, double *c)
{
  int info = 0;
  
  dsyrk_(&(char){'u'}, &(char){'n'}, &m, &n, &alpha, x, &m, &(double){0.0}, c, &m);
  info = pcapack_symmetrize(UPPER, m, m, c);
  
  return info;
}

