// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include "../misc.h"
#include "matlib.h"


/**
 * @file
 * @brief Matrix crossproduct.
 *
 * @details
 * Computes t(x) * x using a rank-k symmetric update.  Uses dsyrk.
 * 
 * @param symmetrize
 * Input.  Should the matrix be symmetrized? If not, only the upper triangle
 * contains useful information.
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
int pcapack_crossprod(const bool symmetrize, const int m, const int n, const double *restrict x, const double alpha, double *restrict c)
{
  int info = 0;
  
  dsyrk_(&(char){'u'}, &(char){'t'}, &n, &m, &alpha, x, &m, &(double){0.0}, c, &n);
  if (info) return info;
  
  if (symmetrize)
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
int pcapack_tcrossprod(const bool symmetrize, const int m, const int n, const double *x, const double alpha, double *c)
{
  int info = 0;
  
  dsyrk_(&(char){'u'}, &(char){'n'}, &m, &n, &alpha, x, &m, &(double){0.0}, c, &m);
  if (info) return info;
  
  if (symmetrize)
    info = pcapack_symmetrize(UPPER, m, m, c);
  
  return info;
}

