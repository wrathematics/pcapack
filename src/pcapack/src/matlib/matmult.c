// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdbool.h>
#include <math.h>

#include "matrix.h"


/**
 * @file
 * @brief Matrix-matrix product.
 *
 * @details
 * Sensible wrapper around dgemm.  Sets: ret_mx*ny = transx(x_mx*nx) * transy(y_mx*ny)
 *
 * @param transx,transy
 * Input.  Use transposed versions of the matrices or not?
 * @param mx,nx
 * Inputs.  Number of rows/cols of x.
 * @param my,ny
 * Input.  Number of rows/cols of y.
 * @param x,y
 * In/Output.  Matrices to multiply.
 * @param ret
 * Output.  The product of x and y.
 * 
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
void matmult(const bool transx, const bool transy, matrix_t *x, matrix_t *y, matrix_t *ret)
{
  // m = # rows of op(x)
  // n = # cols of op(y)
  // k = # cols of op(x)
  int im, in, ik;
  char ctransx, ctransy;
  static const double one = 1., zero = 0.;
  
  if (transx) ctransx = 't';
  else ctransx = 'n';
  if (transy) ctransy = 't';
  else ctransy = 'n';
  
  if (transx)
  {
    im = x->ncols;
    ik = x->nrows;
    
    if (transy == 't')
      in = y->nrows;
    else
      in = y->ncols;
  }
  else
  {
    im = x->nrows;
    ik = x->ncols;
    
    if (transy)
      in = y->nrows;
    else
      in = y->ncols;
  }
  
  dgemm_(&ctransx, &ctransy, &im, &in, &ik, &one, x->data, &(x->nrows), y->data, &(y->nrows), &zero, ret->data, &(ret->nrows));
}


