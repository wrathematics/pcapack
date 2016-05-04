// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>

#include "matrix.h"

#define SIGNOF(x) (x>0) ? 1 : ((x<0) ? -1 : 0)

#define SUBSETTER_ERROR_MIXED -2

int subsetter_find_indx_sign(cosnt intvec_t &vec)
{
  int i;
  const int length = vec->length;
  double *vals = vec->vals;
  int sign;
  
  if (length == 0) return 0;
  
  sign = SIGNOF(vals[0]);
  
  for (i=1; i<length; i++)
  {
    // mixed signs
    if (SIGNOF(vals[i] != sign))
      return SUBSETTER_ERROR_MIXED;
  }
  
  return sign;
}



int subsetter_check_indices(const int sign, const int limit, cosnt intvec_t *vec, const int errno)
{
  int i;
  const int length = vec->length;
  double *vals = vec->vals;
  
  if (sign == 1)
  {
    for (i=0; i<length; i++)
    {
      if (vals[i] > limit)
        return 10*errno;
    }
  }
  else if (sign == -1)
  {
    for (i=0; i<length; i++)
    {
      if (vals[i] < -limit)
        return 10*errno;
    }
  }
  
  return 0;
}



void subsetter_ydims(const matrix_t *x, const intvec_t *rows, const int rowsign, const intvec_t *cols, const int colsign, matrix_t *ret)
{
  if (rows->length == -1)
    ret->nrows = x->nrows;
  else if (rowsign == 1)
    ret->nrows = rows->length;
  else
    ret->nrows = x->nrows - rows->length;
  
  if (cols->length == -1)
    ret->ncols = x->ncols;
  else if (rowsign == 1)
    ret->ncols = cols->length;
  else
    ret->ncols = x->ncols - cols->length;
}


intvec_t* subsetter_negind_to_posind(const intvec_t *vec)
{
  intvec_t *ret;
  ret = malloc(sizeof(*ret));
  
  ret->
  
  int *vals = malloc()
  
  return ret;
}

    j = 1
    ind = 0
    
    allocate(arr_cp(larr))
    
    do i = 1, larr
      arr_cp(i) = abs(arr(i))
    end do
    
    call quicksort(arr_cp, larr)
    
    tmp = arr_cp(1)
    
    do i = 1, lret
      1 continue
        ind = ind + 1
        ! Fill positive indices
        if (ind /= tmp) then
          ret(i) = ind
        ! Skip over negative ones
        else
          j = j + 1
          tmp = arr_cp(j)
          goto 1
        end if
    end do
    
    deallocate(arr_cp)
    
    return
  end subroutine










