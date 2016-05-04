// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "sumstats.h"
#include "../omp.h"


// It ain't pretty, but function pointers are too expensive
#define SWEEP_ROWS_SPECIAL(ASSIGNMENT) \
  for (j=0; j<n; j++){ \
    for (i=0; i<m; i++){ \
      x[i + m*j] ASSIGNMENT vec[i]; }}

#define SWEEP_COLS_SPECIAL(ASSIGNMENT) \
  _Pragma("omp parallel for default(shared) private(i,j,tmp) if(m*n > OMP_MIN_SIZE)") \
  for (j=0; j<n; j++){ \
    tmp = vec[j]; \
    SAFE_SIMD \
    for (i=0; i<m; i++){ \
      x[i + m*j] ASSIGNMENT tmp; }}

#define SWEEP_ROWS(ASSIGNMENT) \
  pos = 1; \
  for (j=0; j<n; j++){ \
    for (i=0; i<m; i++){ \
      x[i + m*j] += vec[pos]; \
      pos = (pos+1) % lvec; }}

#define SWEEP_COLS(ASSIGNMENT) \
  _Pragma("omp parallel for default(shared) private(i,j,tmp) if(m*n > OMP_MIN_SIZE)") \
  for (j=0; j<n; j++){ \
    pos = j%lvec; \
    SAFE_SIMD \
    for (i=0; i<m; i++){ \
      x[i + m*j] += vec[pos]; \
      pos = (pos+n) % lvec; }}

// sweep array out of matrix in-place
int pcapack_sweep(const int m, const int n, double *restrict x, double *restrict vec, const int lvec, const int margin, const int fun)
{
  int i, j;
  int pos;
  double tmp;
  
  if (m == 0 || n == 0) return 0;
  if (margin != ROWS && margin != COLS) return -6;
  if (fun != PLUS && fun != MINUS && fun != TIMES && fun != DIVIDE) return -7;
  
  
  // Special cases --- avoids index checking
  if (margin == ROWS && lvec == m)
  {
    if (fun == PLUS)
    {
      SWEEP_ROWS_SPECIAL(+=);
    }
    else if (fun == MINUS)
    {
      SWEEP_ROWS_SPECIAL(-=);
    }
    else if (fun == TIMES)
    {
      SWEEP_ROWS_SPECIAL(*=);
    }
    else if (fun == DIVIDE)
    {
      SWEEP_ROWS_SPECIAL(/=);
    }
  }
  
  else if (margin == COLS && lvec == n)
  {
    if (fun == PLUS)
    {
      SWEEP_COLS_SPECIAL(+=);
    }
    else if (fun == MINUS)
    {
      SWEEP_COLS_SPECIAL(-=);
    }
    else if (fun == TIMES)
    {
      SWEEP_COLS_SPECIAL(*=);
    }
    if (fun == DIVIDE)
    {
      SWEEP_COLS_SPECIAL(/=);
    }
  }
  
  // General case
  else
  {
    if (fun == PLUS)
    {
      if (margin == ROWS)
      {
        SWEEP_ROWS(+=);
      }
      else if (margin == COLS)
      {
        SWEEP_COLS(+=);
      }
    }
    else if (fun == MINUS)
    {
      if (margin == ROWS)
      {
        SWEEP_ROWS(-=);
      }
      else if (margin == COLS)
      {
        SWEEP_COLS(-=);
      }
    }
    else if (fun == TIMES)
    {
      if (margin == ROWS)
      {
        SWEEP_ROWS(*=);
      }
      else if (margin == COLS)
      {
        SWEEP_COLS(*=);
      }
    }
    else if (fun == DIVIDE)
    {
      if (margin == ROWS)
      {
        SWEEP_ROWS(/=);
      }
      else if (margin == COLS)
      {
        SWEEP_COLS(/=);
      }
    }
  }
  
  return 0;
}

