// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "sumstats.h"
#include "../omp.h"


int pcapack_scale(const bool centerx, const bool scalex, const int m, const int n, double *restrict x)
{
  int i, j;
  double colmean, colvar;
  double dt, tmp;
  
  if (m == 0 || n == 0) return 0;
  
  ALIGNMENT(x, 16);
  
  // Doing both at once, if needed, is more performant
  if (centerx && scalex)
  {
    tmp = 1. / ((double) m-1);
    
    #pragma omp parallel for private(i, j, colmean, colvar, dt) shared(x) if(m*n > OMP_MIN_SIZE)
    for (j=0; j<n; j++)
    {
      colmean = 0;
      colvar = 0;
      
      SAFE_SIMD
      for (i=0; i<m; i++)
      {
        dt = x[i + m*j] - colmean;
        colmean += dt/((double) i+1);
        colvar += dt * (x[i + m*j] - colmean);
      }
      
      colvar = sqrt(colvar * tmp);
      
      // Remove mean and variance
      SAFE_SIMD
      for (i=0; i<m; i++)
        x[i + m*j] = (x[i + m*j] - colmean) / colvar;
    }
  }
  else if (centerx)
  {
    const double div = 1. / ((double) m);
    
    #pragma omp parallel for private(i, j, colmean) shared(x) if(m*n > OMP_MIN_SIZE)
    for (j=0; j<n; j++)
    {
      colmean = 0;
      
      // Get column mean
      SAFE_SIMD
      for (i=0; i<m; i++)
        colmean += x[i   + m*j] * div;
      
      // Remove mean from column
      SAFE_SIMD
      for (i=0; i<m; i++)
        x[i   + m*j] -= colmean;
      
    }
  }
  else if (scalex) // RMSE
  {
    const double div = 1./((double) m-1);
    
    #pragma omp parallel for private(i, j, colvar, tmp) shared(x) if (m*n > OMP_MIN_SIZE)
    for (j=0; j<n; j++)
    {
      colvar = 0;
      
      // Get column variance
      SAFE_SIMD
      for (i=0; i<m; i++)
      {
        tmp = x[i + m*j];
        colvar += tmp*tmp*div;
      }
      
      colvar = sqrt(colvar);
      
      // Remove variance from column
      SAFE_SIMD
      for (i=0; i<m; i++)
        x[i + m*j] /= colvar;
    }
  }
  
  return 0;
}

