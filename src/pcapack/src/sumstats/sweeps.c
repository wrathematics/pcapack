// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "sumstats.h"

#define OMP_MIN_SIZE 2500


// sweep array out of matrix in-place
int pcapack_sweep(const int m, const int n, double *x, double *vec, int lvec, int margin, int fun)
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
      for (j=0; j<n; j++)
      {
        for (i=0; i<m; i++)
          x[i + m*j] += vec[i];
      }
    }
    else if (fun == MINUS)
    {
      for (j=0; j<n; j++)
      {
        for (i=0; i<m; i++)
          x[i + m*j] -= vec[i];
      }
    }
    else if (fun == TIMES)
    {
      for (j=0; j<n; j++)
      {
        for (i=0; i<m; i++)
          x[i + m*j] *= vec[i];
      }
    }
    else if (fun == DIVIDE)
    {
      for (j=0; j<n; j++)
      {
        for (i=0; i<m; i++)
          x[i + m*j] /= vec[i];
      }
    }
  }
  
  else if (margin == COLS && lvec == n)
  {
    if (fun == PLUS)
    {
      for (j=0; j<n; j++)
      {
        tmp = vec[j];
        for (i=0; i<m; i++)
          x[i + m*j] += tmp;
      }
    }
    else if (fun == MINUS)
    {
      for (j=0; j<n; j++)
      {
        tmp = vec[j];
        for (i=0; i<m; i++)
          x[i + m*j] -= tmp;
      }
    }
    else if (fun == TIMES)
    {
      for (j=0; j<n; j++)
      {
        tmp = vec[j];
        for (i=0; i<m; i++)
          x[i + m*j] *= tmp;
      }
    }
    if (fun == DIVIDE)
    {
      for (j=0; j<n; j++)
      {
        tmp = vec[j];
        for (i=0; i<m; i++)
          x[i + m*j] /= tmp;
      }
    }
  }
  
  // General case
  else
  {
    if (fun == PLUS)
    {
      if (margin == ROWS)
      {
        pos = 1;
        for (j=0; j<n; j++)
        {
          for (i=0; i<m; i++)
          {
            x[i + m*j] += vec[pos];
            pos = (pos+1) % lvec;
          }
        }
      }
      else if (margin == COLS)
      {
        for (j=0; j<n; j++)
        {
          pos = j%lvec;
          for (i=0; i<m; i++)
          {
            x[i + m*j] += vec[pos];
            pos = (pos+n) % lvec;
          }
        }
      }
    }
    else if (fun == MINUS)
    {
      
    }
    // TODO etc ...
  }
}



// in-place
int pcapack_scale(bool centerx, bool scalex, const int m, const int n, double *x)
{
  int i, j;
  double colmean, colvar;
  double dt, tmp;
  
  if (m == 0 || n == 0) return 0;
  
  // Doing both at once, if needed, is more performant
  if (centerx && scalex)
  {
    for (j=0; j<n; j++)
    {
      colmean = 0;
      colvar = 0;
      
      for (i=0; i<m; i++)
      {
        dt = x[i + m*j] - colmean;
        colmean += dt/((double) i+1);
        colvar += dt * (x[i + m*j] - colmean);
      }
      
      colvar = sqrt(colvar / ((double) m-1));
      
      // Remove mean and variance
      for (i=0; i<m; i++)
        x[i + m*j] = (x[i   + m*j]- colmean) / colvar;
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
      for (i=0; i<m; i++)
        colmean += x[i   + m*j] * div;
      
      // Remove mean from column
      for (i=0; i<m; i++)
        x[i   + m*j] -= colmean;
      
    }
  }
  else if (scalex) // RMSE
  {
    const double div = 1./((double) m-1);
    for (j=0; j<n; j++)
    {
      colvar = 0;
      
      // Get column variance
      #pragma omp parallel for private(i, tmp) shared(j, x) if(m > OMP_MIN_SIZE) reduction(+:colvar)
      for (i=0; i<m; i++)
      {
        tmp = x[i + m*j];
        colvar += tmp*tmp*div;
      }
      
      colvar = sqrt(colvar);
      
      // Remove variance from column
      #pragma omp parallel for private(i) shared(j, x, colvar)
      for (i=0; i<m; i++)
        x[i + m*j] /= colvar;
    }
  }
  
  return 0;
}

