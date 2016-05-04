// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <math.h>

#include "matlib.h"
#include "lapack.h"
#include "../omp.h"


static inline int idamax(const int n, const double *restrict x)
{
  int i;
  int ret = 0;
  double max = fabs(x[0]);
  
  #pragma omp parallel for simd
  for (i=1; i<n; i++)
  {
    if (fabs(x[i]) > max)
    {
      max = x[i];
      ret = i;
    }
  }
  return ret;
}

static inline int idamin(const int n, const double *restrict x)
{
  int i;
  int ret = 0;
  double min = fabs(x[0]);
  
  #pragma omp parallel for simd default(shared)
  for (i=1; i<n; i++)
  {
    if (fabs(x[i]) < min)
    {
      min = x[i];
      ret = i;
    }
  }
  return ret;
}

// Modification of DNRM2 for more general Minkowski 
// Original code is:
// This version written on 25-October-1982.
//     Modified on 14-October-1993 to inline the call to DLASSQ.
//     Sven Hammarling, Nag Ltd.
// Modified code by Drew Schmidt, 2015
static double dnrm3(const int n, double *restrict x, const int p)
{
  int i;
  const double dp = (double) p;
  double absxi;
  double scale = 0., ssq = 1.;
  double nrm = 0.;
  
  if (n < 1) 
    return 0.;
  else if (n == 1) 
    return fabs(x[0]);
  

  for (i=0; i<n; i++)
  {
    if (x[i] != 0)
    {
      absxi = fabs(x[i]);
      if (scale < absxi)
      {
        ssq = 1. + ssq*pow(scale/absxi, dp);
        scale = absxi;
      }
      else
        ssq += pow(absxi/scale, dp);
    }
  }
  
  nrm = scale * pow(ssq, 1./dp);
  return nrm;
}



double pcapack_distance(const int method, const int n, const double *restrict x, const double *restrict y, const int p)
{
  double dist = 0;
  double tmp;
  int i;
  int ux, uy;
  double *work;
  
  int i1 = 1;
  double neg1 = -1.;
  
  work = malloc(n * sizeof(*work));
  
  
  SAFE_PARALLEL_FOR_SIMD
  for (i=0; i<n; i++)
    work[i] = x[i] - y[i];
  
  
  if (method == DIST_EUCLIDEAN)
    dist = dnrm2(n, work, 1);
  else if (method == DIST_SUPREMUM)
    dist = work[idamax(n, work)];
  else if (method == DIST_INFIMUM)
    dist = work[idamin(n, work)];
  else if (method == DIST_MANHATTAN)
  {
    #pragma omp parallel for simd default(shared) reduction(+:dist)
    for (i=0; i<n; i++)
      dist += fabs(work[i]);
  }
  else if (method == DIST_CANBERRA)
  {
    #pragma omp parallel for simd default(shared) reduction(+:dist)
    for (i=0; i<n; i++)
    {
      tmp = fabs(x[i]) + fabs(y[i]);
      if (tmp != 0.)
        dist += fabs(work[i]) / tmp;
    }
  }
  else if (method == DIST_RCANBERRA)
  {
    #pragma omp parallel for simd default(shared) reduction(+:dist)
    for (i=0; i<n; i++)
    {
      tmp = fabs(x[i] + y[i]);
      if (tmp != 0.)
        dist += fabs(work[i]) / tmp;
    }
  }
  else if (method == DIST_MINKOWSKI)
    dist = dnrm3(n, work, p);
  // The distance is the _proportion_ of bits in which only one is
  // on amongst those in which at least one is on.
  else if (method == DIST_BINARY)
  {
    int on = 0;
    
    #pragma omp parallel for simd default(shared) reduction(+:dist,on)
    for (i=0; i<n; i++)
    {
      if (x[i] || y[i])
      {
        on++;
        if (!(x[i] && y[i]))
          dist += 1.;
      }
    }
    
    if (on != 0)
      dist /= (double) on;
  }
  
  free(work);
  return dist;
}



double pcapack_mahalanobis(int n, const double *restrict x, const double *restrict y, const double *restrict covinv)
{
  int i;
  double mahal = 0.;
  double *work1, *work2;
  
  work1 = malloc(n * sizeof(*work1));
  work2 = malloc(n * sizeof(*work2));
  
  SAFE_PARALLEL_FOR_SIMD
  for (i=0; i<n; i++)
    work1[i] = x[i] - y[i];
  
  // mahal = work1^t * covinv * work2
  dgemv_(&(char){'n'}, &n, &n, &(double){1.}, covinv, &n, work1, &(int){1}, &(double){0.}, work2, &(int){1});
  
  #pragma omp parallel for simd default(shared) reduction(+:mahal)
  for (i=0; i<n; i++)
    mahal += work1[i]*work2[i];
  
  mahal = sqrt(mahal);
  
  free(work1);
  free(work2);
  return mahal;
}

