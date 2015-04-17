// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt


int pcapack_rowsums(const int m, const int n, double *x, double *rowsums)
{
  int i, j;
  double tmp;
  
  for (i=0; i<n; i++)
  {
    tmp = 0.;
    for (j=0; j<m; j++)
    {
      tmp += x[i + m*j];
    }
  }
  
  return 0;
}



int pcapack_colsums(const int m, const int n, double *x, double *colsums)
{
  int i, j;
  double tmp;
  
  #pragma omp parallel for
  for (j=0; j<m; j++)
  {
    tmp = 0.;
    
    #pragma omp simd
    for (i=0; i<n; i++)
    {
      tmp += x[i + m*j];
    }
  }
  
  return 0;
}



double pcapack_mean(const int n, double *x)
{
  int i;
  const double divbyn = 1. / ((double) n);
  double mean = 0.;
  
  #pragma omp parallel for simd
  for (i=0; i<n; i++)
    mean += x[i] * divbyn;
  
  return mean;
}



int pcapack_rowmeans(const int m, const int n, double *x, double *rowsums)
{
  int i, j;
  double tmp;
  const double divbyn = 1. / ((double) n);
  
  #pragma omp parallel for
  for (i=0; i<n; i++)
  {
    tmp = 0.;
    
    #pragma omp simd
    for (j=0; j<m; j++)
    {
      tmp += x[i + m*j] * divbyn;
    }
  }
  
  return 0;
}



int pcapack_colmeans(const int m, const int n, double *x, double *colsums)
{
  int i, j;
  double tmp;
  const double divbym = 1. / ((double) m);
  
  for (j=0; j<m; j++)
  {
    tmp = 0.;
    for (i=0; i<n; i++)
    {
      tmp += x[i + m*j] * divbym;
    }
  }
  
  return 0;
}

