// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <math.h>

double pcapack_variance(const int n, double *x)
{
  int i;
  double mean, var, dt;
  
  mean = 0;
  var = 0;
  
  for (i=0; i<n; i++)
  {
    dt = x[i] - mean;
    mean += dt / ((double) i);
    var += dt*(x[i] - mean);
  }
  
  var /= ((double) n-1);
  
  return var;
}



double pcapack_sdev(const int n, double *x)
{
  double sdev = pcapack_variance(n, x);
  
  return sqrt(sdev);
}



int pcapack_rowvars(const int m, const int n, double *x, double *rowvars)
{
  int i, j;
  
  
}



int pcapack_colvars(const int m, const int n, double *x, double *colvars)
{
  int i, j;
  
  for (j=0; j<n; j++)
    colvars[j] = pcapack_variance(m, x+(m*j));
  
  return 0;
}

// sdevs



int pcapack_colsdev(const int m, const int n, double *x, double *colsdev)
{
  int i, j;
  
  for (j=0; j<n; j++)
    colsdev[j] = pcapack_sdev(m, x+(m*j));
  
  return 0;
}

