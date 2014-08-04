/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "../include/rand.h"
#include "rand_tab.h"
/*#include "../../include/rand_mt.h"*/
/*#include "../../include/rand_uniform.h"*/
/*#include "../../include/shuffler.h"*/


// References
    // http://www.jstatsoft.org/v05/i08



// probability density function for the standard normal distribution
static inline double std_norm_pdf(const double x)
{
  return exp(-x*x*0.5);;
}


// Generate a value from the standard normal pdf using the Ziggurat Algorithm
// Commented numbers refer to those from the reference paper.
double stdnorm()
{
  unsigned int i, j;
  double x, y, U;
  
  // Check that the RNG has been initialized; if not, start it
  mt_check();
  
  while( 1 )
  {
    // 1
    j = mt_extract();
    i = j & 0x000000FF;
    
    // 2
    x = j * zig_tab_w[i];
    
    if (j < zig_tab_k[i])
      return x;
    
    // 3
    else if (i == 0)
    {
      y = -log(j * I32BIT);
      U = runif1(0, 1);
      if (-2.0*log(U) > y*y)
      {
        x += zig_tab_w[1];
        return x;
      }
    }
    
    // 4
    else 
    {
      U = runif1(0, 1);
      y = ( snp_zig_tab_x[i-1] - snp_zig_tab_x[i] )*U + snp_zig_tab_x[i];
      
      if (y < std_norm_pdf(x))
        return x;
    }
    
    // 5 --- return to 1
  }
}





// normal
double rnorm(const double mn, const double sd)
{
  double N;
  
  // Check that the RNG has been initialized; if not, start it
  mt_check();
  
  N = sd*stdnorm() + mn;
  
  return N;
}


double* rnorm_arr(const int n, const double mn, const double sd)
{
  int i, cf, sign[2] = {-1,1};
  double* nrml;
  
  nrml = malloc(n * sizeof(double));
  
  // Check that the RNG has been initialized; if not, start it
  mt_check();
  
  for (i = 0; i < n; i++)
  {
    cf = mt_extract() % 2;
    nrml[i] = sd * sign[cf] * stdnorm() + mn;
  }
  
  return nrml;
}



void rnormn_(int *n, double *mn, double *sd, double *ret)
{
  int i, cf, sign[2] = {-1,1};
  
  // Check that the RNG has been initialized; if not, start it
  mt_check();
  
  for (i = 0; i < *n; i++)
  {
    cf = mt_extract() % 2;
    ret[i] = (*sd) * sign[cf] * stdnorm() + (*mn);
  }
}
