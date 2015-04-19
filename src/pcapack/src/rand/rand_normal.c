/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "rand.h"
#include "constants.h"
#include "rand_tab.h"


// Box-Muller transform method
double rnorm_bm(rng_state_t *rs, double rate)
{
  double u, v;
  
  u = runif(rs, 0, 1);
  v = runif(rs, 0, 1);
  
  return sqrt(-2.*log(u)) * cos(2.*PI*v);
}

// Polar transform method
double rnorm_pol(rng_state_t *rs, double rate)
{
  double u, v, s;
  
  u = runif(rs, 0, 1);
  v = runif(rs, 0, 1);
  
  s = u*u + v*v;
  
  return u * sqrt(-2.*log(s) / s);
}




// References
    // http://www.jstatsoft.org/v05/i08


// probability density function for the standard normal distribution
static inline double std_norm_pdf(const double x)
{
  return exp(-x*x*0.5);;
}


// Generate a value from the standard normal pdf using the Ziggurat Algorithm
// Commented numbers refer to those from the reference paper.
double stdnorm(rng_state_t *rs)
{
  unsigned int i, j;
  double x, y, U;
  
  while(1)
  {
    // 1
    j = rng_extract(rs);
    i = j & 0x000000FF;
    
    // 2
    x = j * zig_tab_w[i];
    
    if (j < zig_tab_k[i])
      return x;
    // 3
    else if (i == 0)
    {
      y = -log(j * I32BIT);
      U = runif(rs, 0, 1);
      if (-2.0*log(U) > y*y){
        x += zig_tab_w[1];
        return x;
      }
    }
    // 4
    else 
    {
      U = runif(rs, 0, 1);
      y = ( snp_zig_tab_x[i-1] - snp_zig_tab_x[i] )*U + snp_zig_tab_x[i];
      
      if (y < std_norm_pdf(x))
        return x;
    }
    // 5 --- return to 1
  }
}





// normal
double rnorm(rng_state_t *rs, const double mn, const double sd)
{
  double N;
  
  rng_check(rs);
  
  N = sd*stdnorm(rs) + mn;
  
  return N;
}


