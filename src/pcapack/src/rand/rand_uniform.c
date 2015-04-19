/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt and Heckendorf

#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "rand.h"


double runif(rng_state_t *rs, const double min, const double max)
{
  double ret;
  
  rng_check(rs);
  
  ret = min + (max-min)*((double) rng_extract(rs))*I32BIT;
  
  return ret;
}



void runif_arr(rng_state_t *rs, double *x, const int len, const double min, const double max)
{
  int i;
  const double mmm = max - min;
  
  rng_check(rs);
  
  
  for (i=0; i<len; i++)
    x[i] = min + mmm*((double) rng_extract(rs))*I32BIT;
}



int runif_int(rng_state_t *rs, const int min, const int max)
{
  int ret;
  const int mmm = max - min + 1;
  
  rng_check(rs);
  
  
  ret = (int) (min + mmm*((double) rng_extract(rs))*I32BIT);
  
  return ret;
}



void runif_int_arr(rng_state_t *rs, int *x, const int len, const int min, const int max)
{
  int i;
  const int mmm = max - min + 1;
  
  rng_check(rs);
  
  
  for (i=0; i<len; i++)
    x[i] = (int) (min + mmm*((double) rng_extract(rs))*I32BIT);
}



