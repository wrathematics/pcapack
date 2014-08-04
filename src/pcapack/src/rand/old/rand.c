/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */


// Copyright 2013, Schmidt and Heckendorf

// This contains a subset of random number generators and samplers
// from a forthcoming HPC library.

#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "../include/rand.h"


// ######################## //
// Random number generators //
// ######################## //

// uniform
double runif1(const double min, const double max)
{
  
  double unif;
  
  // Check that the RNG has been initialized; if not, start it
  mt_check();
  
  unif = min + (max-min)*((double) mt_extract())*I32BIT;
  
  return unif;
}



// sample from min/max
int sample1(const int min, const int max)
{
  int ret;
  
  // Check that the RNG has been initialized; if not, start it
  mt_check();
  
  ret = (int) runif1(min, max+1);
  
  return ret;
}


void sample1_(int *min, int *max, int *ret)
{
  *ret = sample1(*min, *max);
}


