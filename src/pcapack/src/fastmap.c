// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt


#include "pcapack.h"


static inline double sign(double x)
{
  if (x > 0.)
    return 1.;
  else if (x < 0.)
    return -1.;
  else
    return 0.;
}



#define PCAPACK_OOM -2147483647


/**
 * @file
 * @brief Coordinate Mapping Algorithm
 *
 * @details
 * Implementation of the fastmap coordinate mapping algorithm, 
 * based on the paper:
 * "A Matrix Computation View of Fastmap and Robustmap Dimension 
 * Reduction Algorithms", George Ostrouchov, SIAM J. Matrix Anal. 
 * appl.
 * 
 * @param n,p
 * Inputs.  Problem size (dims of x)
 * @param x
 * In/Output.  The data matrix.  
 * @param k
 * TODO
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_cma(const int n, const int p, double *restrict x, const int k)
{
  char trans = 'n';
  int info = 0;
  int i, j, l;
  int ncol = p;
  int intone = 1;
  double one = 1., negone = -1., negtwo = -2., zero = 0.;
  double tmp, best;
  double *a, *b, *y, *work;
  rng_state_t rs;
  
  if (n < 1 || p < 1) return 0;
  if (k > n) return -1;
  if (k > p) return -2;
  
  a = malloc(p * sizeof(*a));
  b = malloc(p * sizeof(*b));
  y = malloc(n * sizeof(*y));
  work = malloc(p * sizeof(*work));
  
  
  // TODO seed handling
  rng_prepare(&rs);
  rng_set_type(&rs, RNG_TYPE_MT);
  rng_init(&rs, 1234);
  
  
  for (i=0; i<k; i++)
  {
    // Select pivot row pair
    pcapack_fastmap(&rs, n, ncol, x+(n*i), a, b, work);
    
    // Translate rows to pivot origin
    pcapack_sweep(n, ncol, x+(n*i), a, ncol, COLS, MINUS);
    
    // Apply householder reflection on right
    daxpy_(&ncol, &negone, a, &intone, b, &intone);
    
    b[0] += sign(b[0]) * dnrm2(ncol, b, 1);
    
    tmp = dnrm2(p, b, 1);
    for (j=0; j<ncol; j++)
      b[j] /= tmp;
    
    dgemv_(&trans, &n, &ncol, &one, x+(n*i), &n, b, &intone, &zero, y, &intone);
    
    dger_(&n, &ncol, &negtwo, y, &intone, b, &intone, x+(n*i), &n);
    
    ncol--;
  }
  
  
  free(work);
  free(y);
  free(b);
  free(a);
  
  return info;
}



// Fastmap pivot selection
    // Working on the right ncol columns of x
    // a is the point (in x) furthest from a random vector in x
    // b is the point (in x) furthest from point a
    // work is a workspace array
static inline void bestdist(int n, int ncol, double *restrict x, double *restrict a, double *restrict b, double *restrict work)
{
  int i, j;
  int ia = 0;
  int intone = 1;
  double tmp, best;
  double one = 1.;
  
  best = 0.;
  
  for (i=0; i<n; i++)
  {
    for (j=0; j<ncol; j++)
      work[j] = -1. * x[i + n*j];
    
    daxpy_(&ncol, &one, b, &intone, work, &intone);
    tmp = dnrm2(ncol, work, 1);
    
    if (tmp > best)
    {
      best = tmp;
      ia = i;
    }
  }
  
  for (j=0; j<ncol; j++)
    a[j] = x[ia + n*j];
  
}


// x[n, ncol]
void pcapack_fastmap(rng_state_t *rs, int n, int ncol, double *restrict x, double *restrict a, double *restrict b, double *restrict work)
{
  // a, b, work of length ncol
  int i, ia, ib;
  
  // Take random row b in x;
  ia = sample(rs, 1, n);
  
  for (i=0; i<ncol; i++)
    b[i] = x[ia + n*i];
  
  // Let a be the most distant point in x from b
  bestdist(n, ncol, x, a, b, work);
  
  // Let b be the most distant point in x from a
  bestdist(n, ncol, x, b, a, work);
}

