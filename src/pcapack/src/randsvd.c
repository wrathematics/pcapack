// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt


#include "pcapack.h"
#include "lapack.h"
#include "misc.h"


// Input x is overwritten by its Q
static inline int get_qr(const int m, const int n, double *x, double *tau, double *work, const int lwork)
{
  int info = 0;
  const int rank = MIN(m, n);
  
  // q = q part or the qr of y
  dgeqrf_(&m, &n, x, &m, tau, work, &lwork, &info);
  if (info != 0) return info;
  dorgqr_(&m, &n, &rank, x, &m, tau, work, &lwork, &info);
  
  return info;
}



//    B <- t(Q) %*% x
//    
//    if (!compute.u)
//      nu <- 0
//    else
//      nu <- min(nrow(B), ncol(B))
//    
//    if (!compute.vt)
//      nv <- 0
//    else
//      nv <- min(nrow(B), ncol(B))
//    
//    svd.B <- La.svd(x=B, nu=nu, nv=nv)
//    call svd(nu, nv, m, n, x, s, u, vt, info)
//    
//  end subroutine



//  k        (input) integer
//           the desired number of singular vectors.
//  u        (output) double precision array
//           mxk

// s = min(m, n)
// u = mxk
// vt = 

// omega = nxl, we can take l=2*k for sizes george used
// maybe try using fastmap to create omega



#define PCAPACK_RANDSVD_RANDINIT    1
#define PCAPACK_RANDSVD_FASTMAPINIT 2
#define PCAPACK_RANDSVD_BADMETHOD   -3




/**
 * @file
 * @brief Randsvd
 *
 * @details
 * The randomized svd method from the paper Halko, Martinsson, and Tropp
 * 2011, "Finding Structure with Randomness: Probabilistic Algorithms 
 * for Constructing Approximate Matrix Decompositions", SIAM Review 
 * 53, 217-288.
 * 
 * We implement the method from the paper in 2 ways.  The first is 
 * "as is", using a random normal matrix generation for the test
 * matrix omega.  Additionally, we use the fastmap algorithm to 
 * generate the matrix omega.
 * 
 * Complexity is O(m*n*k); see the paper for details.
 *
 * @param m,n
 * Inputs.  Problem size (dims of x)
 * @param x
 * 
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
*/
int pcapack_randsvd(const bool retu, const bool retvt, const int k, const int niter, const int method, const int m, const int n, double *restrict x, double *restrict s, double *restrict u, double *restrict vt)
{
  int info = 0;
  int i, j, lwork, rc_min, rc_max;
  int nu, nv;
  double tmp;
  double *omega, *work, *Y, *tau, *B;
  const int twok = 2*k;
  
  
  if (m < 1 || n < 1) return 0;
  if (k > n || k < 1) return -7;
  if (niter < 1) return -9;
  
  if (method != PCAPACK_RANDSVD_RANDINIT && method != PCAPACK_RANDSVD_FASTMAPINIT)
    return PCAPACK_RANDSVD_BADMETHOD;
  
  rc_min = MIN(m, n);
  rc_max = MAX(m, n);
  
  // allocate and initializeomega=nx2k
  omega = malloc(n * twok * sizeof(*omega));
  
  if (method == PCAPACK_RANDSVD_RANDINIT)
  {
    // TODO use custom RNG
    for (i=0; i<n*twok; i++)
      omega[i] = rand()%2;
  }
  else if (method == PCAPACK_RANDSVD_FASTMAPINIT)
  {
    // FIXME
  }
  
  
  // Allocate workspace
  dgeqrf_(&m, &twok, Y, &m, tau, &tmp, &(int){-1}, &info);
  lwork = (int) tmp;
  
  work = malloc(lwork * sizeof(*work));
  Y = malloc(m * twok * sizeof(*Y));
  B = malloc(twok * n * sizeof(*B));
  
  // ----------------------------------------------------
  // Stage A from the paper with orthonormalization step
  // ----------------------------------------------------
  
  // Y = A*omega
  matmult(false, false, m, n, x, n, twok, omega, Y);
  
  // Y is overwritten by its QR's Q
  info = get_qr(m, n, Y, tau, work, lwork);
  if (info != 0) goto cleanup;
  
  
  // subspace iteration
  for (j=1; j<=niter; j++)
  {
    // twY_j = x^t * Q_j-1
    matmult(false, false, m, n, x, m, twok, Y, Y); // FIXME this is wrong...
    
    // qr of twY_j
    info = get_qr(m, n, Y, tau, work, lwork);
    
    // Y = x * Q
    matmult(false, false, m, n, x, m, twok, Y, Y); // FIXME this too...
    
    // Q = the Q part of the qr of Y
    info = get_qr(m, n, Y, tau, work, lwork);
    
    if (info != 0) goto cleanup;
  }
  
  
  // ----------------------------------------------------
  // Stage B from the paper
  // ----------------------------------------------------
  
  
  // B is 2KxN
  // B = t(Y) * x
  matmult(true, false, m, twok, Y, m, n, x, B);
  
  if (!retu) 
    nu = 0;
  else
    nu = MIN(twok, n);
  
  if (!retvt)
    nv = 0;
  else
    nv = MIN(twok, n);
  
  info = pcapack_svd(true, nu, nv, twok, n, B, s, u, vt);
  if (info != 0) goto cleanup;
  
  
  /*
  if (retu)
  {
    u <- svd.B$u
    u <- Q %*% U
    
    d <- svd.B$d
    
    d <- d[1L:k]
    u <- u[, 1L:k]
  }
  
  if (retvt)
  {
    vt <- svd.B$vt[1L:k, ]
  }
  */
  
  
  // ----------------------------------------------------
  // Cleanup and exit
  // ----------------------------------------------------
  
  cleanup:
    free(tau);
    free(Y);
    free(work);
    free(omega);
    free(B);
  
  return info;
}

