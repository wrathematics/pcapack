/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#ifndef __MATLIB_H__
#define __MATLIB_H__


#include <stdbool.h>
#include "matrix.h"


// crossprod.c
int pcapack_crossprod(const bool symmetrize, const int m, const int n, const double *restrict x, const double alpha, double *restrict c);
int pcapack_tcrossprod(const bool symmetrize, const int m, const int n, const double *restrict x, const double alpha, double *restrict c);


// distances.c
#define DIST_EUCLIDEAN  1
#define DIST_SUPREMUM   2
#define DIST_INFIMUM    3
#define DIST_MANHATTAN  4
#define DIST_CANBERRA   5
#define DIST_RCANBERRA  6
#define DIST_MINKOWSKI  7
#define DIST_BINARY     8

static double dnrm3(const int n, double *restrict x, const int p);
double pcapack_distance(const int method, const int n, const double *restrict x, const double *restrict y, const int p);
double pcapack_mahalanobis(int n, const double *restrict x, const double *restrict y, const double *restrict covinv);


// matmult.c
void matmult(const bool transx, const bool transy, matrix_t *x, matrix_t *y, matrix_t *ret);


// qr.c
int get_qr(const int m, const int n, double *x, double *tau, double *work, const int lwork);


// svd.c
int pcapack_svd(const bool inplace, const int nu, const int nv, const int m, const int n, double *restrict x, double *restrict s, double *restrict u, double *restrict vt);
int pcapack_eig(const bool inplace, const bool only_values, const bool symmetric, const int n, double *restrict x, double *restrict values, double *restrict vectors);


// symmetrize.c
#define UPPER 1
#define LOWER 2

bool pcapack_is_symmetric(const int m, const int n, const double *x);
int pcapack_symmetrize(const int triang, const int m, const int n, double *x);





#endif
