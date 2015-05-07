/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#ifndef PCAPACK_H
#define PCAPACK_H


#include <stdlib.h>
#include <stdbool.h>

#include "rand/rand.h"
#include "sumstats/sumstats.h"
#include "utils/utils.h"


// fastmap.c
void pcapack_fastmap(rng_state_t *rs, int n, int ncol, double *x, double *a, double *b, double *work);
int pcapack_cma(int n, int p, double *x, int k);


// pca.c
int pcapack_prcomp_svd(const bool centerx, const bool scalex, const bool retx, const int m, const int n, double *x, double *sdev, double *rotation);
int pcapack_prcomp_eigcov(const bool retx, const int m, const int n, double *x, double *sdev, double *rotation);


// svd.c
int pcapack_svd(const bool inplace, const int nu, const int nv, const int m, const int n, double *x, double *s, double *u, double *vt);
int pcapack_eig(const bool inplace, const bool only_values, const bool symmetric, const int n, double *x, double * values, double *vectors);


#endif
