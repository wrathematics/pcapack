/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2015, Schmidt

#ifndef PCAPACK_H
#define PCAPACK_H


#include <stdlib.h>
#include <stdbool.h>

#include "sumstats/sumstats.h"
#include "utils/utils.h"


// pca.c
int pcapack_prcomp_svd(bool centerx, bool scalex, bool retrot, int m, int n, double *x, double *sdev, double *rotation);
int pcapack_prcomp_eigcov(bool retrot, int m, int n, double *x, double *sdev, double *rotation);


// svd.c
int pcapack_svd(bool inplace, const int nu, const int nv, int m, int n, double *x, double *s, double *u, double *vt);
int pcapack_eig(bool inplace, bool only_values, bool symmetric, int n, double *x, double * values, double *vectors);


#endif
