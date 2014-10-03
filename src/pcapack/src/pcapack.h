/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2014, Schmidt

#ifndef PCAPACK_H
#define PCAPACK_H


#include <stdlib.h>
#include <stdbool.h>


// fastmap_mod.f90
void cma_(int *n, int *p, double *x, int *k, int *info);

// pca_mod.f90
void prcomp_svd_(int *m, int *n, double *x, double *sdev, 
  double *rotation, bool *retrot, bool *centerx, bool *scalex, 
  int *info);

//// svd_mod.f90
//void LA_svd_(int *nu, int *nv, int *m, int *n, double *x, double *s, 
//  double *u, double *vt, int *info);


#endif
