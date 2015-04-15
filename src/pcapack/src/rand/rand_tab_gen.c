/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

// Copyright 2013, Schmidt

#include <math.h>
#include "rand.h"

// Fill the lookup tables needed for the ziggurat algorithm
// Included only for reproducability purposes

double zig_tab_x[ZT_SIZE];
double zig_tab_w[ZT_SIZE];
unsigned long zig_tab_k[ZT_SIZE];


void gen_zig_tabs()
{
  const double r = 3.6541528853610088;
  const double v = 0.00492867323399;
  int i;
  double tmp;

  // x
  zig_tab_x[ZT_SIZE-1] = r;

  for (i = ZT_SIZE-2; i > 0; i--)
    zig_tab_x[i] = sqrt(-2.0*log(v/zig_tab_x[i+1] + std_norm_pdf(zig_tab_x[i+1])));

  // w
  zig_tab_w[0] = pow(0.5, 32.0) * v / std_norm_pdf(r);

  for (i = 1; i < ZT_SIZE; i++)
    zig_tab_w[i] = pow(0.5, 32.0) * zig_tab_x[i];

  // k
  zig_tab_k[0] = pow(2.0, 32.0) * r * std_norm_pdf(r) / v;

  for (i = 1; i < ZT_SIZE; i++)
    zig_tab_k[i] = pow(2.0, 32.0) * zig_tab_x[i-1] / zig_tab_x[i];
}

