// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#include <stdlib.h>
#include <string.h>

#include "sumstats.h"
#include "../utils/rank.h"
#include "../matlib/matlib.h"


int pcapack_cor(const int method, int m, int n, double *x, double *cor)
{
  int info;
  int i, j;
  double *ranks;
  
  
  if (m == 0 || n == 0) return 0;
  if (method != COR_PEARSON && method != COR_SPEARMAN && method != COR_KENDALL) return -1;
  
  if (method == COR_PEARSON)
  {
    // variance(m, n, x, ret);
    // scale(false, true, n, n, ret);
  }
  else if (method == COR_SPEARMAN)
  {
    ranks = malloc(m*n * sizeof(*ranks));
    memcpy(ranks, x, m*n*sizeof(*ranks));
    
    colrank(RANK_MIN, m, n, ranks);
    pcapack_cov(method, m, n, ranks, cor);
    
    free(ranks);
  }
  else if (method == COR_KENDALL)
  {
    // TODO
  }
}

