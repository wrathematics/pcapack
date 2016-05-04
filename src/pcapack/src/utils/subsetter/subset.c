// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt


#include <stdlib.h>

#include "matrix.h"


int pcapack_subsetter_retlen(matrix_t &x, intvec_t &rows, intvec_t &cols, matrix_t &ret)
{
  
}

int pcapack_subsetter_noalloc(matrix_t &x, intvec_t &rows, intvec_t &cols, matrix_t &ret)
{
  int i, j;
  
  
  for (j=0; j<ncol; j++)
  {
    for (i=0; i<nrow; i++)
    {
      
    }
  }
  
  return 0;
}


