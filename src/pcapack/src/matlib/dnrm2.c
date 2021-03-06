// Translation of dnrm2 to C

#include <math.h>
#include "lapack.h"

double dnrm2(int n, double *restrict x, int incx)
{
  int ix;
  double absxi, norm, scale, ssq;
  
  
  if (n < 1 || incx < 1)
    norm = 0.;
  else if (n == 1)
    norm = fabs(x[0]);
  else
  {
    scale = 0.;
    ssq = 1.;
    
    dlassq_(&n, x, &incx, &scale, &ssq);
    norm = scale * sqrt(ssq);
  }
  
  return norm;
}

