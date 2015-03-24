#include <stdbool.h>
#include <math.h>
#include <float.h>


// Inspired by http://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/
#define SIGN(x,y) (y>=0?fabs(x):-fabs(x))

bool fequals(double x, double y)
{
  int ulpsdiff, ux, uy;
  int ulpstol = 2;
  static const double tol = sqrt(DBL_EPSILON);
  
  if (fabs(x-y) < tol)
    return true;
  
  if (SIGN(1.0, x) != SIGN(1.0, y))
  {
    // Check for +0/-0 case
    if (x == y)
      return true;
    else
      return false;
  }
  
  // Difference in ULP's
  ux = transfer(x, 1); //FIXME
  uy = transfer(y, 1); //FIXME
  
  ulpsdiff = abs(ux - uy);
  if (ulpsdiff <= ulpstol)
    return true;
  
}



bool isequal_mats(int m, int n, double *x, double *y, double tol)
{
  int i, j;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<m; i++)
    {
      if (fabs(x[i + m*j] - y[i + m*j]) > tol)
        return false;
    }
  }
  
  return true;
}
