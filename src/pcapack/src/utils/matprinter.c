#include <stdio.h>

void matprinter(int m, int n, double *x)
{
  int i, j;
  
  for (i=0; i<m; i++)
  {
    for (j=0; j<n; j++)
      printf("%f ", x[i+m*j]);
    
    putchar('\n');
  }
}
