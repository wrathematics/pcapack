#include <stdlib.h>
#include <stdio.h>


void LA_svd_(int *nu, int *nv, int *m, int *n, double *x, double *s, 
  double *u, double *vt, int *info);


int main()
{
  int nu = 3, nv = 3, m = 10, n = 3, info = 0;
  double *x = malloc(m*n * sizeof(double));
  double *s = malloc(n * sizeof(double));
  double *u = malloc(m*n * sizeof(double));
  double *vt = malloc(n*n * sizeof(double));
  
  LA_svd_(&nu, &nv, &m, &n, x, s, u, vt, &info);
  
  
  
  free(x);
  free(s);
  free(u);
  free(vt);
  
  return 0;
}
