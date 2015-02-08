#include <stdbool.h>


static inline void swap_ind(int i, int j, double *x)
{
  double tmp;
  tmp = x[i];
  x[i] = x[j];
  x[j] = tmp;
}



// TODO implement a more efficient version eventually...
static void xpose_inplace(const int m, const int n, double *x)
{
  int i, j, k, idx;
  
  for(k=0; k<n*m; k++)
  {
    idx = k;
    
    while (true)
    {
      idx = (idx % n)*m + idx/n;
      if (idx >= k) break;
    }
    
    swap_ind(k, idx, x);
  }
  
  return;
}



static void xpose_square(int m, int n, double *x)
{
  int i, j;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<j; i++)
      swap_ind(i, j, x);
  }
  
  return;
}



void pcapack_xpose(int m, int n, double *x)
{
  if (m == n)
    xpose_square(m, n, x);
  else
    xpose_inplace(m, n, x);
  
  return;
}

