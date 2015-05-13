#include <stdbool.h>


// TODO implement a more efficient version eventually...
static inline void swap_by_ind(int i, int j, double *x)
{
  double tmp;
  tmp = x[i];
  x[i] = x[j];
  x[j] = tmp;
}

static void xpose_inplace_nonsquare(const int m, const int n, double *x)
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
    
    swap_by_ind(k, idx, x);
  }
  
  return;
}



static void xpose_inplace_square(const int n, double *x)
{
  int i, j;
  double tmp;
  
  for (j=0; j<n; j++)
  {
    for (i=0; i<j; i++)
    {
      tmp = x[i + n*j];
      x[i + n*j] = x[j + n*i];
      x[j + n*i] = tmp;
    }
  }
  
  return;
}



void pcapack_xpose(const int m, const int n, double *x)
{
  if (m == n)
    xpose_inplace_square(n, x);
  else
    xpose_inplace_nonsquare(m, n, x);
  
  return;
}

