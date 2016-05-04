#include "../misc.h"
#include "lapack.h"

int get_qr(const int m, const int n, double *x, double *tau, double *work, const int lwork)
{
  int info = 0;
  const int rank = MIN(m, n);
  
  // q = q part or the qr of y
  dgeqrf_(&m, &n, x, &m, tau, work, &lwork, &info);
  if (info != 0) return info;
  dorgqr_(&m, &n, &rank, x, &m, tau, work, &lwork, &info);
  
  return info;
}

