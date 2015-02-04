#include <stdlib.h>
#include "pcapack.h"


int pcapack_svd(const int nu, const int nv, int m, int n, const double *x, double *s, double *u, double *vt)
{
  char jobz;
  int i;
  int info = 0;
  int lwork, *iwork;
  double tmp, *work, *cpx;
  const minmn = m<n ? m : n;
  
  if (nu == 0 && nv == 0)
    jobz = 'n';
  else if ((nu == 0 && m >= n) || (nv == 0 && m < n))
    jobz = 'o';
  else if (nu <= minmn && nv <= minmn)
    jobz = 's';
  else
    jobz = 'a';
  
  
  cpx = malloc(m*n * sizeof(*cpx));
  for (i=0; i<m*n; i++)
    cpx[i] = x[i];
  
  iwork = malloc(8*minmn * sizeof(*iwork));
  
  lwork = -1;
  dgesdd_(&jobz, &m, &n, cpx, &m, s, u, &m, vt, &minmn, &tmp, &lwork, iwork, &info);
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  dgesdd_(&jobz, &m, &n, cpx, &m, s, u, &m, vt, &minmn, work, &lwork, iwork, &info);
  
  free(cpx);
  free(work);
  free(iwork);
  
  return info;
}


