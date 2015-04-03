#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "pcapack.h"
#include "misc.h"
#include "lapack.h"


int pcapack_svd(bool inplace, const int nu, const int nv, int m, int n, double *restrict x, double *restrict s, double *restrict u, double *restrict vt)
{
  char jobz;
  int i;
  int info = 0;
  int lwork, *iwork;
  double tmp, *work, *x_cp;
  int minmn = m<n ? m : n;
  
  if (nu == 0 && nv == 0)
    jobz = 'n';
  else if ((nu == 0 && m >= n) || (nv == 0 && m < n))
    jobz = 'o';
  else if (nu <= minmn && nv <= minmn)
    jobz = 's';
  else
    jobz = 'a';
  
  
  if (likely(inplace))
    x_cp = x;
  else
  {
    x_cp = malloc(m*n * sizeof(*x_cp));
    memcpy(x_cp, x, m*n*sizeof(*x_cp));
  }
  
  iwork = malloc(8*minmn * sizeof(*iwork));
  
  lwork = -1;
  dgesdd_(&jobz, &m, &n, x_cp, &m, s, u, &m, vt, &minmn, &tmp, &lwork, iwork, &info);
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  dgesdd_(&jobz, &m, &n, x_cp, &m, s, u, &m, vt, &minmn, work, &lwork, iwork, &info);
  
  if (likely(!inplace)) free(x_cp);
  free(work);
  free(iwork);
  
  return info;
}



// m == n
int pcapack_eig(bool inplace, bool only_values, bool symmetric, int n, double *restrict x, double *restrict values, double *restrict vectors)
{
  int info = 0;
  double *x_cp;
  // Passed to Fortran
  double worksize;
  int lwork, liwork;
  int *iwork;
  double *work;
  static int neg1 = -1;
  
  
  /*if (!inplace)*/
  /*{*/
    /*x_cp = malloc(n*n * sizeof(*x_cp));*/
    /*memcpy(x_cp, x, n*n*sizeof(double));*/
  /*}*/
  /*else*/
    /*x_cp = x;*/
  if (inplace)
    x_cp = x;
  else
  {
    memcpy(vectors, x, n*n*sizeof(double));
    x_cp = vectors;
  }
  
  
  if (likely(symmetric))
  {
    char uplo = 'u';
    char jobz;
    
    if (only_values)
      jobz = 'n';
    else
      jobz = 'v';
    
    dsyevd_(&jobz, &uplo, &n, x_cp, &n, values, &worksize, &neg1, &liwork, &neg1, &info);

    lwork = (int) worksize;
    work = malloc(lwork * sizeof(*work));
    iwork = malloc(liwork * sizeof(*iwork));
    
    dsyevd_(&jobz, &uplo, &n, x_cp, &n, values, work, &lwork, iwork, &liwork, &info);

    free(work);
    free(iwork);
  }
  #if 0 // TODO
  else
  {
    char jobvl = 'n', jobvr;
    
    if (only_values)
      jobvr = 'v';
    else
      jobvr = 'n';
    
    
    dgeev_(jobvl, jobvr, &n, x, &m, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info);
    
    
    dgeev_(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, &neg1, info);
  }
  #endif
  
  cleanup:
    if (!inplace) free(x_cp);

  return info;
}
/*dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info);*/

