#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "../misc.h"
#include "lapack.h"


/**
 * @file
 * @brief SVD
 *
 * @details
 * Singular value decomposition.  Uses LAPACK's dgesdd.  Right singular values
 * (if requested), are returned transposed.  Use xpose if you need them
 * non-transposed.
 * 
 * @param inplace
 * Input.  Should the svd work on allocated x matrix, or a local copy?
 * @param nu,nv
 * Input.  Number of left/right singular vectors to compute.
 * @maram m,n
 * Input.  Problem size (dim of x).
 * @param x
 * In/Output data matrix.  If inplace==true then the data is destroyed, 
 * otherwise it's not modified.
 * @param s
 * Output.  The singular values.
 * @param u
 * Output.  The left singular vectors.
 * @param vt
 * Output.  The transpose of the right singular vectors.
 *
 * @return
 * The return value indicates that status of the function.  Non-zero values
 * are errors.
 */
int pcapack_svd(const bool inplace, const int nu, const int nv, const int m, const int n, double *restrict x, double *restrict s, double *restrict u, double *restrict vt)
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

