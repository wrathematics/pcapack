// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Copyright 2015, Schmidt

#define MIN(m,n) m<n?m:n;


// make symmetric via copying from one triangle to the other.
#define UPPER 1
#define LOWER 2
int pcapcak_symmetrize(const int triang, const int m, const int n, double *x)
{
  int i, j;
  const int k = MIN(m, n); // TODO FIXME
  
  if (m == 0 || n == 0) return 0;
  if (triang != UPPER && triang != LOWER) return -1;
  
  // Copy upper ONTO lower
  if (triang == UPPER)
  {
    for (j=0; j<n; j++)
    {
      for (i=j+1; i<m; i++)
        x[j + m*i] = x[i + m*j];
    }
  }
  // Copy lower ONTO upper
  else if (triang == LOWER)
  {
    for (j=0; j<n; j++)
    {
      for (i=0; i<j-1; i++)
        x[j + m*i] = x[i + m*j];
    }
  }
  
  return 0;
}



// t(x) * x
int pcakack_crossprod(int m, int n, double *x, double alpha, double *c)
{
  int info = 0;
  int ldx, ldc;
  char nst;
  char triang = 'u';
  double zero = 0.;
  
  nst = 'n';
  ldx = n;
  ldc = m;
  
  dsyrk_(&triang, &nst, &ldc, &ldx, &alpha, x, &m, &zero, c, &ldc);
  info = pcapcak_symmetrize(triang, m, n, x);
  
  return info;
}

// x * t(x)
int pcapack_tcrossprod(int m, int n, double *x, double alpha, double *c)
{
  int info = 0;
  int ldx, ldc;
  char nst;
  char triang = 'u';
  double zero = 0.;
  
  nst = 't'
  ldx = m
  ldc = n
  
  dsyrk_(&triang, &nst, &ldc, &ldx, &alpha, x, &m, &zero, c, &ldc);
  info = pcapcak_symmetrize(triang, m, n, x);
  
  return info;
}



int pcapack_inverse(int n, double *x)
{
  int info = 0;
  int *ipiv;
  int lwork;
  double tmp;
  double *work;
  
  
  // Factor x = LU
  ipiv = malloc(n * sizeof(*ipiv));
  dgetrf_(&n, &n, x, &n, ipiv, &info);
  if (info != 0) goto cleanup;
  
  
  // Invert
  lwork = -1;
  dgetri_(&n, x, &n, ipiv, &tmp, &lwork, &info);
  if (info != 0) goto cleanup;
  
  lwork = (int) tmp;
  work = malloc(lwork * sizeof(*work));
  dgetri_(&n, x, &n, ipiv, work);
  
  
  free(work);
  cleanup:
  free(ipiv);
  
  return info;
}
