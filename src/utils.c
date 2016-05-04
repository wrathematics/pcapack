#include <math.h>
#include "pcapack.h"


#define PREBUFLEN 2

SEXP make_pca_default_colnames(const int n)
{
  R_INIT;
  int i;
  int buflen;
  SEXP ret;
  
  buflen = (int) (ceil(log10((double)n)) + PREBUFLEN);
  char *buf = malloc(buflen * sizeof(*buf));
  buf[0] = 'P';
  buf[1] = 'C';
  
  newRlist(ret, n);
  
  for (i=0; i<n; i++)
  {
    sprintf(buf+PREBUFLEN, "%d", i+1);
    buflen = (int) (ceil(log10((double)i+PREBUFLEN)) + PREBUFLEN);
    buflen = MAX(buflen, PREBUFLEN+1);
    SET_VECTOR_ELT(ret, i, mkCharLen(buf, buflen));
  }
  
  free(buf);
  
  R_END;
  return ret;
}

