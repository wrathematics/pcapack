#include <stdint.h>
#include <string.h>
#include <stdbool.h>


void r_set_na_real_(double *val)
{
  int64_t x = 0x7FF00000000007A2LL;
  memcpy((void *) val, (void *) &x, 8);
}

void r_set_nan_real_(double *val)
{
  int64_t x = 0x7FF0000000000000LL;
  memcpy((void *) val, (void *) &x, 8);
}


void r_set_na_int_(int *val)
{
  *val = INT32_MIN;
}



bool pcapack_anyna(const int n, const double *x)
{
  int i;
  const double NA;
  r_set_na_real(&NA);
  
  for (i=0; i<n; i++)
  {
    if (x[i] == NA)
      return true;
  }
  
  return false;
}

