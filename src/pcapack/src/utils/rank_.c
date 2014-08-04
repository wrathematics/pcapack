// Rank statistics
// Copyright 2013, Heckendorf and Schmidt

#include <stdlib.h>
#include "rand.h"

#define RANK_INDEX(i) (ptr[i]-data)


// tiefighters
static double tb_avg(const double **ptr, const double *data, const double *rank, const int index, const int len){
  int i;
  double ret=0;
  const int end=len+index;
  
  for(i=index; i<end; i++)
    ret += rank[RANK_INDEX(i)]/(double)len;
  
  return ret;
}

static double tb_min(const double **ptr, const double *data, const double *rank, const int index, const int len){
  int i;
  double ret=rank[RANK_INDEX(index)];
  const int end=len+index;
  
  for(i=index; i<end; i++){
    if (rank[RANK_INDEX(i)] < ret)
      ret = rank[RANK_INDEX(i)];
  }
  
  return ret;
}

static double tb_max(const double **ptr, const double *data, const double *rank, const int index, const int len){
  int i;
  double ret=rank[RANK_INDEX(index)];
  const int end=len+index;
  
  for(i=index; i<end; i++){
    if (rank[RANK_INDEX(i)] > ret)
      ret = rank[RANK_INDEX(i)];
  }
  
  return ret;
}


// Rank functions
static int runlength(const double **ptr, const int len){
    int ret=1;
    int i;
        
    for(i=1;i<len && *ptr[i]==*ptr[0];i++)
        ret++;

    return ret;
}

static int cmp_data(const void *a, const void *b){
    const double da = **(double**)a;
    const double db = **(double**)b;

    if(da<db) return -1; 
    else if(da>db) return 1;
    else return 0;
}


static void rank1(const double *data, double *rank, const int len)
{
  int i;
  int step;
  const double **ptr=malloc(sizeof(double*)*len);
  
  
  for(i=0;i<len;i++)
    ptr[i]=data+i;
  
  qsort(ptr,len,sizeof(*ptr),cmp_data);
  
  for(i=0;i<len;i++)
    rank[(ptr[i]-data)]=i+1;
  
  for(i=0;i<len;i++)
    step=runlength(data+i,len-i);
  
  free(ptr);
}

// *data = input array to be ranked
// *rank = output array, allocated by user
// len = length of data
// (other crap) = callback for tiebreaker, e.g. tb_avg
static void rank2(const double *data, double *rank, const int len, double (*tiefighter)(const double **ptr, const double *data, const double *rank, const int index, const int len)){
  int i;
  int run;
  int step;
  const double **ptr=malloc(sizeof(double*)*len);
  
  for(i=0;i<len;i++)
    ptr[i]=data+i;
    
  qsort(ptr,len,sizeof(*ptr),cmp_data);
    
  for(i=0;i<len;i++)
    rank[RANK_INDEX(i)]=i+1;
    
  for(i=0;i<len;i++){
    step=run=runlength(ptr+i,len-i);
    if(run>1){
      rank[RANK_INDEX(i)]=tiefighter(ptr,data,rank,i,run);
      while(run>1){
        run--;
        rank[RANK_INDEX(i+run)]=rank[RANK_INDEX(i)];
      }
      i+=step-1;
    }
  }
  
  free(ptr);
}


void genrank(const double *data, double *rank, const int len, int method)
{
  double (*cb)(const double*,const int);
  switch(method){
    case 1:
      rank2(data, rank, len, tb_min);
      break;
    case 2:
      rank2(data, rank, len, tb_max);
      break;
    case 3:
      rank2(data, rank, len, tb_avg);
      break;
    default: 
      rank1(data, rank, len);
  }
}

