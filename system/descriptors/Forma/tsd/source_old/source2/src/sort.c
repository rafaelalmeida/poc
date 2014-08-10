#include "sort.h"

void SelectionSort(int *val, int nelems, char order)
{
  int i,j,jmax;

  if (order == INCREASING){
    for (i=nelems-1; i > 0; i--){ 
      jmax = i;
      for (j=0; j < i; j++){
	if (val[j] > val[jmax])
	  jmax = j;
      }
      if (jmax != i)
	Change(&val[i],&val[jmax]);
    }
  } else { /* DECREASING */
    for (i=0; i < nelems-1; i++){ 
      jmax = i;
      for (j=i+1; j < nelems; j++){
	if (val[j] > val[jmax])
	  jmax = j;
      }
      if (jmax != i)
	Change(&val[i],&val[jmax]);
    }
  }
}

int *BucketSort(int *val, int nelems, char order)
{
  int i,j,maxval=INT_MIN,*sval=NULL;
  Queue *Q=NULL;

  for(i=0; i < nelems; i++) 
    if (val[i] > maxval)
      maxval = val[i];
  
  Q = CreateQueue(maxval+1,nelems);
  for(i=0; i < nelems; i++)
    InsertQueue(Q,val[i],i);

  sval = AllocIntArray(nelems);
  if (order==INCREASING){
    j = 0;
    while(!EmptyQueue(Q)) {
      i = RemoveQueue(Q);
      sval[j] = val[i];
      j++;
    }
  } else { /* order = DECREASING */
    j = nelems-1;
    while(!EmptyQueue(Q)) {
      i = RemoveQueue(Q);
      sval[j] = val[i];
      j--;
    }  
  }

  DestroyQueue(&Q);
  
  return(sval);
}

