#ifndef _SORT_H_
#define _SORT_H_

#include "queue.h"
#include "realheap.h"

int   *BucketSort(int *val, int nelems, char order);
void   SelectionSort(int *val, int nelems, char order); /* in place */

#endif
