#ifndef _REALHEAP_H_
#define _REALHEAP_H_

#include "common.h"
#include "heap.h"
#include "gqueue.h"


typedef struct _realheap {
  real *cost;
  char *color;
  int *pixel;
  int *pos;
  int last;
  int n;
  char removal_policy; /* 0 is MINVALUE and 1 is MAXVALUE */
} RealHeap;


void SetRemovalPolicyRealHeap(RealHeap *H, char policy);
bool IsFullRealHeap(RealHeap *H);
bool IsEmptyRealHeap(RealHeap *H);
RealHeap *CreateRealHeap(int n, real *cost);
void DestroyRealHeap(RealHeap **H);
bool InsertRealHeap(RealHeap *H, int pixel);
bool RemoveRealHeap(RealHeap *H, int *pixel);
void UpdateRealHeap(RealHeap *H, int p, real value);
void GoUpRealHeap(RealHeap *H, int i);
void GoDownRealHeap(RealHeap *H, int i);
void ResetRealHeap(RealHeap *H);

#endif



