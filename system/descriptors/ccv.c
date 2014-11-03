#include "ccv.h"

#include "adjacency.h"
#include "common.h"
#include "image.h"
#include "queue.h"

#define MIN_AREA 1

#define LOW  0
#define HIGH 1

#define SIZE 64

VisualFeature *CCVCreateVisualFeature(int n)
{
  VisualFeature *vf=NULL;

  vf = (VisualFeature *) calloc(1,sizeof(VisualFeature));
  if (vf != NULL) {
    vf->lowH = AllocULongArray(n);
    vf->highH = AllocULongArray(n);
    vf->n = n;
  } else {
    Error(MSG1,"CCVCreateVisualFeature");
  }
  return(vf);
}

void CCVDestroyVisualFeature(VisualFeature **vf)
{
  VisualFeature *aux;

  aux = *vf;
  if (aux != NULL) {
    if (aux->lowH != NULL) free(aux->lowH);
    if (aux->highH != NULL) free(aux->highH);
    free(aux);
    *vf = NULL;
  }
}

CompressedVisualFeature *CCVCreateCompressedVisualFeature(int n)
{
  CompressedVisualFeature *cvf=NULL;

  cvf = (CompressedVisualFeature *) calloc(1,sizeof(CompressedVisualFeature));
  if (cvf != NULL) {
    cvf->lowH = AllocUCharArray(n);
    cvf->highH = AllocUCharArray(n);
    cvf->n = n;
  } else {
    Error(MSG1,"CCVCreateCompressedVisualFeature");
  }
  return(cvf);
}

void CCVDestroyCompressedVisualFeature(CompressedVisualFeature **cvf)
{
  CompressedVisualFeature *aux;

  aux = *cvf;
  if (aux != NULL) {
    if (aux->lowH != NULL) free(aux->lowH);
    if (aux->highH != NULL) free(aux->highH);
    free(aux);
    *cvf = NULL;
  }
}

int *CCVQuantizeColors(CImage *cimg, int color_dim)
{
  ulong i;
  ulong r, g, b;
  ulong fator_g, fator_b;
  int *color, n;
  
  n = cimg->C[0]->nrows * cimg->C[0]->ncols;  

  color = AllocIntArray(n);
  
  fator_g = color_dim;
  fator_b = fator_g*color_dim;
  
  for(i=0; i<n; i++){
    r = color_dim*cimg->C[0]->val[i]/256;
    g = color_dim*cimg->C[1]->val[i]/256;
    b = color_dim*cimg->C[2]->val[i]/256;
    
    color[i] = (r + fator_g*g + fator_b*b);
  }
  return(color);
}

void CCVCompressHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeNorm((float) h[i] / (float) max);
    ch[i] = (uchar)(v);
  }
}

void CCVComputeFrequencyProperty(Image *img, Property *ppt, int npixels)
{ 
  int i, x, y, p, q;
  Image *label;
  int nlabels;
  ulong *area;
  Pixel u, v;
  AdjRel *A;
  Queue *Q;
  
  label = CreateImage(img->ncols, img->nrows);
  for (p=0;p<npixels;p++)
    label->val[p] = NIL;

  A = Circular(1.0);

  Q = CreateQueue(SIZE, npixels);
  nlabels = 0;
  for (y=0; y<img->nrows; y++)
    for (x=0; x<img->ncols; x++) {
      p = x + img->tbrow[y];
      if (label->val[p] == NIL) {
        label->val[p] = nlabels++;
        InsertQueue(Q, ppt[p].color, p);
        while (!EmptyQueue(Q)) {
          p = RemoveQueue(Q);
          u.x = p % img->ncols;
          u.y = p / img->ncols;
          for (i=1; i < A->n; i++){
            v.x = u.x + A->dx[i];
            v.y = u.y + A->dy[i];
            if (ValidPixel(img, v.x, v.y)) {
              q = v.x + img->tbrow[v.y];
              if (label->val[q] == NIL && ppt[p].color == ppt[q].color) {
                label->val[q] = label->val[p];
                InsertQueue(Q, ppt[q].color, q);
              }
            }                   
          } 
        }
      }
    }
  DestroyQueue(&Q);
  DestroyAdjRel(&A);

  area = AllocULongArray(nlabels);
  for (i=0;i<nlabels;i++)
    area[i] = 0;

  for (p=0;p<npixels;p++)
    area[label->val[p]]++;

  for (p=0;p<npixels;p++)
    if (100 * area[label->val[p]] < MIN_AREA * npixels)
      ppt[p].frequency = LOW;
    else
      ppt[p].frequency = HIGH;
  DestroyImage(&label);
  free(area);
}

Property *CCVComputePixelsProperties(CImage *cimg)
{
  Property *p=NULL;
  int *color, i, n;
  
  n = cimg->C[0]->nrows * cimg->C[0]->ncols;  
  
  p = AllocPropertyArray(n);
  
  color = CCVQuantizeColors(cimg, 4);
  for(i=0; i<n; i++) 
    p[i].color=color[i];
  CCVComputeFrequencyProperty(cimg->C[0], p, n);
  
  free(color);
  return(p);
}

VisualFeature *CCVComputeHistograms(Property *p, Image *mask, 
                                 int npixels, int *npoints)
{
  VisualFeature *vf=NULL;
  ulong i;
  
  vf = CCVCreateVisualFeature(SIZE);
  for(i=0; i<SIZE; i++){
    vf->lowH[i] = 0;
    vf->highH[i] = 0;
  }
  
  *npoints = 0;
  for(i=0; i<npixels; i++)
    if (mask->val[i]) {
      (*npoints)++;
      if(p[i].frequency==LOW) 
        vf->lowH[p[i].color]++;
      else 
        vf->highH[p[i].color]++;
    }
  return(vf);
}

CompressedVisualFeature *CCVCompressHistograms(VisualFeature *vf, int npixels)
{
  CompressedVisualFeature *cvf=NULL;
  
  cvf = CCVCreateCompressedVisualFeature(SIZE);
  CCVCompressHistogram(cvf->lowH, vf->lowH, npixels, SIZE);
  CCVCompressHistogram(cvf->highH, vf->highH, npixels, SIZE);
  
  return(cvf);
}

CompressedVisualFeature *CCVReadCompressedVisualFeatures(char *filename)
{
  CompressedVisualFeature *cvf=NULL;
  FILE *fp;
  int i, n;
  uchar l, h;

  fp = fopen(filename,"r");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  fscanf(fp,"%d\n",&n);
  cvf = CCVCreateCompressedVisualFeature(n);
  for (i=0; i<n; i++) {
    fscanf(fp,"%c%c",&l,&h);
    cvf->lowH[i] = l;
    cvf->highH[i] = h;
  }
  fclose(fp);
  return(cvf);
}

void CCVWriteCompressedVisualFeatures(CompressedVisualFeature *cvf,char *filename)
{
  FILE *fp;
  int i;
  
  fp = fopen(filename,"w");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  
  fprintf(fp,"%d\n",SIZE);
  for(i=0;i<SIZE;i++)
    fprintf(fp, "%c%c", cvf->lowH[i], cvf->highH[i]);
  
  fclose(fp);
}

CompressedVisualFeature *CCVExtractCompressedVisualFeatures(CImage *cimg, 
                                                          Image *mask)
{
  CompressedVisualFeature *cvf=NULL;
  VisualFeature *vf;
  Property *p;
  int npixels;
  int npoints;
  
  npixels = cimg->C[0]->nrows * cimg->C[0]->ncols;
  
  p = CCVComputePixelsProperties(cimg);
  vf = CCVComputeHistograms(p, mask, npixels, &npoints);
  cvf = CCVCompressHistograms(vf, npoints);
  
  free(p);
  CCVDestroyVisualFeature(&vf);
  return(cvf);
}

Histogram *CCV(CImage *cimg, Image *mask)
{
  Histogram *histogram = NULL;
  CompressedVisualFeature *cvf;
  int i;
  
  cvf = CCVExtractCompressedVisualFeatures(cimg, mask);

  histogram = CreateHistogram(2*SIZE);  
  for (i=0; i<SIZE; i++){
    histogram->v[i] = cvf->lowH[i];
    histogram->v[i+SIZE] = cvf->highH[i];
  }
  CCVDestroyCompressedVisualFeature(&cvf);

  return(histogram);
}

