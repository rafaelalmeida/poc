#include "lch.h"

#include "common.h"
#include "image.h"

#define BINS  4
#define SIZE 64

#define NBINS (BINS*BINS*SIZE)

typedef struct _LCHProperty {
  int color;
} LCHProperty;

typedef struct _LCHVisualFeature {
  ulong *colorH;
  int n;
} LCHVisualFeature;

typedef struct _LCHCompressedVisualFeature {
  uchar *colorH;
  int n;
} LCHCompressedVisualFeature;

LCHProperty *LCHAllocPropertyArray(int n)
{
  LCHProperty *v=NULL;
  v = (LCHProperty *) calloc(n,sizeof(LCHProperty));
  if (v==NULL)
    Error(MSG1,"LCHAllocPropertyArray");
  return(v);
}

LCHVisualFeature *CreateVisualFeature(int n)
{
  LCHVisualFeature *vf=NULL;

  vf = (LCHVisualFeature *) calloc(1,sizeof(LCHVisualFeature));
  if (vf != NULL) {
    vf->colorH = AllocULongArray(n);
    vf->n = n;
  } else {
    Error(MSG1,"CreateVisualFeature");
  }
  return(vf);
}

void DestroyVisualFeature(LCHVisualFeature **vf)
{
  LCHVisualFeature *aux;

  aux = *vf;
  if (aux != NULL) {
    if (aux->colorH != NULL) free(aux->colorH);
    free(aux);
    *vf = NULL;
  }
}

LCHCompressedVisualFeature *CreateCompressedVisualFeature(int n)
{
  LCHCompressedVisualFeature *cvf=NULL;

  cvf = (LCHCompressedVisualFeature *) calloc(1,sizeof(LCHCompressedVisualFeature));
  if (cvf != NULL) {
    cvf->colorH = AllocUCharArray(n);
    cvf->n = n;
  } else {
    Error(MSG1,"CreateCompressedVisualFeature");
  }
  return(cvf);
}

void DestroyCompressedVisualFeature(LCHCompressedVisualFeature **cvf)
{
  LCHCompressedVisualFeature *aux;

  aux = *cvf;
  if (aux != NULL) {
    if (aux->colorH != NULL) free(aux->colorH);
    free(aux);
    *cvf = NULL;
  }
}

int *QuantizeColors(CImage *cimg, int color_dim)
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

void CompressHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeNorm((float) h[i] / (float) max);
    ch[i] = (uchar)(v);
  }
}

LCHProperty *ComputePixelsProperties(CImage *cimg)
{
  LCHProperty *p=NULL;
  int *color, i, n;
  
  n = cimg->C[0]->nrows * cimg->C[0]->ncols;  
  
  p = LCHAllocPropertyArray(n);
  
  color = QuantizeColors(cimg, 4);
  for(i=0; i<n; i++) 
    p[i].color=color[i];
  
  free(color);
  return(p);
}

LCHVisualFeature *ComputeHistograms(LCHProperty *p, Image *mask, int *npoints)
{
  LCHVisualFeature *vf=NULL;
  ulong fator_x, fator_y;
  ulong rows, cols;
  ulong x, y;
  ulong r, c;
  ulong i;

  rows = mask->nrows;
  cols = mask->ncols;
  
  vf = CreateVisualFeature(NBINS);
  for(i=0; i<NBINS; i++){
    vf->colorH[i] = 0;
  }
  
  fator_x = SIZE;
  fator_y = BINS * SIZE;

  *npoints = 0;
  for(r=0; r<rows; r++)
    for(c=0; c<cols; c++)
      if (mask->val[mask->tbrow[r]+c]) {
        x = BINS*r/rows;
        y = BINS*c/cols;

        vf->colorH[p[r*cols+c].color + fator_x*x + fator_y*y]++;
        (*npoints)++;
      }

  return(vf);
}

LCHCompressedVisualFeature *CompressHistograms(LCHVisualFeature *vf, int npixels)
{
  LCHCompressedVisualFeature *cvf=NULL;

  cvf = CreateCompressedVisualFeature(NBINS);
  CompressHistogram(cvf->colorH, vf->colorH, npixels, NBINS);
  
  return(cvf);
}

LCHCompressedVisualFeature *ReadCompressedVisualFeatures(char *filename)
{
  LCHCompressedVisualFeature *cvf=NULL;
  FILE *fp;
  int i, n;
  uchar c;

  fp = fopen(filename,"r");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  fscanf(fp,"%d\n",&n);
  cvf = CreateCompressedVisualFeature(n);
  for (i=0; i<n; i++) {
    fscanf(fp,"%c",&c);
    cvf->colorH[i] = c;
  }
  fclose(fp);
  return(cvf);
}

void LCHWriteCompressedVisualFeatures(LCHCompressedVisualFeature *cvf,char *filename)
{
  FILE *fp;
  int i;
  
  fp = fopen(filename,"w");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  
  fprintf(fp,"%d\n",NBINS);
  for(i=0;i<NBINS;i++)
    fprintf(fp, "%c", cvf->colorH[i]);
  
  fclose(fp);
}

LCHCompressedVisualFeature *ExtractCompressedVisualFeaturesLCH(CImage *cimg,
                                                          Image *mask)
{
  LCHCompressedVisualFeature *cvf=NULL;
  LCHVisualFeature *vf;
  LCHProperty *p;
  int npixels;
  int npoints;

  npixels = cimg->C[0]->ncols * cimg->C[0]->nrows;
  
  p = ComputePixelsProperties(cimg);
  vf = ComputeHistograms(p, mask, &npoints);
  cvf = CompressHistograms(vf, npoints);
  
  free(p);
  DestroyVisualFeature(&vf);
  return(cvf);
}

Histogram *LCH(CImage *cimg, Image *mask)
{
  Histogram *histogram = NULL;
  LCHCompressedVisualFeature *cvf;
  int i;
  
  cvf = ExtractCompressedVisualFeaturesLCH(cimg, mask);

  histogram = CreateHistogram(NBINS);  
  for (i=0; i<NBINS; i++)
    histogram->v[i] = cvf->colorH[i];
  DestroyCompressedVisualFeature(&cvf);

  return(histogram);
}

