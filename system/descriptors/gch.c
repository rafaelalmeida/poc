#include "gch.h"

#include "common.h"
#include "image.h"

#define SIZE 64

typedef struct _GCHProperty {
  int color;
} GCHProperty;

typedef struct _GCHVisualFeature {
  ulong *colorH;
  int n;
} GCHVisualFeature;

typedef struct _GCHCompressedVisualFeature {
  uchar *colorH;
  int n;
} GCHCompressedVisualFeature;

GCHProperty *GCHAllocPropertyArray(int n)
{
  GCHProperty *v=NULL;
  v = (GCHProperty *) calloc(n,sizeof(GCHProperty));
  if (v==NULL)
    Error(MSG1,"GCHAllocPropertyArray");
  return(v);
}

GCHVisualFeature *GCHCreateVisualFeature(int n)
{
  GCHVisualFeature *vf=NULL;

  vf = (GCHVisualFeature *) calloc(1,sizeof(GCHVisualFeature));
  if (vf != NULL) {
    vf->colorH = AllocULongArray(n);
    vf->n = n;
  } else {
    Error(MSG1,"GCHCreateVisualFeature");
  }
  return(vf);
}

void GCHDestroyVisualFeature(GCHVisualFeature **vf)
{
  GCHVisualFeature *aux;

  aux = *vf;
  if (aux != NULL) {
    if (aux->colorH != NULL) free(aux->colorH);
    free(aux);
    *vf = NULL;
  }
}

GCHCompressedVisualFeature *GCHCreateCompressedVisualFeature(int n)
{
  GCHCompressedVisualFeature *cvf=NULL;

  cvf = (GCHCompressedVisualFeature *) calloc(1,sizeof(GCHCompressedVisualFeature));
  if (cvf != NULL) {
    cvf->colorH = AllocUCharArray(n);
    cvf->n = n;
  } else {
    Error(MSG1,"GCHCreateCompressedVisualFeature");
  }
  return(cvf);
}

void GCHDestroyCompressedVisualFeature(GCHCompressedVisualFeature **cvf)
{
  GCHCompressedVisualFeature *aux;

  aux = *cvf;
  if (aux != NULL) {
    if (aux->colorH != NULL) free(aux->colorH);
    free(aux);
    *cvf = NULL;
  }
}

int *GCHQuantizeColors(CImage *cimg, int color_dim)
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

void GCHCompressHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeNorm((float) h[i] / (float) max);
    ch[i] = (uchar)(v);
  }
}

GCHProperty *GCHComputePixelsProperties(CImage *cimg)
{
  GCHProperty *p=NULL;
  int *color, i, n;
  
  n = cimg->C[0]->nrows * cimg->C[0]->ncols;  
  
  p = GCHAllocPropertyArray(n);
  
  color = GCHQuantizeColors(cimg, 4);
  for(i=0; i<n; i++) 
    p[i].color=color[i];
  
  free(color);
  return(p);
}

GCHVisualFeature *GCHComputeHistograms(GCHProperty *p, Image *mask, 
                                 int npixels, int *npoints)
{
  GCHVisualFeature *vf=NULL;
  ulong i;
  
  vf = GCHCreateVisualFeature(SIZE);
  for(i=0; i<SIZE; i++){
    vf->colorH[i] = 0;
  }
  
  *npoints = 0;
  for(i=0; i<npixels; i++)
    if (mask->val[i]) {
      vf->colorH[p[i].color]++;
      (*npoints)++;
    }
  return(vf);
}

GCHCompressedVisualFeature *GCHCompressHistograms(GCHVisualFeature *vf, int npixels)
{
  GCHCompressedVisualFeature *cvf=NULL;
  
  cvf = GCHCreateCompressedVisualFeature(SIZE);
  GCHCompressHistogram(cvf->colorH, vf->colorH, npixels, SIZE);
  
  return(cvf);
}

GCHCompressedVisualFeature *GCHReadCompressedVisualFeatures(char *filename)
{
  GCHCompressedVisualFeature *cvf=NULL;
  FILE *fp;
  int i, n;
  uchar c;

  fp = fopen(filename,"r");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  fscanf(fp,"%d\n",&n);
  cvf = GCHCreateCompressedVisualFeature(n);
  for (i=0; i<n; i++) {
    fscanf(fp,"%c\n",&c);
    cvf->colorH[i] = c;
  }
  fclose(fp);
  return(cvf);
}

void WriteCompressedVisualFeatures(GCHCompressedVisualFeature *cvf,char *filename)
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
    fprintf(fp, "%c", cvf->colorH[i]);
  
  fclose(fp);
}

GCHCompressedVisualFeature *ExtractCompressedVisualFeatures(CImage *cimg, 
                                                          Image *mask)
{
  GCHCompressedVisualFeature *cvf=NULL;
  GCHVisualFeature *vf;
  GCHProperty *p;
  int npixels;
  int npoints;
  
  npixels = cimg->C[0]->nrows * cimg->C[0]->ncols;
  
  p = GCHComputePixelsProperties(cimg);
  vf = GCHComputeHistograms(p, mask, npixels, &npoints);
  cvf = GCHCompressHistograms(vf, npoints);
  
  free(p);
  GCHDestroyVisualFeature(&vf);
  return(cvf);
}

Histogram *GCH(CImage *cimg, Image *mask)
{
  Histogram *histogram = NULL;
  GCHCompressedVisualFeature *cvf;
  int i;
  
  cvf = ExtractCompressedVisualFeatures(cimg, mask);

  histogram = CreateHistogram(SIZE);  
  for (i=0; i<SIZE; i++)
    histogram->v[i] = cvf->colorH[i];
  GCHDestroyCompressedVisualFeature(&cvf);

  return(histogram);
}

