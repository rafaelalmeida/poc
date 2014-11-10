#include "bic.h"

#include "adjacency.h"
#include "common.h"
#include "image.h"

#define LOW  0
#define HIGH 1

#define SIZE 64

VisualFeature *BICCreateVisualFeature(int n)
{
  VisualFeature *vf=NULL;

  vf = (VisualFeature *) calloc(1,sizeof(VisualFeature));
  if (vf != NULL) {
    vf->lowH = AllocULongArray(n);
    vf->highH = AllocULongArray(n);
    vf->n = n;
  } else {
    Error(MSG1,"BICCreateVisualFeature");
  }
  return(vf);
}

void BICDestroyVisualFeature(VisualFeature **vf)
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

CompressedVisualFeature *BICCreateCompressedVisualFeature(int n)
{
  CompressedVisualFeature *cvf=NULL;

  cvf = (CompressedVisualFeature *) calloc(1,sizeof(CompressedVisualFeature));
  if (cvf != NULL) {
    cvf->lowH = AllocUCharArray(n);
    cvf->highH = AllocUCharArray(n);
    cvf->n = n;
  } else {
    Error(MSG1,"BICCreateCompressedVisualFeature");
  }
  return(cvf);
}

void BICDestroyCompressedVisualFeature(CompressedVisualFeature **cvf)
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

int *BICQuantizeColors(CImage *cimg, int color_dim)
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

void BICCompressHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeLog((float) h[i] / (float) max);
    ch[i] = (uchar)(48 + v);
  }
}

void BICComputeFrequencyProperty(Image *img, Property *ppt)
{ 
  ulong x, y, p, q;
  int i, border;
  AdjRel *A;
  Pixel v;
  
  A = Circular(1.0);
  
  for(y=0; y<img->nrows; y++){
    for(x=0; x<img->ncols; x++){
      p = x + img->tbrow[y];
      border=false;
      for (i=1; i < A->n; i++){
	v.x = x + A->dx[i];
	v.y = y + A->dy[i];
	if (ValidPixel(img,v.x,v.y)){
	  q = v.x + img->tbrow[v.y];
	  if(ppt[p].color!=ppt[q].color){ 
	    border=true;
	    break;
	  }
	}
      }
      if(border==false) 
	ppt[p].frequency=LOW;
      else 
        ppt[p].frequency=HIGH;
    }
  }
  DestroyAdjRel(&A);
}

Property *BICComputePixelsProperties(CImage *cimg)
{
  Property *p=NULL;
  int *color, i, n;
  
  n = cimg->C[0]->nrows * cimg->C[0]->ncols;  
  
  p = AllocPropertyArray(n);
  
  color = BICQuantizeColors(cimg, 4);
  for(i=0; i<n; i++) 
    p[i].color=color[i];
  BICComputeFrequencyProperty(cimg->C[0], p);
  
  free(color);
  return(p);
}

VisualFeature *BICComputeHistograms(Property *p, Image *mask, 
                                 int npixels, int *npoints)
{
  VisualFeature *vf=NULL;
  ulong i;
  
  vf = BICCreateVisualFeature(SIZE);
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

CompressedVisualFeature *BICCompressHistograms(VisualFeature *vf, int npixels)
{
  CompressedVisualFeature *cvf=NULL;
  
  cvf = BICCreateCompressedVisualFeature(SIZE);
  BICCompressHistogram(cvf->lowH, vf->lowH, npixels, SIZE);
  BICCompressHistogram(cvf->highH, vf->highH, npixels, SIZE);
  
  return(cvf);
}

CompressedVisualFeature *BICReadCompressedVisualFeatures(char *filename)
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
  cvf = BICCreateCompressedVisualFeature(n);
  for (i=0; i<n; i++) {
    fscanf(fp,"%c%c",&l,&h);
    cvf->lowH[i] = l;
    cvf->highH[i] = h;
  }
  fclose(fp);
  return(cvf);
}

void BICWriteCompressedVisualFeatures(CompressedVisualFeature *cvf,char *filename)
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

CompressedVisualFeature *ExtractCompressedVisualFeaturesBIC(CImage *cimg, 
                                                          Image *mask)
{
  CompressedVisualFeature *cvf=NULL;
  VisualFeature *vf;
  Property *p;
  int npixels;
  int npoints;
  
  npixels = cimg->C[0]->nrows * cimg->C[0]->ncols;
  
  p = BICComputePixelsProperties(cimg);
  vf = BICComputeHistograms(p, mask, npixels, &npoints);
  cvf = BICCompressHistograms(vf, npoints);
  
  free(p);
  BICDestroyVisualFeature(&vf);
  return(cvf);
}

Histogram *BIC(CImage *cimg, Image *mask)
{
  Histogram *histogram = NULL;
  CompressedVisualFeature *cvf;
  int i;
  
  cvf = ExtractCompressedVisualFeaturesBIC(cimg, mask);

  histogram = CreateHistogram(2*SIZE);  
  for (i=0; i<SIZE; i++){
    histogram->v[i] = cvf->lowH[i];
    histogram->v[i+SIZE] = cvf->highH[i];
  }
  BICDestroyCompressedVisualFeature(&cvf);

  return(histogram);
}

int BICDimensions() {
  return 2*SIZE;
}
