#include "acc.h"

#include "adjacency.h"
#include "common.h"
#include "image.h"

#define CDIM 4
#define SIZE (4 * CDIM * CDIM * CDIM)

#define BINS  4

#define NBINS (BINS*BINS*SIZE)

typedef struct _ACCProperty {
  int color;
  int frequency[4];
} ACCProperty;

typedef struct _ACCVisualFeature {
  ulong *colorH;
  int n;
} ACCVisualFeature;

typedef struct _ACCCompressedVisualFeature {
  uchar *colorH;
  int n;
} ACCCompressedVisualFeature;

ACCCompressedVisualFeature *ACCCreateCompressedVisualFeature(int n);

void ACCCompressHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeNorm((float) h[i] / (float) max);
    ch[i] = (uchar)(v);
  }
}

ACCProperty *ACCAllocPropertyArray(int n)
{
  ACCProperty *v=NULL;
  v = (ACCProperty *) calloc(n,sizeof(ACCProperty));
  if (v==NULL)
    Error(MSG1,"ACCAllocPropertyArray");
  return(v);
}

ACCCompressedVisualFeature *ACCCompressHistograms(ACCVisualFeature *vf, int npixels)
{
  ACCCompressedVisualFeature *cvf=NULL;

  cvf = ACCCreateCompressedVisualFeature(NBINS);
  ACCCompressHistogram(cvf->colorH, vf->colorH, npixels, NBINS);
  
  return(cvf);
}

ACCVisualFeature *ACCCreateVisualFeature(int n)
{
  ACCVisualFeature *vf=NULL;

  vf = (ACCVisualFeature *) calloc(1,sizeof(ACCVisualFeature));
  if (vf != NULL) {
    vf->colorH = AllocULongArray(n);
    vf->n = n;
  } else {
    Error(MSG1,"ACCCreateVisualFeature");
  }
  return(vf);
}

void ACCDestroyVisualFeature(ACCVisualFeature **vf)
{
  ACCVisualFeature *aux;

  aux = *vf;
  if (aux != NULL) {
    if (aux->colorH != NULL) free(aux->colorH);
    free(aux);
    *vf = NULL;
  }
}

ACCCompressedVisualFeature *ACCCreateCompressedVisualFeature(int n)
{
  ACCCompressedVisualFeature *cvf=NULL;

  cvf = (ACCCompressedVisualFeature *) calloc(1,sizeof(ACCCompressedVisualFeature));
  if (cvf != NULL) {
    cvf->colorH = AllocUCharArray(n);
    cvf->n = n;
  } else {
    Error(MSG1,"ACCCreateCompressedVisualFeature");
  }
  return(cvf);
}

void ACCDestroyCompressedVisualFeature(ACCCompressedVisualFeature **cvf)
{
  ACCCompressedVisualFeature *aux;

  aux = *cvf;
  if (aux != NULL) {
    if (aux->colorH != NULL) free(aux->colorH);
    free(aux);
    *cvf = NULL;
  }
}

int *ACCQuantizeColors(CImage *cimg, int color_dim)
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

uchar ComputeNormACC(float value)
{
  return ((uchar)(255. * value));
}

void LinearNormalizeHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeNormACC((float) h[i] / (float) max);
    ch[i] = (uchar)(v);
  }
}

void NonLinearNormalizeHistogram(uchar *ch, ulong *h, ulong max, int size)
{
  int i;
  uchar v;
  
  for(i=0; i<size; i++){
    v = ComputeLog((float) h[i] / (float) max);
    ch[i] = (uchar)('0' + v);
  }
}

void ComputeFrequencyProperty(Image *img, ACCProperty *ppt)
{ 
  ulong x, y, p, q;
  uchar d, r;
  AdjRel *A;
  Pixel v;
  int i;

  for (p=0; p<img->nrows*img->ncols; p++)
    for (d=0; d<4; d++)
      ppt[p].frequency[d] = 0;

  A = Circular(1.0);
  for(y=0; y<img->nrows; y++)
    for(x=0; x<img->ncols; x++){
      p = x + img->tbrow[y];
      for(r=1,d=0; r<=7; r+=2,d++)
        for (i=1; i < A->n; i++){
          v.x = x + r * A->dx[i]; 
          v.y = y + r * A->dy[i];
          if (ValidPixel(img,v.x,v.y)){
            q = v.x + img->tbrow[v.y];
            if(ppt[p].color == ppt[q].color)
	      ppt[p].frequency[d]++;
	  }
        }
    }
  DestroyAdjRel(&A);
}

ACCProperty *ACCComputePixelsProperties(CImage *cimg)
{
  ACCProperty *p=NULL;
  int *color, i, n;
  
  n = cimg->C[0]->nrows * cimg->C[0]->ncols;  
  
  p = ACCAllocPropertyArray(n);
  
  color = ACCQuantizeColors(cimg, CDIM);
  for(i=0; i<n; i++)
    p[i].color=color[i];
  ComputeFrequencyProperty(cimg->C[0], p);
  
  free(color);
  return(p);
}

ACCVisualFeature *ComputeHistogramsACC(ACCProperty *p, Image *mask, 
                                 int npixels, int *npoints)
{
  ACCVisualFeature *vf=NULL;
  ulong i, d;
  
  vf = ACCCreateVisualFeature(SIZE);
  for(i=0; i<SIZE; i++)
    vf->colorH[i] = 0;

  *npoints = 0;
  for(d=0; d<4; d++)
    for(i=0; i<npixels; i++) {
             if (mask->val[i]) {
                vf->colorH[4 * p[i].color + d] += p[i].frequency[d];
                (*npoints)++;
             }
    }
  return(vf);
}

ACCCompressedVisualFeature *NormalizeHistograms(ACCVisualFeature *vf, int npixels)
{
  ACCCompressedVisualFeature *cvf=NULL;
  
  cvf = ACCCreateCompressedVisualFeature(SIZE);
  NonLinearNormalizeHistogram(cvf->colorH, vf->colorH, 4 * npixels, SIZE);
  
  return(cvf);
}

ACCCompressedVisualFeature *ACCReadCompressedVisualFeatures(char *filename)
{
  ACCCompressedVisualFeature *cvf=NULL;
  FILE *fp;
  int i, n;
  uchar c;

  fp = fopen(filename,"r");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }
  fscanf(fp,"%d\n",&n);
  cvf = ACCCreateCompressedVisualFeature(n);
  for (i=0; i<n; i++) {
    fscanf(fp,"%c\n",&c);
    cvf->colorH[i] = c;
  }
  fclose(fp);
  return(cvf);
}

void ACCWriteCompressedVisualFeatures(ACCCompressedVisualFeature *cvf,char *filename)
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

ACCCompressedVisualFeature *ExtractCompressedVisualFeaturesACC(CImage *cimg, 
                                                          Image *mask)
{
  ACCCompressedVisualFeature *cvf=NULL;
  ACCVisualFeature *vf;
  ACCProperty *p;
  int npixels;
  int npoints;
  
  npixels = cimg->C[0]->nrows * cimg->C[0]->ncols;
  
  p = ACCComputePixelsProperties(cimg);
  vf = ComputeHistogramsACC(p, mask, npixels, &npoints);
  cvf = ACCCompressHistograms(vf, npoints);
  
  free(p);
  ACCDestroyVisualFeature(&vf);
  return(cvf);
}

Histogram *ACC(CImage *cimg, Image *mask)
{
  Histogram *h = NULL;
  ACCCompressedVisualFeature *cvf;
  int i;

  cvf = ExtractCompressedVisualFeaturesACC(cimg, mask);

  h = CreateHistogram(SIZE);  
  for (i=0; i<SIZE; i++)
    h->v[i] = cvf->colorH[i];
  ACCDestroyCompressedVisualFeature(&cvf);

  return(h);
}

int ACCDimensions() {
  return SIZE;
}
