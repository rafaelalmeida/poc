#ifndef TENSORSCALE_H
#define TENSORSCALE_H 1

#include "common.h"
#include "image.h"
#include "cimage.h"
#include "dimage.h"
#include "geometry.h"
#include "color.h"
#include "adjacency.h"
#include "queue.h"
#include "segmentation.h"
#include "descriptor.h"

#define HISTOGRAMSIZE 180
#define DEBUG 0

/*
typedef struct _tensorscale {
  DImage *anisotropy;
  DImage *orientation;
  DImage *thickness;

  int m_pairs;
} TensorScale;
*/
typedef struct _tensorscale{
  DImage *anisotropy;
  DImage *orientation;
  DImage *thickness;

  int m_pairs;
  //int L;
}TensorScale;

typedef struct _imageHSV {
  float H,S,V;
} ImageHSV;

/* Tensor Scale */
TensorScale *CreateBinaryTensorScale(Image *bin, int m_pairs);
void        DestroyTensorScale(TensorScale **ts);
CImage      *ConvertTS2CImage(TensorScale *ts);
void        OutputTSColorSpace(char *filename, int size);

/* Tensor Scale Descriptor (TSD) */
float  *TSOrientationHistogram(TensorScale *ts);
CImage *TSShowHistograms(float *hist1, float *hist2, int offset);
float  TSHistogramMatch(float *hist1, float *hist2, int *offset);

/* Tensor Scale Descriptor with Influence Zones (TSDIZ) */
FeatureVector1D *TSDIZ_ExtractionAlgorithm(Image *bin, int nsamples);
void            WriteTSDIZ(FeatureVector1D *desc,char *filename);
double          TSDIZ_SimilarityAlgorithm(FeatureVector1D *desc1, FeatureVector1D *desc2);

/* Tensor Scale Contour Salience (TSCS) */
FeatureVector2D *TSCS_ExtractionAlgorithm(Image *bin, double threshold); //Salience detector

/* Tensor Scale Sibgrapi 2005 - NAO TENHO CERTEZA - COPIEI APENAS ESTA FUNCAO DA IFT DA FERNANDA */
TensorScale *CreateTensorScale(Image *img, Image *mask, int m_pairs, int L, float stdDev);
/*CONFERIR AR FUNCOES:
ok- void DestroyTensorScale(TensorScale **ts);
D- CImage *ConvertTS2CImage(TensorScale *ts);
     -> possui diferencas no ratio
ok- void OutputTSColorSpace(char *filename, int size);
D- float *TSOrientationHistogram(TensorScale *ts);
     -> histograma eh float
     -> if para verificar se incrementa ou nao eh diferente
     -> incrementa em "an" o bin do histograma (e nao em 1)
     -> normaliza o histograma
     -> alguns tipos de variaveis estao diferentes
D- CImage *TSShowHistograms(int *hist1, int *hist2, int offset);
     -> muitas diferencas, mas esta funcao parece ser usada apenas para exibir os histogramas
- float TSHistogramMatch(int *hist1, int *hist2, int *offset);
     -> ultima posicao de newh1 é dividida por 4 (e nao por 3*n1)
     -> nao possui variaveis n1, n2

*/

void RGB2HSV_tsd(CImage *RGB, ImageHSV **HSV);
void WriteTSD(float *hist, char *filename);
float *ReadTSD(char *filename);

#endif
