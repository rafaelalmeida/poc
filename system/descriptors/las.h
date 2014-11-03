#ifndef _LAS_H_
#define _LAS_H_

#include "image.h"
#include "cimage.h"
#include "common.h"

#define QTD_BINS 4  //qtd de bins para cada gi

typedef struct _imageHSV {
  float H,S,V;
} ImageHSV;

/*Funcoes auxiliares*/
void RGB2HSV_las(CImage *RGB, ImageHSV **HSV);
void CalculaGradientes(Image *img, Image *g1, Image *g2, Image *g3, Image *g4);
float *HistogramaGradiente(Image *g1, Image *g2, Image *g3, Image *g4);
void NormalizaHistograma(float *hist, int hist_size, int qtd_pixels);

#endif
