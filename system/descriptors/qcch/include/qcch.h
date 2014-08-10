#ifndef _QCCH_H_
#define _QCCH_H_

#include "image.h"
#include "cimage.h"
#include "common.h"

#define RAIO     1
#define QTD_VIZ  (1+(4*RAIO)+(4*RAIO*RAIO))  //esse valor depende do raio=a formula funcionou para RAIO ate 4 (nao testei para valores maiores)
#define QTD_BINS 40  //qtde de bins de t(i,j) no espaco quantizado nao uniformemente

Image *RGB2Cinza(CImage *image);
long int SomaJanela(int centroI, int centroJ, int raio, Image *img);

void Rate_Change_Gray(Image *image, Image *t, int r);
void Change_Histogram(Image *t, double h[], int r);
double *QCCH(CImage *cimg, int raio);

double Dist_Measure(double Ha[],double Hb[]);

void ReadFile(char *filename, double h[]);
double *ReadFileBin(char *filename, int size);

void WriteStream(double *h,FILE *fp);
void WriteFile(double h[],char *filename);
void WriteFileBin(double *h, char *filename, int size);

#endif
