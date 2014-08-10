#include "libtexturedescriptors.h"

void Extraction(char *img_path, char *fv_path)
{

/* ----------------------------------------------------------------------
      Declarações e início para verificação do uso da memória dinâmica.
      ---------------------------------------------------------------------- *
   struct mallinfo info = mallinfo();
   int MemDinInicial = info.uordblks, MemDinFinal;
   // ---------------------------------------------------------------------- */

  CImage *cimg=NULL;
  Features *img=NULL;
  Features *steerable=NULL;
  FeatureVector1D *featurevector=NULL;
  int i,j;

  cimg = ReadCImage(img_path);
  img = CreateFeatures(cimg->C[0]->ncols, cimg->C[0]->nrows, 1);
  
  //translate to HSV space

  for (i=0; i<(img->ncols*img->nrows); i++)
    img->elem[i].feat[0] = (cimg->C[0]->val[i])+(cimg->C[1]->val[i])+(cimg->C[2]->val[i])/3;

  steerable = SteerablePyramidFeats(img);
  
  featurevector = CreateFeatureVector1D(2*steerable->nfeats);
  
  //mean of the steerable features
  for (i=0; i<steerable->nelems; i++)
  for (j=0; j<steerable->nfeats; j++)
  		featurevector->X[j] += (double) (img->elem[i].feat[j] / steerable->nelems);
  		
  //standad deviation of the steerable features
  for (i=0; i<steerable->nelems; i++)
  for (j=0; j<steerable->nfeats; j++){
  		featurevector->X[j+steerable->nfeats] += (double) 
  		((img->elem[i].feat[j]-featurevector->X[j])*(img->elem[i].feat[j]-featurevector->X[j]) / steerable->nelems);
  	}
  for (j=0; j<steerable->nfeats; j++)
  		featurevector->X[j+steerable->nfeats] = sqrt(featurevector->X[j+steerable->nfeats]);

  //
  // Sort feature vector to achive rotation and scale invariance
  // 
/*
  int maxpos;
  double max, aux;

  for (i=0; i<steerable->nfeats-1; i++){
	max = featurevector->X[i];
	maxpos = i;
     	for (j=i+1; j<steerable->nfeats; j++){
		if (featurevector->X[j]>max){
			max = featurevector->X[j];
			maxpos = j;
		}
	}
	aux = featurevector->X[maxpos+steerable->nfeats];
	featurevector->X[maxpos] = featurevector->X[i];
	featurevector->X[maxpos+steerable->nfeats] = featurevector->X[i+steerable->nfeats];
	featurevector->X[i] = max;
	featurevector->X[i+steerable->nfeats] = aux;
  }

*/
  //
  // End sort. Comment lines above if rotation or scale invariance isn't desired.
  //

  WriteFileFeatureVector1D(featurevector,fv_path);
  
  DestroyFeatureVector1D(&featurevector);
  DestroyCImage(&cimg);
  DestroyFeatures(&img);
  DestroyFeatures(&steerable);

/* --------------------------------------------------------------------- */
   /* Verificação final do uso da memória dinâmica                          */
   /* --------------------------------------------------------------------- *
   info = mallinfo();
   MemDinFinal = info.uordblks;
   if (MemDinInicial!=MemDinFinal)
      printf("\n\nMemória dinâmica não foi totalmente liberada (%d, %d)\n",
	     MemDinInicial,MemDinFinal);
   // --------------------------------------------------------------------- */

}


int main(int argc, char** argv)
{

  if (argc != 3) {
    fprintf(stderr,"usage: gch_extraction <image_path> <fv_path>\n");
    exit(-1);
  }

  Extraction(argv[1], argv[2]);

  return(0);
}
