#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Ap_FeatureVector1D *cbc1=NULL;
  Ap_FeatureVector1D *cbc2=NULL;
  FeatureVector1D *wcbc1=NULL;
  FeatureVector1D *wcbc2=NULL;
  double distance;
  int n1, n2;
  int i;

  if (argc != 3) {
    fprintf(stderr,"usage: matchcbc <region1> <region2>\n");
    exit(-1);
  }  

  cbc1 = ReadAp_FeatureVector1D(argv[1], &n1);
  cbc2 = ReadAp_FeatureVector1D(argv[2], &n2);

  wcbc1 = CreateFeatureVector1D(n1);
  for (i=0; i<n1; i++)
    wcbc1->X[i] = cbc1[i]->X[5];

  wcbc2 = CreateFeatureVector1D(n2);
  for (i=0; i<n2; i++)
    wcbc2->X[i] = cbc2[i]->X[5];

  tic = Tic();
  distance=IRM(cbc1, wcbc1, cbc2, wcbc2, CBCRegionDistance);
  toc = Toc();
  printf("%f\n",distance);

  for (i=0; i<n1; i++)
    DestroyFeatureVector1D(&cbc1[i]);
  free(cbc1);

  for (i=0; i<n2; i++)
    DestroyFeatureVector1D(&cbc2[i]);
  free(cbc2);

  DestroyFeatureVector1D(&wcbc1);
  DestroyFeatureVector1D(&wcbc2);

  return(0);
}
