#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  CImage *cimg=NULL;
  Image *mask=NULL;
  Ap_FeatureVector1D *cbc=NULL;
  int i, n;

  if (argc != 4) {
    fprintf(stderr,"usage: generatecbc <image> <mask> <region>\n");
    exit(-1);
  }  

  cimg = ReadJPEGFile(argv[1]);
  mask = ReadImage(argv[2]);

  tic = Tic();
  cbc = CBC(cimg, mask, &n);
  toc = Toc();
  printf("CreateCbc in %f milliseconds\n\n",CTime(tic, toc));

  WriteAp_FeatureVector1D(cbc, n, argv[3]);

  for (i=0;i<n;i++)
    DestroyFeatureVector1D(&cbc[i]);
  free(cbc);

  DestroyCImage(&cimg);
  DestroyImage(&mask);

  return(0);
}
