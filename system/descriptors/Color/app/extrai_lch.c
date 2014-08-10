#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  CImage *cimg=NULL;
  Image *mask=NULL;
  Histogram *hist=NULL;
  //int i, n;

  if (argc != 4) {
    fprintf(stderr,"usage: extrai_bic <image> <mask> <outputfile>\n");
    exit(-1);
  }  

  cimg = ReadJPEGFile(argv[1]);
  mask = ReadImage(argv[2]);
  
  tic = Tic();
  hist = LCH(cimg, mask);  
  toc = Toc();
  printf("LCH extracted in %f milliseconds\n\n",CTime(tic, toc));
  WriteFileHistogram(hist,argv[3]);  
  
  DestroyCImage(&cimg);
  DestroyHistogram(&hist);
  DestroyImage(&mask);

  return(0);
}
