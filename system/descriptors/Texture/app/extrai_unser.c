#include "texturedescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Image *img=NULL;
  Image *mask=NULL;
  Histogram *hist=NULL;
  //int i, n;

  if (argc != 4) {
    fprintf(stderr,"usage: extrai_unser <image> <mask> <outputfile>\n");
    exit(-1);
  }  

  img = ReadImage(argv[1]);
  mask = ReadImage(argv[2]);

  tic = Tic();
  hist = Unser(img, mask);
  toc = Toc();
  printf("Unser extracted in %f milliseconds\n\n",CTime(tic, toc));
  
  WriteFileHistogram(hist, argv[3]);

  DestroyImage(&img);
  DestroyHistogram(&hist);
  DestroyImage(&mask);

  return(0);
}
